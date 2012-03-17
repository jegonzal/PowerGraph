/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */




/*
  fscope engine in the fully synchronous setting  */

#ifndef GRAPHLAB_DISTRIBUTED_SYNCHRONOUS_ENGINE_HPP
#define GRAPHLAB_DISTRIBUTED_SYNCHRONOUS_ENGINE_HPP

#include <deque>
#include <boost/bind.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/context/icontext.hpp>
#include <graphlab/context/context.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/update_functor/iupdate_functor.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>


#include <graphlab/util/tracepoint.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/rpc/async_consensus.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {
  
  

  template<typename Graph, typename UpdateFunctor>
  class distributed_synchronous_engine : public iengine<Graph, UpdateFunctor> {
 
  public:
    // Include parent types
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef distributed_synchronous_engine<Graph, UpdateFunctor> engine_type;
    
    typedef typename iengine_base::graph_type graph_type;
    typedef typename graph_type::local_graph_type local_graph_type;
    typedef typename iengine_base::update_functor_type update_functor_type;
    
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef typename graph_type::local_edge_list_type local_edge_list_type;
    typedef typename graph_type::edge_type edge_type;
    typedef typename graph_type::lvid_type lvid_type;

    typedef ischeduler<local_graph_type, update_functor_type > ischeduler_type;
    
    typedef typename iengine_base::icontext_type  icontext_type;
    typedef context<distributed_fscope_engine>           context_type;

    
    
    consistency_model default_consistency;
    
    
  private:
    dc_dist_object<distributed_synchronous_engine<Graph, UpdateFunctor> > rmi;

    //! The local engine options
    graphlab_options opts; 
    //! the distributed graph
    graph_type& graph;
    //! local threads object
    thread_pool threads;

    //! Engine state
    bool started;
    std::vector<thread_local_data> thrlocal;

    std::vector<spinlock> vlocks;
    std::vector<UpdateFunctor> workingset;
    vertex_functor_set scheduleset;
    dense_bitset hastask;
    
    atomic<uint64_t> issued_tasks;
    atomic<uint64_t> completed_tasks;
    
    PERMANENT_DECLARE_DIST_EVENT_LOG(eventlog);
    DECLARE_TRACER(disteng_internal_task_queue);
    DECLARE_TRACER(disteng_scheduler_task_queue);


  public:
    distributed_synchronous_engine(distributed_control &dc, graph_type& graph,
                              size_t ncpus) : 
      rmi(dc, this), graph(graph), threads(ncpus) {
      rmi.barrier();
    }

    /**
     * Unique to the distributed engine. This must be called prior
     * to any schedule() call. 
     */
    void initialize() {
      graph.finalize();
      vlocks.resize(graph.num_local_vertices());
      workingset.resize(graph.num_local_vertices());
      hastask.resize(graph.num_local_vertices());
      scheduleset.resize(graph.num_local_vertices());
      
      logstream(LOG_INFO) 
        << rmi.procid() << ": Initializing..." << std::endl;
      rmi.barrier();
    }
    
    ~distributed_synchronous_engine() {
    }
    
    /**
     * \brief Force engine to terminate immediately.
     *
     * This function is used to stop the engine execution by forcing
     * immediate termination.  Any existing update tasks will finish
     * but no new update tasks will be started and the call to start()
     * will return.
     */
    void stop() { }

    
    /**
     * \brief Describe the reason for termination.
     *
     * Return the reason for the last termination.
     */
    execution_status::status_enum last_exec_status() const { 
      return execution_status::UNSET ; }
   
    /**
     * \brief Get the number of updates executed by the engine.
     *
     * This function returns the numbe of updates executed by the last
     * run of this engine.
     * 
     * \return the total number of updates
     */
    size_t last_update_count() const { 
      return completed_tasks.value;
    }
           
    void schedule_local(vertex_id_type local_vid ,
                        const update_functor_type& update_functor) {
      scheduleset.add(local_vid, update_functor);
      hastask.set(local_vid);
    }

    /**
     * \brief Adds an update task with a particular priority.
     * This function is forwarded to the scheduler.
     */
    void schedule(vertex_id_type vid,
                  const update_functor_type& update_functor) {
      schedule_local(graph.local_vid(vid), update_functor);
    } // end of schedule


    /**
     * \brief Creates a collection of tasks on all the vertices in the
     * graph, with the same update function and priority This function
     * is forwarded to the scheduler. Must be called by all machines
     * simultaneously
     */
    void schedule_all(const update_functor_type& update_functor,
                      const std::string& order = "shuffle") {
      for(lvid_type lvid = 0; lvid < graph.get_local_graph().num_vertices(); ++lvid) {
        if (graph.l_get_vertex_record(lvid).owner == rmi.procid()) 
          workingset[lvid] = update_functor;
      } 
      rmi.barrier();
    } // end of schedule all


    /**state
     * Schedule an update on all the neighbors of a particular vertex
     */
    void schedule_in_neighbors(const vertex_id_type& vertex, 
                               const update_functor_type& update_fun) {
      assert(false); //TODO: IMPLEMENT
    } // end of schedule in neighbors

    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    void schedule_out_neighbors(const vertex_id_type& vertex, 
                                const update_functor_type& update_fun) {
      assert(false); //TODO: IMPLEMENT
    } // end of schedule out neighbors
                                                  
    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    virtual void schedule_neighbors(const vertex_id_type& vertex, 
                                    const update_functor_type& update_fun) {
      assert(false); //TODO: IMPLEMENT
    } // end of schedule neighbors


    // /**
    //  * \brief associate a termination function with this engine.
    //  *
    //  * An engine can typically have many termination functions
    //  * associated with it. A termination function is a function which
    //  * takes a constant reference to the shared data and returns a
    //  * boolean which is true if the engine should terminate execution.
    //  *
    //  * A termination function has the following type:
    //  * \code
    //  * bool term_fun(const ishared_data_type* shared_data)
    //  * \endcode
    //  */
    // void add_termination_condition(termination_function_type term) { 
    //   rmi.barrier();
    // }
    // //!  remove all associated termination functions
    // void clear_termination_conditions() { 
    //   rmi.barrier();
    // };
    
    /**
     *  \brief The timeout is the total
     *  ammount of time in seconds that the engine may run before
     *  exeuction is automatically terminated.
     */
    void set_timeout(size_t timeout_secs) { assert(false); };

    /**
     * elapsed time in milliseconds since start was called
     */
    size_t elapsed_time() const { assert(false); return 0; }
    
    /**
     * \brief set a limit on the number of tasks that may be executed.
     * 
     * By once the engine has achived the max_task parameter execution
     * will be terminated. If max_tasks is set to zero then the
     * task_budget is ignored.  If max_tasks is greater than zero than
     * the value of max tasks is used.  Note that if max_task is
     * nonzero the engine encurs the cost of an additional atomic
     * operation in the main loop potentially reducing the overall
     * parallel performance.
     */
    void set_task_budget(size_t max_tasks) { }


    /** \brief Update the engine options.  */
    void set_options(const graphlab_options& new_opts) {
      opts = new_opts;
    } 

    /** \brief get the current engine options. */
    const graphlab_options& get_options() { return opts; }

    void do_apply(lvid_type lvid) { 
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, 1);
      BEGIN_TRACEPOINT(disteng_evalfac);
      const vertex_id_type vid = graph.global_vid(lvid);
      update_functor_type& ufun = vstate[lvid].current;
      context_type context(this, &graph, vid, VERTEX_CONSISTENCY);
      write_lock_vertex(lvid);
      ufun.apply(context);
      release_lock_vertex(lvid);
      END_TRACEPOINT(disteng_evalfac);
    }
 

    void do_init_gather(lvid_type lvid) {
      update_functor_type& ufun = vstate[lvid].current;
      const vertex_id_type vid = graph.global_vid(lvid);
      context_type context(this, &graph, vid, ufun.gather_consistency());
      ufun.init_gather(context);
    }
    
    void do_gather(lvid_type lvid) { // Do gather
      BEGIN_TRACEPOINT(disteng_evalfac);
      update_functor_type& ufun = vstate[lvid].current;
      const vertex_id_type vid = graph.global_vid(lvid);
      context_type context(this, &graph, vid, ufun.gather_consistency());
      if(ufun.gather_edges() == graphlab::IN_EDGES || 
         ufun.gather_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_in_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          const lvid_type lneighbor = edges[i].l_source();
          if(i == 0) lock_single_edge(lvid, lneighbor);
          else swap_single_edge(lvid, edges[i-1].l_source(), lneighbor);
          ufun.gather(context, edges[i]);
          if(i+1 == edges.size()) release_single_edge(lvid, lneighbor);
        }
      }
      if(ufun.gather_edges() == graphlab::OUT_EDGES ||
         ufun.gather_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_out_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          const lvid_type lneighbor = edges[i].l_target();
          if(i == 0) lock_single_edge(lvid, lneighbor);
          else swap_single_edge(lvid, edges[i-1].l_target(), lneighbor);
          ufun.gather(context, edges[i]);
          if(i+1 == edges.size()) release_single_edge(lvid, lneighbor);
        }
      }
      END_TRACEPOINT(disteng_evalfac);
    } // end of do_gather
    
    void do_scatter(lvid_type lvid) {
      BEGIN_TRACEPOINT(disteng_evalfac);
      const vertex_id_type vid = graph.global_vid(lvid);
      update_functor_type& ufun = vstate[lvid].current;
      context_type context(this, &graph, vid, ufun.scatter_consistency());
      if(ufun.scatter_edges() == graphlab::IN_EDGES || 
         ufun.scatter_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_in_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          const lvid_type lneighbor = edges[i].l_source();
          if(i == 0) lock_single_edge(lvid, lneighbor);
          else swap_single_edge(lvid, edges[i-1].l_source(), lneighbor);
          ufun.scatter(context, edges[i]);
          if(i+1 == edges.size()) release_single_edge(lvid, lneighbor);
        }
      }
      if(ufun.scatter_edges() == graphlab::OUT_EDGES ||
         ufun.scatter_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_out_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          const lvid_type lneighbor = edges[i].l_target();
          if(i == 0) lock_single_edge(lvid, lneighbor);
          else swap_single_edge(lvid, edges[i-1].l_target(), lneighbor);
          ufun.scatter(context, edges[i]);
          if(i+1 == edges.size()) release_single_edge(lvid, lneighbor);
        }
      }
      END_TRACEPOINT(disteng_evalfac);
    } // end of do scatter
    
    void transmit_schedule_and_shift_to_active() {
      
    }


    void parallel_gather() {
      
    }
    
    /**
     * \brief Start the engine execution.
     *
     * This \b blocking function starts the engine and does not
     * return until either one of the termination conditions evaluate
     * true or the scheduler has no tasks remaining.
     */
    void start() {
      logstream(LOG_INFO) 
        << "Spawning " << threads.size() << " threads" << std::endl;
      ASSERT_TRUE(scheduler_ptr != NULL);
      // start the scheduler
      scheduler_ptr->start();
      started = true;
      threads_alive.value = threads.size();
      default_consistency = 
        string_to_consistency_model(opts.get_scope_type());

      ASSERT_TRUE(default_consistency == EDGE_CONSISTENCY || 
                  default_consistency == VERTEX_CONSISTENCY);
      if (default_consistency == EDGE_CONSISTENCY) {
        std::cout << "Consistency Model is: Edge" << std::endl;
      } else {
        std::cout << "Consistency Model is: Vertex" << std::endl;      
      }
      rmi.barrier();
      for (size_t i = 0; i < threads.size(); ++i) {
        const boost::function<void (void)> run_function = 
          boost::bind(&engine_type::thread_start, this, i);
        threads.launch(run_function);
      }
      threads.join();
      size_t ctasks = completed_tasks.value;
      rmi.all_reduce(ctasks);
      completed_tasks.value = ctasks;
      
      ctasks = issued_tasks.value;
      rmi.all_reduce(ctasks);
      issued_tasks.value = ctasks;
      
      ctasks = blocked_issues.value;
      rmi.all_reduce(ctasks);
      blocked_issues.value = ctasks;
      if (rmi.procid() == 0) {
        std::cout << "Completed Tasks: " << completed_tasks.value << std::endl;
        std::cout << "Issued Tasks: " << issued_tasks.value << std::endl;
        std::cout << "Blocked Issues: " << blocked_issues.value << std::endl;
      }
      // test if all schedulers are empty
      for (size_t i = 0; i < threads.size(); ++i) {
        if (thrlocal[i].npending) {
          logstream(LOG_WARNING) << rmi.procid() << ": Found " 
                                 << thrlocal[i].npending 
                                 << " Internal Tasks on CPU " << i << std::endl;
        }
        vertex_id_type sched_lvid;
        update_functor_type task;
        sched_status::status_enum stat = 
          scheduler_ptr->get_next(i, sched_lvid, task);
        if (stat == sched_status::NEW_TASK) {
          logstream(LOG_WARNING) << rmi.procid() << ": Found Scheduler Tasks" 
                                 << " on CPU " << i << std::endl;
        }
      }
      /*for (size_t i = 0;i < vstate.size(); ++i) {
        if(vstate[i].state != NONE) {
        std::cout << "Vertex: " << i << ": " << vstate[i].state << " " << (int)(cmlocks->philosopherset[i].state) << " " << cmlocks->philosopherset[i].num_edges << " " << cmlocks->philosopherset[i].forks_acquired << "\n";
          
        foreach(typename local_graph_type::edge_type edge, cmlocks->graph.in_edges(i)) {
        std::cout << (int)(cmlocks->forkset[cmlocks->graph.edge_id(edge)]) << " ";
        }
        std::cout << "\n";
        foreach(typename local_graph_type::edge_type edge, cmlocks->graph.out_edges(i)) {
        std::cout << (int)(cmlocks->forkset[cmlocks->graph.edge_id(edge)]) << " ";
        }
        std::cout << "\n";
        getchar();
        }
        }*/
    }
  


    /////////////////////////// Global Variabes ////////////////////////////
    //! Get the global data and lock
    void get_global(const std::string& key,
                    graphlab::any_vector*& ret_values_ptr,
                    bool& ret_is_const) { }
  
    //! Get the global data and lock
    void acquire_global_lock(const std::string& key,
                             size_t index = 0) { }
    //! Release the global data lock
    void release_global_lock(const std::string& key,
                             size_t index = 0) { }
  
  public:

    /**
     * Define a global mutable variable (or vector of variables).
     *
     * \param key the name of the variable (vector)
     * \param value the initial value for the variable (vector)
     * \param size the initial size of the global vector (default = 1)
     *
     */
    template< typename T >
    void add_global(const std::string& key, const T& value,
                    size_t size = 1) {
      logstream(LOG_DEBUG) << "Add Global: " << key << std::endl;
    }

    /**
     * Define a global constant.
     */
    template< typename T >
    void add_global_const(const std::string& key, const T& value,
                          size_t size = 1) {
      logstream(LOG_DEBUG) << "Add Global Const: " << key << std::endl;
    }


    //! Change the value of a global entry
    template< typename T >
    void set_global(const std::string& key, const T& value,
                    size_t index = 0) {
      logstream(LOG_DEBUG) << "Set Global: " << key << std::endl;
    }

    //! Get a copy of the value of a global entry
    template< typename T >
    T get_global(const std::string& key, size_t index = 0) {
      logstream(LOG_DEBUG) << "Get Global : " << key << std::endl;
      return T();
    }

    //! \brief Registers an aggregator with the engine
    template<typename Aggregate>
    void add_aggregator(const std::string& key,
                        const Aggregate& zero,
                        size_t interval,
                        bool use_barrier = false,
                        vertex_id_type begin_vid = 0,
                        vertex_id_type end_vid =
                        std::numeric_limits<vertex_id_type>::max()) {
      logstream(LOG_DEBUG) << "Add Aggregator : " << key << std::endl;
    }


    //! Performs a sync immediately.
    void aggregate_now(const std::string& key) {
      logstream(LOG_DEBUG) << "Aggregate Now : " << key << std::endl;
    }
  }; // end of class
}; // namespace

#include <graphlab/macros_undef.hpp>

#endif // GRAPHLAB_DISTRIBUTED_FSCOPE_ENGINE_HPP

