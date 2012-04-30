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
  The distributed factorized engine implements the factorized
  consistency model.  */

#ifndef GRAPHLAB_DISTRIBUTED_FSCOPE_ENGINE_HPP
#define GRAPHLAB_DISTRIBUTED_FSCOPE_ENGINE_HPP

#include <deque>
#include <boost/bind.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/context/icontext.hpp>
#include <graphlab/context/context.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/update_functor/iupdate_functor.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>
#include <graphlab/aggregation/distributed_aggregator.hpp>

#include <graphlab/util/tracepoint.hpp>
#include <graphlab/util/memory_info.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/rpc/async_consensus.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {
  
  

  template<typename Graph, typename UpdateFunctor>
  class distributed_fscope_engine : public iengine<Graph, UpdateFunctor> {
 
  public:
    // Include parent types
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef distributed_fscope_engine<Graph, UpdateFunctor> engine_type;
    
    typedef typename iengine_base::graph_type graph_type;
    typedef typename graph_type::local_graph_type local_graph_type;
    typedef typename iengine_base::update_functor_type update_functor_type;
    
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef typename graph_type::local_edge_list_type local_edge_list_type;
    typedef typename graph_type::edge_type edge_type;
    typedef typename graph_type::lvid_type lvid_type;

    typedef ischeduler<local_graph_type, update_functor_type > ischeduler_type;
    typedef scheduler_factory<local_graph_type, update_functor_type> 
    scheduler_factory_type;

    typedef distributed_aggregator<distributed_fscope_engine> aggregator_type;

    
    typedef typename iengine_base::icontext_type  icontext_type;
    typedef context<distributed_fscope_engine>           context_type;

    
    
    consistency_model default_consistency;
    
    enum vertex_execution_state {
      NONE = 0,
      GATHERING,           // state on owner
      APPLYING,            // state on owner
      SCATTERING,          // state on owner
      WAIT_FOR_SCATTER_COMPLETE, // state on owner
      MIRROR_GATHERING,    // state on mirror
      MIRROR_SCATTERING,   // state on mirror
    }; // end of vertex execution state


    struct vertex_state {
      int32_t apply_scatter_count_down;    // used to count down the gathers and for the number of scatter responses
      simple_spinlock lock;
      bool hasnext;
      vertex_execution_state state; // current state of the vertex 
      update_functor_type current;  // What is currently being executed
                                    //  accumulated
      update_functor_type next;     // next is set if the vertex is being
                                    // executed, but for whatever reason
                                    // it got popped from the scheduler
                                    // again
      vertex_state(): apply_scatter_count_down(0), hasnext(false), 
                      state(NONE) { }
      std::ostream& operator<<(std::ostream& os) const {
        switch(state) {
        case NONE: { os << "NONE"; break; }
        case GATHERING: { os << "GATHERING: " << apply_scatter_count_down; break; }
        case APPLYING: { os << "APPLYING"; break; }
        case SCATTERING: { os << "SCATTERING"; break; }
        case MIRROR_GATHERING: { os << "MIRROR_GATHERING"; break; }
        case MIRROR_SCATTERING: { os << "MIRROR_SCATTERING"; break; }
        default: { }
        }
        return os;
      }
    }; // end of vertex_state
    
    
    struct thread_local_data {
      mutex lock;
      size_t npending;
      std::deque<vertex_id_type> pending_vertices;
      thread_local_data() : npending(0) { }       
      void add_task(vertex_id_type v) {
        lock.lock();
        ++npending;
        pending_vertices.push_back(v);
        lock.unlock();
      }
      void add_task_priority(vertex_id_type v) {
        lock.lock();
        ++npending;
        pending_vertices.push_front(v);
        lock.unlock();
      }
      bool get_task(std::deque<vertex_id_type> &v) {
        v = std::deque<vertex_id_type>();
        lock.lock();
        if (npending == 0) { lock.unlock(); return false; }
        npending = 0;
        v.swap(pending_vertices);
        lock.unlock();
        return true;
      }
    }; // end of thread local data

  private:
    dc_dist_object< distributed_fscope_engine > rmi;

    //! The local engine options
    graphlab_options opts; 
    //! the distributed graph
    graph_type& graph;
    //! local threads object
    thread_pool threads;
    //! the context manager
    std::vector<mutex> vlocks;
    //! The scheduler
    ischeduler_type* scheduler_ptr;
    //! The aggregator type
    aggregator_type aggregator;


    //! the vertex state
    std::vector<vertex_state> vstate;
    //! Engine state
    bool started;
    async_consensus consensus;
    atomic<size_t> threads_alive;
    std::vector<thread_local_data> thrlocal;
    
    atomic<uint64_t> blocked_issues; // issued tasks which cannot start
                                     // and have to be reinjected into the 
                                     // scheduler
    atomic<uint64_t> issued_tasks;
    atomic<uint64_t> completed_tasks;
    
    size_t max_pending_tasks;
    size_t recv_throttle_threshold;
    size_t send_throttle_threshold;
    size_t timed_termination;
    float engine_start_time;
    bool timeout;
    
    PERMANENT_DECLARE_DIST_EVENT_LOG(eventlog);
    DECLARE_TRACER(disteng_eval_sched_task);
    DECLARE_TRACER(disteng_init_gathering); 
    DECLARE_TRACER(disteng_init_scattering);
    DECLARE_TRACER(disteng_waiting_for_vstate_locks);
    DECLARE_TRACER(disteng_evalfac);
    DECLARE_TRACER(disteng_internal_task_queue);
    DECLARE_TRACER(disteng_scheduler_task_queue);

    
    bool try_to_quit(size_t threadid,
                     bool& has_internal_task,
                     std::deque<lvid_type>& internal_lvid,
                     bool& has_sched_task,
                     lvid_type& sched_lvid,
                     update_functor_type &task) {

      static size_t ctr = 0;
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, NO_WORK_EVENT, 1);
      
      if (lowres_time_seconds() - engine_start_time > timed_termination) {
        timeout = true;
      }
      if (!timeout &&
          issued_tasks.value != completed_tasks.value + blocked_issues.value) {
        ++ctr;
        if (ctr % 10 == 0) usleep(1);
        return false;
      }
      has_internal_task = false;
      has_sched_task = false;
      threads_alive.dec();
      consensus.begin_done_critical_section(threadid);
      
      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      if (thrlocal[threadid].get_task(internal_lvid)) {
        has_internal_task = true;
        consensus.cancel_critical_section(threadid);
        threads_alive.inc();
        END_TRACEPOINT(disteng_internal_task_queue);
        return false;
      }
      END_TRACEPOINT(disteng_internal_task_queue);
      
      sched_status::status_enum stat = 
        scheduler_ptr->get_next(threadid, sched_lvid, task);
      if (stat == sched_status::EMPTY) {
        bool ret = consensus.end_done_critical_section(threadid);
        threads_alive.inc();
        return ret;
      } else {
        consensus.cancel_critical_section(threadid);
        has_sched_task = true;
        threads_alive.inc();
        return false;
      }
    } // end of try to quit



    /// Locking code --------------------------------------------------------
    inline void read_lock_vertex(lvid_type vid) {
      ASSERT_LT(vid, vlocks.size()); 
      vlocks[vid].lock();
      // vlocks[vid].readlock();
    }

    inline bool try_read_lock_vertex(lvid_type vid) {
      ASSERT_LT(vid, vlocks.size()); 
      return vlocks[vid].try_lock();
      // return vlocks[vid].try_readlock();
    }

    inline void write_lock_vertex(lvid_type vid) {
      ASSERT_LT(vid, vlocks.size()); 
      vlocks[vid].lock();
      // vlocks[vid].writelock();
    }

    inline void release_lock_vertex(lvid_type vid) {
      ASSERT_LT(vid, vlocks.size()); 
      vlocks[vid].unlock();
    }

    inline void lock_single_edge(lvid_type center, lvid_type neighbor) {
      if(center < neighbor) { 
        write_lock_vertex(center); read_lock_vertex(neighbor); 
      } else { 
        read_lock_vertex(neighbor); write_lock_vertex(center); 
      }
    }

    inline void release_single_edge(lvid_type center, lvid_type neighbor) {
      release_lock_vertex(center); release_lock_vertex(neighbor);
    }

    inline void swap_single_edge(lvid_type center, lvid_type old_neighbor,
                                 lvid_type new_neighbor) {
      release_lock_vertex(old_neighbor);
      if( !try_read_lock_vertex(new_neighbor) ) {
        release_lock_vertex(center); lock_single_edge(center, new_neighbor);
      }      
    } // end of swap_single_edge




    inline void ASSERT_I_AM_OWNER(const lvid_type lvid) const {
      ASSERT_EQ(graph.l_get_vertex_record(lvid).owner, rmi.procid());
    }
    inline void ASSERT_I_AM_NOT_OWNER(const lvid_type lvid) const {
      ASSERT_NE(graph.l_get_vertex_record(lvid).owner, rmi.procid());
    }

    enum {
      SCHEDULE_EVENT = 0,
      UPDATE_EVENT = 1,
      WORK_ISSUED_EVENT = 2,
      INTERNAL_TASK_EVENT = 3,
      NO_WORK_EVENT = 4,
      SCHEDULE_FROM_REMOTE_EVENT = 5,
      USER_OP_EVENT = 6,
      ENGINE_START_EVENT = 7,
      ENGINE_STOP_EVENT = 8
    };
    
  public:
    distributed_fscope_engine(distributed_control& dc, graph_type& graph, 
                              size_t ncpus) : 
      rmi(dc, this), graph(graph), threads(ncpus),
      vlocks(graph.num_local_vertices()), 
      scheduler_ptr(NULL), 
      aggregator(dc, *this, graph),
      consensus(dc, ncpus), 
      max_pending_tasks(-1),
      recv_throttle_threshold(-1),
      send_throttle_threshold(-1),
      timed_termination(-1){
      rmi.barrier();
      aggregator.get_threads().resize(2);
#ifdef USE_EVENT_LOG
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 3000, 
                                          dist_event_log::RATE_BAR);
#else
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 3000, 
                                          dist_event_log::LOG_FILE);
#endif
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, SCHEDULE_EVENT, "Schedule");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, UPDATE_EVENT, "Updates");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, WORK_ISSUED_EVENT, "Work Issued");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, INTERNAL_TASK_EVENT, "Internal");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, NO_WORK_EVENT, "No Work");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, SCHEDULE_FROM_REMOTE_EVENT, 
                                    "Remote Schedule");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, USER_OP_EVENT, "User Ops");
      PERMANENT_ADD_IMMEDIATE_DIST_EVENT_TYPE(eventlog, ENGINE_START_EVENT, "Engine Start");
      PERMANENT_ADD_IMMEDIATE_DIST_EVENT_TYPE(eventlog, ENGINE_STOP_EVENT, "Engine Stop");
      
      INITIALIZE_TRACER(disteng_eval_sched_task, 
                        "distributed_fscope_engine: Evaluate Scheduled Task");
      INITIALIZE_TRACER(disteng_init_gathering,
                        "distributed_fscope_engine: Initialize Gather");
      INITIALIZE_TRACER(disteng_init_scattering,
                        "distributed_fscope_engine: Initialize Scattering");
      INITIALIZE_TRACER(disteng_waiting_for_vstate_locks,
                        "distributed_fscope_engine: vstate Lock Contention");
      INITIALIZE_TRACER(disteng_evalfac,
                        "distributed_fscope_engine: "
                        "Time in Factorized Update user code");
      INITIALIZE_TRACER(disteng_internal_task_queue,
                        "distributed_fscope_engine: Time in Internal Task Queue");
      INITIALIZE_TRACER(disteng_scheduler_task_queue,
                        "distributed_fscope_engine: Time in Scheduler Task Queue");
    }

    /**
     * Unique to the distributed engine. This must be called prior
     * to any schedule() call. 
     */
    void initialize() {      
      graph.finalize();
      if (rmi.procid() == 0) memory_info::print_usage("Before Engine Initialization");
      logstream(LOG_INFO) 
        << rmi.procid() << ": Initializing..." << std::endl;
      // currently this code wipes out any exisiting data structures
      if(scheduler_ptr != NULL) {  delete scheduler_ptr; scheduler_ptr = NULL; }
      
      // construct the scheduler
      scheduler_ptr = scheduler_factory_type::
        new_scheduler(opts.scheduler_type,
                      opts.scheduler_args,
                      graph.get_local_graph(),
                      threads.size());      
      // construct the context manager
      vlocks.resize(graph.num_local_vertices());
      vstate.resize(graph.num_local_vertices());
      
      thrlocal.resize(threads.size());
      if (rmi.procid() == 0) memory_info::print_usage("After Engine Initialization");
      rmi.barrier();
    }
    
    ~distributed_fscope_engine() {
      if(scheduler_ptr != NULL) {  delete scheduler_ptr; scheduler_ptr = NULL; }
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
           
    void schedule_from_remote(vertex_id_type vid,
                              const update_functor_type& update_functor) {
      if (timeout) return;
      const lvid_type local_vid = graph.local_vid(vid);
      BEGIN_TRACEPOINT(disteng_scheduler_task_queue);
      scheduler_ptr->schedule(local_vid, update_functor);
      END_TRACEPOINT(disteng_scheduler_task_queue);
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_FROM_REMOTE_EVENT, 1);
      consensus.cancel();
    }
    
    void schedule_local(vertex_id_type local_vid ,
                        const update_functor_type& update_functor) {
      if (timeout) return;
      if (started) {
        BEGIN_TRACEPOINT(disteng_scheduler_task_queue);
        scheduler_ptr->schedule_from_execution_thread(thread::thread_id(),
                                                      local_vid, update_functor);
        END_TRACEPOINT(disteng_scheduler_task_queue);
      }
      else {
        scheduler_ptr->schedule(local_vid, update_functor);
      }
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_EVENT, 1);
      consensus.cancel();
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
      std::vector<vertex_id_type> vtxs;
      vtxs.reserve(graph.get_local_graph().num_vertices());
      for(lvid_type lvid = 0; lvid < graph.get_local_graph().num_vertices(); 
          ++lvid) {
        if (graph.l_get_vertex_record(lvid).owner == rmi.procid()) 
          vtxs.push_back(lvid);        
      } 
      if(order == "shuffle") 
        graphlab::random::shuffle(vtxs.begin(), vtxs.end());
      foreach(lvid_type lvid, vtxs)
        scheduler_ptr->schedule(lvid, update_functor);    
      if (started) consensus.cancel();
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
      if(opts.engine_args.get_option("max_pending", max_pending_tasks)) {
        std::cout << "Max Pending: " << max_pending_tasks << std::endl;
      }
      if(opts.engine_args.get_option("timed_termination", timed_termination)) {
        std::cout << "Timed Termination (s): " << timed_termination << std::endl;
      }
      opts.engine_args.get_option("recv_throttle_threshold", recv_throttle_threshold);
      std::cout << "Throttle Threshold (calls): " << recv_throttle_threshold << std::endl;

      opts.engine_args.get_option("send_throttle_threshold", send_throttle_threshold);
      std::cout << "Send Throttle Threshold: " << send_throttle_threshold << std::endl;
    } 

    /** \brief get the current engine options. */
    const graphlab_options& get_options() { return opts; }
        
    void get_a_task(size_t threadid, 
                    bool& has_internal_task,
                    std::deque<lvid_type>& internal_lvid,
                    bool& has_sched_task,
                    lvid_type& sched_lvid,
                    update_functor_type &task) {
      // reset return task booleans
      has_internal_task = false;
      has_sched_task = false;
      if (lowres_time_seconds() - engine_start_time > timed_termination) {
        return;
      }
      
      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      if (thrlocal[threadid].get_task(internal_lvid)) {
        has_internal_task = true;
        END_TRACEPOINT(disteng_internal_task_queue);
        return;
      }
      END_TRACEPOINT(disteng_internal_task_queue);
      // Assert: no internal tasks
      
      // Determine the number of pending task (ensure less than
      // max_pending)
      const size_t pending_tasks = 
        issued_tasks - (completed_tasks + blocked_issues);        
      if( pending_tasks > max_pending_tasks ) { return; }
      // Assert the number of pending tasks is less than max pending
      // so proceed to start another task

      if (rmi.dc().recv_queue_length() > recv_throttle_threshold) {
        // std::cout << rmi.procid() << ": Throttle: " << rmi.dc().pending_queue_length() << std::endl;
        usleep(1000);
      }

      
      size_t sq = rmi.dc().send_queue_length();
      if (sq > send_throttle_threshold) {
        usleep(1);
        return;
      }
      // Get a fresh task from the scheduler
      sched_status::status_enum stat = 
        scheduler_ptr->get_next(threadid, sched_lvid, task);
      has_sched_task = stat != sched_status::EMPTY;
    } // end of get a task

    void gather_complete(const lvid_type lvid) {
      ASSERT_I_AM_OWNER(lvid);
      ASSERT_TRUE(vstate[lvid].state == GATHERING);
      ASSERT_GT(vstate[lvid].apply_scatter_count_down, 0);
      vstate[lvid].apply_scatter_count_down--;
      if (vstate[lvid].apply_scatter_count_down == 0) {
        vstate[lvid].state = APPLYING;
        add_internal_task(lvid);
      }
    } // end of gather complete
    
    void rpc_gather_complete(vertex_id_type vid, const update_functor_type& uf) {
      const vertex_id_type lvid = graph.local_vid(vid);
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate[lvid].lock.lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate[lvid].current.merge(uf);
      gather_complete(lvid);
      vstate[lvid].lock.unlock();
    } // end of rpc gather complete


    void scatter_complete(const lvid_type lvid) {
      ASSERT_I_AM_OWNER(lvid);
      int32_t ct = __sync_sub_and_fetch(&vstate[lvid].apply_scatter_count_down, 1);
      ASSERT_GE(ct, 0);
      
      if (ct == 0) {
        ASSERT_EQ(vstate[lvid].state, WAIT_FOR_SCATTER_COMPLETE);
        update_functor_type uf;
        bool has_new_task = false;
      
        vstate[lvid].lock.lock();
        if (vstate[lvid].hasnext) {
          uf = vstate[lvid].next;
          has_new_task = true;
          vstate[lvid].next = update_functor_type();
          vstate[lvid].current = update_functor_type();
          vstate[lvid].hasnext = false;
        }
        vstate[lvid].state = NONE;
        vstate[lvid].lock.unlock();

        if (has_new_task) eval_sched_task(lvid, uf);
      }
    } // end of gather complete

    void rpc_scatter_complete(vertex_id_type vid) {
      const vertex_id_type lvid = graph.local_vid(vid);
      scatter_complete(lvid);
    } // end of rpc gather complete

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
    
    void process_gather(lvid_type lvid) {
      // This function is called from within a vstate[lvid].lock;
      ASSERT_TRUE(vstate[lvid].state == GATHERING || 
                  vstate[lvid].state == MIRROR_GATHERING);      
      do_gather(lvid);
      const procid_t vowner = graph.l_get_vertex_record(lvid).owner;
      if (vowner == rmi.procid()) {
        gather_complete(lvid);
      } else {
        rmi.remote_call(vowner,
                        &engine_type::rpc_gather_complete,
                        graph.global_vid(lvid),
                        vstate[lvid].current);
        vstate[lvid].current = update_functor_type();
      }
    } // end of process gather

    
    void eval_internal_task(lvid_type lvid) {
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, INTERNAL_TASK_EVENT, 1);
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate[lvid].lock.lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      // move some calls out of the lock by setting flags;
      bool call_scatter_complete = false;
      bool remote_call_scatter_complete = false;
      
      switch(vstate[lvid].state) {
      case NONE: { logstream(LOG_FATAL) << "Empty Internal Task"; break; }// break;
      case GATHERING: { process_gather(lvid); break; }
      case MIRROR_GATHERING: { 
        do_init_gather(lvid); process_gather(lvid); 
        // mirrors are dumb upon complete of a gather, they go straight to NONE
        vstate[lvid].state = NONE;
        break; 
      }
      case APPLYING: { 
        do_apply(lvid);
        vstate[lvid].state = SCATTERING;
        if (vstate[lvid].current.scatter_edges() == graphlab::NO_EDGES) {
          vstate[lvid].apply_scatter_count_down = 1;
        }
        else {
          vstate[lvid].apply_scatter_count_down = graph.l_get_vertex_record(lvid).num_mirrors() + 1;
        }
        master_broadcast_scattering(lvid,
                                    vstate[lvid].current,
                                    graph.get_local_graph().vertex_data(lvid));
        // fall through to scattering
      }
      case SCATTERING: {
        do_scatter(lvid);
        vstate[lvid].state = WAIT_FOR_SCATTER_COMPLETE;
        call_scatter_complete = true;
        completed_tasks.inc();
        break;
      }
      case MIRROR_SCATTERING: {
        do_scatter(lvid);
        vstate[lvid].current = update_functor_type();
        vstate[lvid].state = NONE;
        // remote call scatter complete to the owner outside of the locked region
        remote_call_scatter_complete = true;
        break;
      }
      case WAIT_FOR_SCATTER_COMPLETE: ASSERT_TRUE(false);
      } // end of switch statement
      vstate[lvid].lock.unlock();
      // complete some calls outside of the locked region
      if (call_scatter_complete) {
        scatter_complete(lvid);
      }
      else if (remote_call_scatter_complete) {
        rmi.remote_call(graph.l_get_vertex_record(lvid).owner,
                        &distributed_fscope_engine::rpc_scatter_complete,
                        graph.global_vid(lvid));
      }
    } // end of eval internal task


    void add_internal_task(lvid_type lvid) {
      if (timeout) return;
      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      size_t i = lvid % threads.size();

      thrlocal[i].add_task(lvid);
      consensus.cancel_one(i);
      END_TRACEPOINT(disteng_internal_task_queue);
    }


    // If I receive the call I am a mirror of this vid
    void rpc_begin_gathering(vertex_id_type sched_vid, 
                             const update_functor_type& task) {
      ASSERT_NE(graph.get_vertex_record(sched_vid).owner, rmi.procid());
      // immediately begin issuing the lock requests
      vertex_id_type sched_lvid = graph.local_vid(sched_vid);
      // set the vertex state
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate[sched_lvid].lock.lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      ASSERT_EQ(vstate[sched_lvid].state, NONE);
      vstate[sched_lvid].state = MIRROR_GATHERING;
      vstate[sched_lvid].current = task;
      vstate[sched_lvid].lock.unlock();
      add_internal_task(sched_lvid);
    } // end of rpc begin gathering


    /**
     * Task was added to the vstate. Now to begin scheduling the gathers
     */
    void master_broadcast_gathering(lvid_type sched_lvid,
                                    const update_functor_type& task) {
      // check to see if there are no edges to gather on.  If this is
      // the case we can skip the broadcast 

      // Assert that we do want to proceed with a gather
      BEGIN_TRACEPOINT(disteng_init_gathering);
      ASSERT_I_AM_OWNER(sched_lvid);
      const vertex_id_type sched_vid = graph.global_vid(sched_lvid);
      const typename graph_type::vertex_record& vrec = 
        graph.l_get_vertex_record(sched_lvid);

      if (task.gather_edges() == graphlab::NO_EDGES) {
        vstate[sched_lvid].apply_scatter_count_down = 1;
        gather_complete(sched_lvid);
      }
      else {
        const unsigned char prevkey = 
        rmi.dc().set_sequentialization_key(sched_vid % 254 + 1);
        rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                      &engine_type::rpc_begin_gathering, sched_vid, task);
        rmi.dc().set_sequentialization_key(prevkey);
        add_internal_task(sched_lvid);
      }

      END_TRACEPOINT(disteng_init_gathering);
    }


    void rpc_begin_scattering(vertex_id_type vid, update_functor_type task,
                              const vertex_data_type &central_vdata) {
      const vertex_id_type lvid = graph.local_vid(vid);
      vstate[lvid].lock.lock();
      ASSERT_I_AM_NOT_OWNER(lvid);
      graph.get_local_graph().vertex_data(lvid) = central_vdata;
      if (task.scatter_edges() != NO_EDGES) {
        vstate[lvid].state = MIRROR_SCATTERING;
        vstate[lvid].current = task;
      }
      vstate[lvid].lock.unlock();
      add_internal_task(lvid);

    }


    void rpc_scatter_data_only(vertex_id_type vid, 
                              const vertex_data_type &central_vdata) {
      const vertex_id_type lvid = graph.local_vid(vid);
      vstate[lvid].lock.lock();
      ASSERT_I_AM_NOT_OWNER(lvid);
      graph.get_local_graph().vertex_data(lvid) = central_vdata;
      vstate[lvid].lock.unlock();
    }
    /**
     * Task was added to the vstate. Now to begin scheduling the gathers
     */
    void master_broadcast_scattering(lvid_type sched_lvid,
                                     const update_functor_type& task,
                                     const vertex_data_type &central_vdata) {
      BEGIN_TRACEPOINT(disteng_init_scattering);
      ASSERT_I_AM_OWNER(sched_lvid);
      const vertex_id_type sched_vid = graph.global_vid(sched_lvid);
      const typename graph_type::vertex_record& vrec = 
        graph.l_get_vertex_record(sched_lvid);
      const unsigned char prevkey = 
        rmi.dc().set_sequentialization_key(sched_vid % 254 + 1);
      if (task.scatter_edges() != graphlab::NO_EDGES) {
        rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                        &engine_type::rpc_begin_scattering,
                        sched_vid, task, central_vdata);
      }
      else {
        rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                        &engine_type::rpc_scatter_data_only,
                        sched_vid, central_vdata);
      }
      rmi.dc().set_sequentialization_key(prevkey);
      END_TRACEPOINT(disteng_init_scattering);
    } // end of master_broadcast_scattering



    void eval_sched_task(const lvid_type sched_lvid, 
                         const update_functor_type& task) {

      BEGIN_TRACEPOINT(disteng_eval_sched_task);
      // If I am not the owner just forward the task to the other
      // scheduler and return
      const typename graph_type::vertex_record& rec = graph.l_get_vertex_record(sched_lvid);
      const procid_t owner = rec.owner;
      if (owner != rmi.procid()) {
        const vertex_id_type vid = rec.gvid;
        rmi.remote_call(owner, &engine_type::schedule_from_remote, vid, task);
        return;
      }
      ASSERT_I_AM_OWNER(sched_lvid);

      // this is in local VIDs
      issued_tasks.inc();
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate[sched_lvid].lock.lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      bool initiate_gathering = false;
      if (vstate[sched_lvid].state == NONE) {
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, UPDATE_EVENT, 1);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, WORK_ISSUED_EVENT, 
                                        rec.num_in_edges + rec.num_out_edges);
        // we start gather right here.
        // set up the state
        vstate[sched_lvid].state = GATHERING;
        vstate[sched_lvid].current = task;
        vstate[sched_lvid].apply_scatter_count_down =   
          graph.l_get_vertex_record(sched_lvid).num_mirrors() + 1;
        initiate_gathering = true;
        // we are going to broadcast after unlock
      } else {
        blocked_issues.inc();
        if (vstate[sched_lvid].hasnext) {
          vstate[sched_lvid].next += task;
        } else {
          vstate[sched_lvid].hasnext = true;
          vstate[sched_lvid].next = task;
        }
      }
      END_TRACEPOINT(disteng_eval_sched_task);
      if (initiate_gathering) do_init_gather(sched_lvid);
      vstate[sched_lvid].lock.unlock();

      // begin gathering
      if (initiate_gathering) {
        // Broadcast the uninitialized task to all mirrors
        master_broadcast_gathering(sched_lvid, task);
      }
    } // eval sched task
    
    void thread_start(size_t threadid) {
      bool has_internal_task = false;
      bool has_sched_task = false;
      std::deque<vertex_id_type> internal_lvid;
      vertex_id_type sched_lvid;
      update_functor_type task;
      while(1) {
        aggregator.evaluate_queue();
        get_a_task(threadid, 
                   has_internal_task, internal_lvid,
                   has_sched_task, sched_lvid, task);
        // if we managed to get a task..
        if (has_internal_task) {
          while(!internal_lvid.empty()) {
            eval_internal_task(internal_lvid.front());
            internal_lvid.pop_front();
          }
        } else if (has_sched_task) { 
          eval_sched_task(sched_lvid, task);
        } else if (!try_to_quit(threadid,
                                has_internal_task, internal_lvid,
                                has_sched_task, sched_lvid, task)) {
          if (has_internal_task) {
            while(!internal_lvid.empty()) {
              eval_internal_task(internal_lvid.front());
              internal_lvid.pop_front();
            }
          } else if (has_sched_task) {
            eval_sched_task(sched_lvid, task);
          }
        } else { break; }
      }
    } // end of thread start
    
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
      aggregator.initialize_queue();
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

      size_t allocatedmem = memory_info::allocated_bytes();
      rmi.all_reduce(allocatedmem);
 
      engine_start_time = lowres_time_seconds();
      timeout = false;
      rmi.dc().flush_counters();
      if (rmi.procid() == 0) {
        logstream(LOG_INFO) << "Total Allocated Bytes: " << allocatedmem << std::endl;
        PERMANENT_IMMEDIATE_DIST_EVENT(eventlog, ENGINE_START_EVENT);
      }

      for (size_t i = 0; i < threads.size(); ++i) {
        const boost::function<void (void)> run_function = 
          boost::bind(&engine_type::thread_start, this, i);
        threads.launch(run_function);
      }
      // TODO: This should be in a while loop to catch all exceptions
      join_threads(threads);
      join_threads(aggregator.get_threads());
      
      rmi.dc().flush_counters();
      if (rmi.procid() == 0) {
        PERMANENT_IMMEDIATE_DIST_EVENT(eventlog, ENGINE_STOP_EVENT);
      }
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
    } // end of start



    void join_threads(thread_pool& threads) {
      threads.join();          
      // // Join all the threads (looping while failed)
      // bool join_successful = false;   
      // while(!join_successful) {
      //   try { 
      //     // Join blocks until all exection threads return.  However if
      //     // a thread terminates early due to exception in user code
      //     // (e.g., update function throws an exception) this join may
      //     // fail before all threads are joined.  If this happens we
      //     // catch the exception and then try to kill the rest of the
      //     // threads and proceed to join again.
      //     threads.join();
      //     join_successful = true; 
      //   } catch (const char* error_str) {
      //     logstream(LOG_ERROR) 
      //       << "Exception in execution thread caught: " << error_str 
      //       << std::endl;     
      //     // // killall the running threads
      //     exception_message = error_str;
      //     exec_status = execution_status::EXCEPTION;
      //   } 
      // }
    } // end of join_threads

  


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
                        size_t interval) {
      logstream(LOG_DEBUG) << "Add Aggregator : " << key << std::endl;
      aggregator.add_aggregator(key, zero, interval);
    }


    //! Performs a sync immediately.
    void aggregate_now(const std::string& key) {
      logstream(LOG_DEBUG) << "Aggregate Now : " << key << std::endl;
      aggregator.aggregate_now(key);
    }
  }; // end of class
}; // namespace

#include <graphlab/macros_undef.hpp>

#endif // GRAPHLAB_DISTRIBUTED_FSCOPE_ENGINE_HPP

