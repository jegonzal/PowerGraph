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




/* \file iengine.hpp
   \brief The file containing the iengine description
   
   This file contains the description of the engine interface.  All
   graphlab engines (single_threaded, multi_threaded, distributed, ...)
   should satisfy the functionality described below.
*/

#ifndef GRAPHLAB_DISTRIBUTED_ENGINE_HPP
#define GRAPHLAB_DISTRIBUTED_ENGINE_HPP

#include <queue>
#include <boost/bind.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/context/icontext.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/update_functor/iupdate_functor.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>
#include <graphlab/util/chandy_misra.hpp>

#include <graphlab/rpc/async_consensus.hpp>

#include <graphlab/macros_def.hpp>

#define VERTEX_LOCK_COUNT 131072

namespace graphlab {
  
  /**
   * Used to lie to the scheduler about the graph type
   */
  template <typename Graph, typename UpdateFunctor>
  class pseudo_engine {
   public:
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef typename iengine_base::update_functor_type update_functor_type;

    typedef typename Graph::local_graph_type graph_type;
    
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef typename graph_type::edge_type edge_type;
    
    typedef ischeduler<pseudo_engine<Graph, UpdateFunctor> > ischeduler_type;
  };
  

  template<typename Graph, typename UpdateFunctor>
  class distributed_engine: public iengine<Graph, UpdateFunctor> {
    
  public:
    // Include parent types
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef distributed_engine<Graph, UpdateFunctor> engine_type;
    
    typedef typename iengine_base::graph_type graph_type;
    typedef typename iengine_base::update_functor_type update_functor_type;
    
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef typename graph_type::edge_type edge_type;

    typedef ischeduler<pseudo_engine<Graph, UpdateFunctor> > ischeduler_type;
    
    typedef typename iengine_base::icontext_type  icontext_type;
    //typedef context<distributed_engine>         context_type;
    //typedef context_manager<distributed_engine> context_manager_type;
   
    
    typedef typename iengine_base::termination_function_type termination_function_type;

    enum vertex_execution_state {
      NONE = 0,
      GATHERING,   // state on owner
      APPLYING,    // state on owner
      SCATTERING,  // state on owner
      MIRROR_GATHERING, // state on mirror
      MIRROR_SCATTERING, // state on mirror
    };
    struct vertex_state {
      uint32_t apply_count_down; // used to count down the gathers
      bool hasnext;
      vertex_execution_state state; // current state of the vertex 
      update_functor_type current; // What is currently being executed / accumulated
      update_functor_type next; // next is set if the vertex is being executed, 
                                // but for whatever reason it got popped from the
                                // scheduler again
      vertex_state(): apply_count_down(0), hasnext(false), state(NONE) { }
      
      std::ostream& operator<<(std::ostream& os) const {
        switch(state) {
          case NONE:
            os << "NONE";
            break;
          case GATHERING:
            os << "GATHERING: " << apply_count_down;
            break;
          case APPLYING:
            os << "APPLYING";
            break;
          case SCATTERING:
            os << "SCATTERING";
            break;
          case MIRROR_GATHERING:
            os << "MIRROR_GATHERING";
            break;
          case MIRROR_SCATTERING:
            os << "MIRROR_SCATTERING";
            break;
        }
        return os;
      }
    };
    
    
    struct thread_local_data {
      mutex lock;
      size_t npending;
      std::queue<vertex_id_type> pending_vertices;
      
      thread_local_data():npending(0) { } 
      
      void add_task(vertex_id_type v) {
        lock.lock();
        ++npending;
        pending_vertices.push(v);
        lock.unlock();
      }
      bool get_task(vertex_id_type &v) {
        lock.lock();
        if (npending == 0) {
          lock.unlock();
          return false;
        }
        --npending;
        v = pending_vertices.front();
        pending_vertices.pop();
        lock.unlock();
        return true;
      }
    };
  private:
    dc_dist_object<distributed_engine<Graph, UpdateFunctor> > rmi;

    //! The local engine options
    graphlab_options opts; 

    graph_type& graph;
    chandy_misra<typename graph_type::local_graph_type>* cmlocks;
    
    thread_group thrgroup;
    
    //! The scheduler
    ischeduler_type* scheduler_ptr;
    
    std::vector<vertex_state> vstate;
    std::vector<mutex> vstate_locks;
    
    size_t ncpus;
    bool started;
    async_consensus* consensus;
    atomic<size_t> threads_alive;
    std::vector<thread_local_data> thrlocal;
    
  bool try_to_quit(size_t threadid, 
                   bool& has_internal_task,
                   vertex_id_type& internal_lvid,
                   bool& has_sched_task,
                   vertex_id_type& sched_lvid,
                   update_functor_type &task) {
    logstream(LOG_DEBUG) << rmi.procid() << ": " << "Termination Attempt" << std::endl;
    has_internal_task = false;
    has_sched_task = false;
    threads_alive.dec();
    consensus->begin_done_critical_section();
    if (thrlocal[threadid].get_task(internal_lvid)) {
      logstream(LOG_DEBUG) << rmi.procid() << ": " << "\tCancelled by Internal Task" << std::endl;
      has_internal_task = true;
      consensus->end_done_critical_section(false);
      return false;
    }
    sched_status::status_enum stat = scheduler_ptr->get_next(threadid, sched_lvid, task);
    if (stat == sched_status::EMPTY) {
      logstream(LOG_DEBUG) << rmi.procid() << ": " << "\tTermination Success" << std::endl;
      bool ret = consensus->end_done_critical_section(true);
      threads_alive.inc();
      return ret;
    } else {
      logstream(LOG_DEBUG) << rmi.procid() << ": " << "\tCancelled by Scheduler Task" << std::endl;
      consensus->end_done_critical_section(false);
      has_sched_task = true;
      threads_alive.inc();
      return false;
    }
  } // end of try to qui

    
  public:
    distributed_engine(distributed_control &dc, graph_type& graph, size_t ncpus): 
                                                        rmi(dc, this), graph(graph),
                                                        scheduler_ptr(NULL),
                                                        ncpus(ncpus){
      rmi.barrier();
    }

    /**
     * Unique to the distributed engine. This must be called prior
     * to any schedule() call. 
     */
    void initialize() {
      logstream(LOG_INFO) << rmi.procid() << ": Initializing..." << std::endl;
      scheduler_ptr = scheduler_factory<pseudo_engine<Graph, UpdateFunctor> >::
                          new_scheduler(opts.scheduler_type,
                                      opts.scheduler_args,
                                      graph.get_local_graph(),
                                      ncpus);
      vstate.resize(graph.num_local_vertices());
      vstate_locks.resize(VERTEX_LOCK_COUNT);
      consensus = new async_consensus(rmi.dc(), ncpus);
      cmlocks = new chandy_misra<typename graph_type::local_graph_type>(graph.get_local_graph());
      thrlocal.resize(ncpus);
      rmi.barrier();
    }
    
    ~distributed_engine() {
      delete scheduler_ptr;
      delete consensus;
      delete cmlocks;
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
    execution_status::status_enum last_exec_status() const { return execution_status::UNSET ; }
   
    /**
     * \brief Get the number of updates executed by the engine.
     *
     * This function returns the numbe of updates executed by the last
     * run of this engine.
     * 
     * \return the total number of updates
     */
    size_t last_update_count() const { return 0; }
           
    
    /**
     * \brief Adds an update task with a particular priority.
     * This function is forwarded to the scheduler.
     */
    void schedule(vertex_id_type vid,
                  const update_functor_type& update_functor) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule " << vid << std::endl;
      procid_t owner = graph.get_vertex_record(vid).owner;
      if (owner == rmi.procid()) {
        const size_t cpuid = random::fast_uniform<size_t>(0, ncpus - 1);
        scheduler_ptr->schedule(cpuid, graph.local_vid(vid), update_functor);
        if (started && threads_alive.value < ncpus) {
          consensus->cancel_one();
        }
      }
      else {
        logstream(LOG_DEBUG) << rmi.procid() << ": Forwarding Schedule " << vid << " to " << owner << std::endl;
        rmi.remote_call(owner, &engine_type::schedule, vid, update_functor);
      }
    }


    /**
     * \brief Creates a collection of tasks on all the vertices in the
     * graph, with the same update function and priority This function
     * is forwarded to the scheduler. Must be called by all machines 
     * simultaneously
     */
    void schedule_all(const update_functor_type& update_functor) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule All" << std::endl;
      for (vertex_id_type lvid = 0; 
           lvid < graph.get_local_graph().num_vertices(); 
           ++lvid) {
        if (graph.l_get_vertex_record(lvid).owner == rmi.procid()) {
          scheduler_ptr->schedule(lvid % ncpus, lvid, update_functor);
        }
      }
      
      if (started && threads_alive.value < ncpus) {
        consensus->cancel();
      }
      rmi.barrier();
    }

    /**
     * \brief associate a termination function with this engine.
     *
     * An engine can typically have many termination functions
     * associated with it. A termination function is a function which
     * takes a constant reference to the shared data and returns a
     * boolean which is true if the engine should terminate execution.
     *
     * A termination function has the following type:
     * \code
     * bool term_fun(const ishared_data_type* shared_data)
     * \endcode
     */
    void add_termination_condition(termination_function_type term) { 
      rmi.barrier();
    }

    //!  remove all associated termination functions
    void clear_termination_conditions() { 
      rmi.barrier();
    };
    
    /**
     *  \brief The timeout is the total
     *  ammount of time in seconds that the engine may run before
     *  exeuction is automatically terminated.
     */
    void set_timeout(size_t timeout_secs) { };
    
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
    const graphlab_options& get_options() {
      return opts;
    }
    
    
    void get_a_task(size_t threadid, 
                    bool& has_internal_task,
                    vertex_id_type& internal_lvid,
                    bool& has_sched_task,
                    vertex_id_type& sched_lvid,
                    update_functor_type &task) {
      has_internal_task = false;
      has_sched_task = false;
      if (thrlocal[threadid].get_task(internal_lvid)) {
        has_internal_task = true;
        return;
      }
      sched_status::status_enum stat = scheduler_ptr->get_next(threadid, sched_lvid, task);
      if (stat != sched_status::EMPTY) {
        has_sched_task = true;
      }
    }
    
    
    void eval_internal_task(vertex_id_type lvid) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Internal Task: " << graph.global_vid(lvid) << std::endl;
      vstate_locks[lvid % VERTEX_LOCK_COUNT].lock();
      
      std::vector<vertex_id_type> ready_vertices;
      switch(vstate[lvid].state) {
        case NONE:
          break;
        case GATHERING:
        case MIRROR_GATHERING:
          {
            vertex_id_type tmp = cmlocks->make_philosopher_hungry(lvid);
            if (tmp != cmlocks->invalid_vid()) {
              ready_vertices.push_back(tmp);
            }
          }
          break;
        case APPLYING:
          break;
        case SCATTERING:
        case MIRROR_SCATTERING:
          ready_vertices = cmlocks->philosopher_stops_eating(lvid);
          break;
      }
      
      for (size_t i = 0;i < ready_vertices.size(); ++i) {
        std::cout << ready_vertices[i] << " ";
      }
      std::cout << std::endl;
      vstate_locks[lvid % VERTEX_LOCK_COUNT].unlock();
    }


    void add_internal_task(vertex_id_type lvid) {
      size_t i = random::rand() % ncpus;
      size_t j = random::rand() % ncpus;
      if (thrlocal[i].npending < thrlocal[j].npending) {
        thrlocal[i].add_task(lvid);
      }
      else {
        thrlocal[i].add_task(lvid);
      }
    }
    
    // If I receive the call I am a mirror of this vid
    void rpc_begin_gathering(vertex_id_type sched_vid, const update_functor_type& task) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Mirror Begin Gathering: " << sched_vid << std::endl;
      ASSERT_NE(graph.get_vertex_record(sched_vid).owner, rmi.procid());
      // immediately begin issuing the lock requests
      vertex_id_type sched_lvid = graph.local_vid(sched_vid);
      // set the vertex state
      vstate_locks[sched_lvid % VERTEX_LOCK_COUNT].lock();
      ASSERT_TRUE(vstate[sched_lvid].state == NONE);
      vstate[sched_lvid].state = MIRROR_GATHERING;
      vstate[sched_lvid].current = task;
      vstate_locks[sched_lvid % VERTEX_LOCK_COUNT].unlock();
      // lets go
      // add_internal_task(sched_lvid);
    }

    /**
     * Task was added to the vstate. Now to begin scheduling the gathers
     */
    void master_broadcast_gathering(vertex_id_type sched_lvid, 
                                    const update_functor_type& task) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Broadcast Gathering: " << graph.global_vid(sched_lvid) << std::endl;
      ASSERT_EQ(graph.l_get_vertex_record(sched_lvid).owner, rmi.procid());
      vertex_id_type sched_vid = graph.global_vid(sched_lvid);
      const typename graph_type::vertex_record& vrec = graph.get_vertex_record(sched_vid);
      foreach(procid_t pid, vrec.get_replicas()) {
        if (pid != rmi.procid()) {
          rmi.remote_call(pid, &engine_type::rpc_begin_gathering, sched_vid, task);
        }
      }
      add_internal_task(sched_lvid);
    }

    void eval_sched_task(size_t sched_lvid, const update_functor_type& task) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule Task: " << graph.global_vid(sched_lvid) << std::endl;
      ASSERT_EQ(graph.l_get_vertex_record(sched_lvid).owner, rmi.procid());
      // this is in local VIDs
      bool begin_gathering_vertex = false;
      vstate_locks[sched_lvid % VERTEX_LOCK_COUNT].lock();
      if (vstate[sched_lvid].state == NONE) {
        // we start gather right here.
        // set up the state
        vstate[sched_lvid].state = GATHERING;
        vstate[sched_lvid].current = task;
        vstate[sched_lvid].apply_count_down =   
                    graph.l_get_vertex_record(sched_lvid).get_replicas().size();
        // we are going to broadcast after unlock
        begin_gathering_vertex = true;

      }
      else {
        if (vstate[sched_lvid].hasnext) {
          vstate[sched_lvid].next += task;
        }
        else {
          vstate[sched_lvid].hasnext = true;
          vstate[sched_lvid].next = task;
        }
      }
      vstate_locks[sched_lvid % VERTEX_LOCK_COUNT].unlock();
      if (begin_gathering_vertex) master_broadcast_gathering(sched_lvid, task);
    }
    
    void thread_start(size_t threadid) {
      bool has_internal_task = false;
      bool has_sched_task = false;
      vertex_id_type internal_lvid;
      vertex_id_type sched_lvid;
      update_functor_type task;
      
      while(1) {
        get_a_task(threadid, 
                   has_internal_task, internal_lvid,
                   has_sched_task, sched_lvid, task);
        // if we managed to get a task..
        if (has_internal_task) {
          eval_internal_task(internal_lvid);
        }
        else if (has_sched_task) {
          eval_sched_task(sched_lvid, task);
        }
        /*
         * We failed to obtain a task, try to quit 
         */
        else if (!try_to_quit(threadid,
                              has_internal_task, internal_lvid,
                              has_sched_task, sched_lvid, task)) {
          if (has_internal_task) {
            eval_internal_task(internal_lvid);
          }
          else if (has_sched_task) {
            eval_sched_task(sched_lvid, task);
          }
        }
        else {
          break;
        }
      }
    }
    
    /**
     * \brief Start the engine execution.
     *
     * This \b blocking function starts the engine and does not
     * return until either one of the termination conditions evaluate
     * true or the scheduler has no tasks remaining.
     */
    void start() {
      ASSERT_TRUE(scheduler_ptr != NULL);
      // start the scheduler
      scheduler_ptr->start();
      started = true;
      threads_alive.value = ncpus;
      
      rmi.barrier();
      for (size_t i = 0; i < ncpus; ++i) {
        thrgroup.launch(boost::bind(&engine_type::thread_start, this, i));
      }
      thrgroup.join();
    }
  };

}

#include <graphlab/macros_undef.hpp>

#endif // GRAPHLAB_DISTRIBUTED_ENGINE_HPP

