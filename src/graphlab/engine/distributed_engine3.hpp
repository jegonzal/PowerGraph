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

#ifndef GRAPHLAB_DISTRIBUTED_ENGINE3_HPP
#define GRAPHLAB_DISTRIBUTED_ENGINE3_HPP

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
#include <graphlab/rpc/distributed_chandy_misra.hpp>

#include <graphlab/util/tracepoint.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/rpc/async_consensus.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {



  template<typename Graph, typename UpdateFunctor>
  class distributed_engine3: public iengine<Graph, UpdateFunctor> {

  public:
    // Include parent types
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef distributed_engine3<Graph, UpdateFunctor> engine_type;

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

    typedef typename iengine_base::icontext_type  icontext_type;
    typedef context<distributed_engine3>           context_type;
    typedef distributed_chandy_misra<Graph> lock_set_type;

    //typedef context_manager<distributed_engine3> context_manager_type;


    // typedef typename iengine_base::termination_function_type
    // termination_function_type;

    consistency_model context_range;

    enum vertex_execution_state {
      NONE = 0,
      LOCKING,     // state on owner
      GATHERING,   // state on owner
      APPLYING,    // state on owner
      SCATTERING,  // state on owner
      MIRROR_GATHERING, // state on mirror
      MIRROR_SCATTERING, // state on mirror
      MIRROR_SCATTERING_AND_NEXT_LOCKING, // state on mirror
    }; // end of vertex execution state


    struct vertex_state {
      uint32_t apply_count_down; // used to count down the gathers
      bool hasnext;
      vertex_execution_state state; // current state of the vertex
      update_functor_type current; // What is currently being executed
                                   //  accumulatedn
      vertex_state(): apply_count_down(0), hasnext(false), state(NONE) { }
      std::ostream& operator<<(std::ostream& os) const {
        switch(state) {
        case NONE: { os << "NONE"; break; }
        case GATHERING: { os << "GATHERING: " << apply_count_down; break; }
        case APPLYING: { os << "APPLYING"; break; }
        case SCATTERING: { os << "SCATTERING"; break; }
        case MIRROR_GATHERING: { os << "MIRROR_GATHERING"; break; }
        case MIRROR_SCATTERING: { os << "MIRROR_SCATTERING"; break; }
        case MIRROR_SCATTERING_AND_NEXT_LOCKING: { os << "MIRROR_SCATTERING_AND_NEXT_LOCKING"; break; }
        }
        return os;
      }
    }; // end of vertex_state

    // preferred order in the scheduling buffer. order from least preferred
    // to most preferred.
    struct schedule_preference_generator{
      local_graph_type* lgraph;
      lock_set_type* lockset;
      schedule_preference_generator():lgraph(NULL), lockset(NULL) { }
      schedule_preference_generator(local_graph_type* lgraph,
                                     lock_set_type* lockset):lgraph(lgraph), lockset(lockset) { }

      float operator()(vertex_id_type a) {
        float pa = (float)lockset->philosopherset[a].forks_acquired /
                   (float)lockset->philosopherset[a].num_edges;
        return pa;
        //return lgraph->color(a.first) < lgraph->color(b.first);
      }
    };
    
    struct thread_local_data {
      mutex lock;
      mutex schedlock;
      size_t extracted_sched_size;
      size_t npending;
      
      struct taskentry {
        vertex_id_type vid;
        update_functor_type uf;
        float pr;
        bool operator<(const taskentry &te) const {
          return pr < te.pr;
        }
      };
      std::vector<taskentry> extracted_schedule;
      std::deque<vertex_id_type> pending_vertices;

      schedule_preference_generator generator;
      
      thread_local_data() : extracted_sched_size(0), npending(0) { }

    
      void fill_extracted_schedule_unsync(size_t threadid, ischeduler_type* sched) {
        // fill schedule from the scheduler
        while(extracted_sched_size < extracted_schedule.size()) {
          vertex_id_type vid;
          update_functor_type task;
          sched_status::status_enum stat =
                        sched->get_next(threadid, vid, task);
          if (stat != sched_status::EMPTY) {
            extracted_schedule[extracted_sched_size].vid = vid;
            extracted_schedule[extracted_sched_size].uf = task;
            extracted_schedule[extracted_sched_size].pr = generator(vid);
            ++extracted_sched_size;
          }
          else {
            break;
          }
        }
        std::sort(extracted_schedule.begin(),
                  extracted_schedule.begin() + extracted_sched_size);
      }
    
      void init_reorder_buffer(size_t bufsize, schedule_preference_generator comp) {
        extracted_schedule.resize(bufsize);
        extracted_sched_size = 0;
        generator = comp;
      }
      
      bool get_schedule(size_t threadid,
                        ischeduler_type* sched,
                        vertex_id_type& vid,
                        update_functor_type &task) {
        schedlock.lock();
        bool retval = false;
        if (extracted_sched_size > 0) {
          // if there is something in the extracted schedule, use it
          extracted_sched_size--;
          vid = extracted_schedule[extracted_sched_size].vid;
          task = extracted_schedule[extracted_sched_size].uf;
          retval = true;
        }
        else {
          // Otherwise, try to fill it from the scheduler,
          // and try again
          fill_extracted_schedule_unsync(threadid, sched);
          if (extracted_sched_size > 0) {
            extracted_sched_size--;
            vid = extracted_schedule[extracted_sched_size].vid;
            task = extracted_schedule[extracted_sched_size].uf;
            retval = true;
          }
        }
        schedlock.unlock();
        return retval;
      }
      
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
    dc_dist_object<distributed_engine3<Graph, UpdateFunctor> > rmi;

    //! The local engine options
    graphlab_options opts;

    graph_type& graph;

    lock_set_type* cmlocks;

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

    atomic<uint64_t> joined_tasks;
    atomic<uint64_t> blocked_issues; // issued tasks which either
                                     // 1: cannot start and have to be
                                     //    reinjected into the
                                     //    scheduler.
                                     // 2: issued but is combined into a
                                     //    task which is currently locking
    atomic<uint64_t> issued_tasks;
    atomic<uint64_t> completed_tasks;

    size_t max_pending_edges;
    atomic<size_t> pending_edges;
    size_t max_clean_forks;
    size_t max_thread_reorder_buffer;

    PERMANENT_DECLARE_DIST_EVENT_LOG(eventlog);
    DECLARE_TRACER(disteng_eval_sched_task);
    DECLARE_TRACER(disteng_chandy_misra);
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
      if (issued_tasks.value != completed_tasks.value + blocked_issues.value) {
        ++ctr;
        if (ctr % 10 == 0) usleep(1);
        return false;
      }
      logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid << ": " << "Termination Attempt "
                           << completed_tasks.value << "/" << issued_tasks.value << std::endl;
      has_internal_task = false;
      has_sched_task = false;
      threads_alive.dec();
      consensus->begin_done_critical_section(threadid);

      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      if (thrlocal[threadid].get_task(internal_lvid)) {
        logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tCancelled by Internal Task"  << std::endl;
        has_internal_task = true;
        consensus->cancel_critical_section(threadid);
        threads_alive.inc();
        END_TRACEPOINT(disteng_internal_task_queue);
        return false;
      }
      END_TRACEPOINT(disteng_internal_task_queue);

      has_sched_task = thrlocal[threadid].get_schedule(threadid, scheduler_ptr,
                                                       sched_lvid, task);
      if (!has_sched_task) {
        logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tTermination Double Checked" << std::endl;
        bool ret = consensus->end_done_critical_section(threadid);
        threads_alive.inc();
        if (ret == false) {
          logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tCancelled" << std::endl;
        }
        return ret;
      } else {
        logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tCancelled by Scheduler Task" << std::endl;
        consensus->cancel_critical_section(threadid);
        threads_alive.inc();
        return false;
      }
    } // end of try to quit

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
      USER_OP_EVENT = 6
    };



    void schedule_local_next(vertex_id_type local_vid) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule " << local_vid << std::endl;
      if (started) {
        BEGIN_TRACEPOINT(disteng_scheduler_task_queue);
        scheduler_ptr->schedule_from_execution_thread(thread::thread_id(),
                                                      local_vid);
        END_TRACEPOINT(disteng_scheduler_task_queue);
      }
      else {
        scheduler_ptr->schedule(local_vid);
      }
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_EVENT, 1);
      consensus->cancel();
    }


  public:
    distributed_engine3(distributed_control &dc, graph_type& graph,
                       size_t ncpus) :
      rmi(dc, this), graph(graph), scheduler_ptr(NULL), ncpus(ncpus),
      max_pending_edges((size_t)(-1)),max_clean_forks((size_t)(-1)), max_thread_reorder_buffer(100) {
      rmi.barrier();
      // TODO: Remove context creation.
      // Added context to force compilation.
      context_type context;

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
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, SCHEDULE_FROM_REMOTE_EVENT, "Remote Schedule");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, USER_OP_EVENT, "User Ops");


      INITIALIZE_TRACER(disteng_eval_sched_task,
                        "distributed_engine3: Evaluate Scheduled Task");
      INITIALIZE_TRACER(disteng_init_gathering,
                        "distributed_engine3: Initialize Gather");
      INITIALIZE_TRACER(disteng_init_scattering,
                        "distributed_engine3: Initialize Scattering");
      INITIALIZE_TRACER(disteng_waiting_for_vstate_locks,
                      "distributed_engine3: vstate Lock Contention");
      INITIALIZE_TRACER(disteng_evalfac,
                      "distributed_engine3: Time in Factorized Update user code");
      INITIALIZE_TRACER(disteng_internal_task_queue,
                      "distributed_engine3: Time in Internal Task Queue");
      INITIALIZE_TRACER(disteng_scheduler_task_queue,
                      "distributed_engine3: Time in Scheduler Task Queue");
      INITIALIZE_TRACER(disteng_chandy_misra,
                      "distributed_engine3: Time in Chandy Misra");
    }

    /**
     * Unique to the distributed engine. This must be called prior
     * to any schedule() call.
     */
    void initialize() {
      logstream(LOG_INFO) << rmi.procid() << ": Initializing..." << std::endl;
      scheduler_ptr = scheduler_factory_type::
        new_scheduler(opts.scheduler_type,
                      opts.scheduler_args,
                      graph.get_local_graph(),
                      ncpus);
      // create initial fork arrangement based on the alternate vid mapping
      cmlocks = new lock_set_type(rmi.dc(), graph,
                                  boost::bind(&engine_type::lock_ready, this, _1),
                                  boost::bind(&engine_type::forward_cached_schedule, this, _1));
      cmlocks->compute_initial_fork_arrangement();

      vstate.resize(graph.num_local_vertices());
      vstate_locks.resize(graph.num_local_vertices());
      consensus = new async_consensus(rmi.dc(), ncpus);

      thrlocal.resize(ncpus);
      rmi.barrier();
    }

    ~distributed_engine3() {
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
    size_t last_update_count() const {
      return completed_tasks.value;
    }

    void schedule_from_remote(vertex_id_type vid,
                              const update_functor_type& update_functor) {
      const lvid_type local_vid = graph.local_vid(vid);
      BEGIN_TRACEPOINT(disteng_scheduler_task_queue);
      bool direct_injection = false;
      if (vstate[local_vid].state == LOCKING) {
        vstate_locks[local_vid].lock();
        if (vstate[local_vid].state == LOCKING) {
          vstate[local_vid].current += update_functor;
          direct_injection = true;
          joined_tasks.inc();
        }
        vstate_locks[local_vid].unlock();
      }
      if (direct_injection == false) {
        scheduler_ptr->schedule(local_vid, update_functor);
      }
      END_TRACEPOINT(disteng_scheduler_task_queue);
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_FROM_REMOTE_EVENT, 1);
      consensus->cancel();
    }

    void schedule_local(vertex_id_type local_vid ,
                        const update_functor_type& update_functor) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule " << local_vid << std::endl;
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
      consensus->cancel();
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
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule All" << std::endl;
      std::vector<vertex_id_type> vtxs;
      vtxs.reserve(graph.get_local_graph().num_vertices());
      for(lvid_type lvid = 0; lvid < graph.get_local_graph().num_vertices(); ++lvid) {
        if (graph.l_get_vertex_record(lvid).owner == rmi.procid())
          vtxs.push_back(lvid);
      }
      if(order == "shuffle")
        graphlab::random::shuffle(vtxs.begin(), vtxs.end());
      foreach(lvid_type lvid, vtxs)
        scheduler_ptr->schedule(lvid, update_functor);
      if (started) consensus->cancel();
      rmi.barrier();
    } // end of schedule all


    /**state
     * Schedule an update on all the neighbors of a particular vertex
     */
    void schedule_in_neighbors(const vertex_id_type& vertex,
                               const update_functor_type& update_fun) {
      assert(false); //TODO: IMPLEMENTh
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
    void set_timeout(size_t timeout_secs) { };

    /**
     * elapsed time in milliseconds since start was called
     */
    size_t elapsed_time() const { return 0; }

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
      if(opts.engine_args.get_option("max_pending_edges", max_pending_edges)) {
        std::cout << "Max Pending Edges: " << max_pending_edges << std::endl;
      }
      if(opts.engine_args.get_option("max_clean_forks", max_clean_forks)) {
        std::cout << "Max Clean Forks: " << max_clean_forks << std::endl;
      }

      float fraction = 0;
      if(opts.engine_args.get_option("max_pending_fraction", fraction)) {
        max_pending_edges = fraction * graph.num_local_edges();
        std::cout << "Max Pending Edges: " << max_pending_edges << std::endl;
      }
      if(opts.engine_args.get_option("max_clean_fraction", fraction)) {
        max_clean_forks = fraction * graph.num_local_edges();
        std::cout << "Max Clean Forks: " << max_clean_forks << std::endl;
      }
    
      opts.engine_args.get_option("max_thread_reorder_buffer", max_thread_reorder_buffer);
      std::cout << "Max Reorder Buffer: " << max_thread_reorder_buffer << std::endl;
    }

    /** \brief get the current engine options. */
    const graphlab_options& get_options() {
      return opts;
    }

    void lock_ready(vertex_id_type lvid) {
      logstream(LOG_DEBUG) << "Lock ready on " << "L" << lvid << std::endl;
      vstate_locks[lvid].lock();
//      ASSERT_EQ(vstate[lvid].state, (int)LOCKING);
      vstate[lvid].state = GATHERING;
      update_functor_type uf = vstate[lvid].current;
      do_init_gather(lvid);
      vstate_locks[lvid].unlock();
      master_broadcast_gathering(lvid, uf);
    }

    void get_a_task(size_t threadid,
                    bool& has_internal_task,
                    std::deque<lvid_type>& internal_lvid,
                    bool& has_sched_task,
                    lvid_type& sched_lvid,
                    update_functor_type &task) {
      has_internal_task = false;
      has_sched_task = false;

      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      if (thrlocal[threadid].get_task(internal_lvid)) {
        has_internal_task = true;
        END_TRACEPOINT(disteng_internal_task_queue);
        return;
      }
      END_TRACEPOINT(disteng_internal_task_queue);

      if( pending_edges.value > max_pending_edges ) { return; }
      if (cmlocks->num_clean_forks() >= max_clean_forks) {
          return;
        }

      has_sched_task = thrlocal[threadid].get_schedule(threadid, scheduler_ptr,
                                                      sched_lvid, task);
    } // end of get a task

    void locked_gather_complete(const lvid_type lvid) {
      // make sure that I am the owner
//      ASSERT_I_AM_OWNER(lvid);
//      ASSERT_TRUE(vstate[lvid].state == GATHERING);
//      ASSERT_GT(vstate[lvid].apply_count_down, 0);
      vstate[lvid].apply_count_down--;
      logstream(LOG_DEBUG) << rmi.procid() << ": Partial Gather Complete: "
                    << graph.global_vid(lvid) << "(" << vstate[lvid].apply_count_down << ")" << std::endl;
      if (vstate[lvid].apply_count_down == 0) {
        logstream(LOG_DEBUG) << rmi.procid() << ": Gather Complete "
                             << graph.global_vid(lvid) << std::endl;
        vstate[lvid].state = APPLYING;
        add_internal_task(lvid);
      }
    }

    void rpc_gather_complete(vertex_id_type vid, const update_functor_type& uf) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Receiving Gather Complete of "
                           << vid << std::endl;

      vertex_id_type lvid = graph.local_vid(vid);
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate_locks[lvid].lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate[lvid].current.merge(uf);
      locked_gather_complete(lvid);
      vstate_locks[lvid].unlock();
    }


    void do_apply(lvid_type lvid) {
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, 1);
      BEGIN_TRACEPOINT(disteng_evalfac);
      const vertex_id_type vid = graph.global_vid(lvid);
      logstream(LOG_DEBUG) << rmi.procid() << ": Apply On " << vid << std::endl;
      update_functor_type& ufun = vstate[lvid].current;
      context_type context(this, &graph, vid, VERTEX_CONSISTENCY);
      ufun.apply(context);
      END_TRACEPOINT(disteng_evalfac);
    }

    void do_init_gather(lvid_type lvid) {
      update_functor_type& ufun = vstate[lvid].current;
      vertex_id_type vid = graph.global_vid(lvid);
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
        foreach(const edge_type& edge, edges) ufun.gather(context, edge);
      }
      if(ufun.gather_edges() == graphlab::OUT_EDGES ||
         ufun.gather_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_out_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        foreach(const edge_type& edge, edges) ufun.gather(context, edge);
      }
      END_TRACEPOINT(disteng_evalfac);
    }

    void do_scatter(lvid_type lvid) {
      BEGIN_TRACEPOINT(disteng_evalfac);
      const vertex_id_type vid = graph.global_vid(lvid);
      update_functor_type& ufun = vstate[lvid].current;
      context_type context(this, &graph, vid, ufun.scatter_consistency());
      if(ufun.scatter_edges() == graphlab::IN_EDGES ||
         ufun.scatter_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_in_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        foreach(const edge_type& edge, edges) ufun.scatter(context, edge);
      }
      if(ufun.scatter_edges() == graphlab::OUT_EDGES ||
         ufun.scatter_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_out_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        foreach(const edge_type& edge, edges) ufun.scatter(context, edge);
      }
      END_TRACEPOINT(disteng_evalfac);
    } // end of do scatter

    void process_gather(lvid_type lvid) {
      // in theory I do not need a lock here.
      // but what the hell
//      ASSERT_TRUE(vstate[lvid].state == GATHERING ||
//                  vstate[lvid].state == MIRROR_GATHERING);
      const vertex_id_type vid = graph.global_vid(lvid);
      logstream(LOG_DEBUG) << rmi.procid() << ": Gathering on " << vid
                           << std::endl;
      do_gather(lvid);

      const procid_t vowner = graph.l_get_vertex_record(lvid).owner;
      if (vowner == rmi.procid()) {
        locked_gather_complete(lvid);
      } else {
        vstate[lvid].state = MIRROR_SCATTERING;
        logstream(LOG_DEBUG) << rmi.procid() << ": Send Gather Complete of " << vid
                             << " to " << vowner << std::endl;

        rmi.remote_call(vowner,
                        &engine_type::rpc_gather_complete,
                        graph.global_vid(lvid),
                        vstate[lvid].current);

        vstate[lvid].current = update_functor_type();
      }
    }


    void eval_internal_task(lvid_type lvid) {
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, INTERNAL_TASK_EVENT, 1);
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate_locks[lvid].lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      const typename graph_type::vertex_record& rec = graph.l_get_vertex_record(lvid);
      switch(vstate[lvid].state) {
      case NONE:
        ASSERT_MSG(false, "Empty Internal Task");
      case LOCKING: {
          BEGIN_TRACEPOINT(disteng_chandy_misra);
          cmlocks->make_philosopher_hungry_per_replica(lvid);
          END_TRACEPOINT(disteng_chandy_misra);
          break;
      }
      case GATHERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Internal Task: "
                              << graph.global_vid(lvid) << ": GATHERING(" << vstate[lvid].apply_count_down << ")" << std::endl;

          process_gather(lvid);
          break;
      }
      case MIRROR_GATHERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Internal Task: "
                              << graph.global_vid(lvid) << ": MIRROR_GATHERING" << std::endl;
          do_init_gather(lvid);
          process_gather(lvid);
          break;
        }
      case APPLYING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Internal Task: "
                              << graph.global_vid(lvid) << ": APPLYING" << std::endl;

          do_apply(lvid);
          vstate[lvid].state = SCATTERING;
          master_broadcast_scattering(lvid,
                                      vstate[lvid].current,
                                      graph.get_local_graph().vertex_data(lvid));
          // fall through to scattering
        }
      case SCATTERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Scattering: "
                              << graph.global_vid(lvid) << ": SCATTERING" << std::endl;

          do_scatter(lvid);
          completed_tasks.inc();
          pending_edges.dec(rec.num_in_edges + rec.num_out_edges);
          BEGIN_TRACEPOINT(disteng_chandy_misra);
          cmlocks->philosopher_stops_eating_per_replica(lvid);
          END_TRACEPOINT(disteng_chandy_misra);

          if (vstate[lvid].hasnext) {
            // stick next back into the scheduler
            schedule_local_next(lvid);
            vstate[lvid].hasnext = false;
          }
          vstate[lvid].state = NONE;
          vstate[lvid].current = update_functor_type();
          break;
        }
      case MIRROR_SCATTERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Scattering: "
                              << graph.global_vid(lvid) << ": MIRROR_SCATTERING" << std::endl;
          do_scatter(lvid);
          vstate[lvid].state = NONE;
          vstate[lvid].current = update_functor_type();
          cmlocks->philosopher_stops_eating_per_replica(lvid);
//          ASSERT_FALSE(vstate[lvid].hasnext);
          break;
        }
      case MIRROR_SCATTERING_AND_NEXT_LOCKING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Scattering: "
                              << graph.global_vid(lvid) << ": MIRROR_SCATTERING_AND_LOCKING" << std::endl;
          do_scatter(lvid);
          vstate[lvid].current = update_functor_type();
          vstate[lvid].state = LOCKING;
//          ASSERT_FALSE(vstate[lvid].hasnext);
          cmlocks->philosopher_stops_eating_per_replica(lvid);
          cmlocks->make_philosopher_hungry_per_replica(lvid);
          break;
        }
      }
      vstate_locks[lvid].unlock();
    } // end of eval internal task


    void add_internal_task(lvid_type lvid) {
      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      size_t i = lvid % ncpus;
      if (vstate[lvid].state == APPLYING || vstate[lvid].state == SCATTERING ||
          vstate[lvid].state == MIRROR_SCATTERING ||
          vstate[lvid].state == MIRROR_SCATTERING_AND_NEXT_LOCKING) {
        thrlocal[i].add_task_priority(lvid);
      }
      else {
        thrlocal[i].add_task(lvid);
      }
      consensus->cancel_one(i);
      END_TRACEPOINT(disteng_internal_task_queue);
    }

    void forward_cached_schedule(vertex_id_type lvid) {
      update_functor_type uf;
      const typename graph_type::vertex_record& rec = graph.l_get_vertex_record(lvid);
      if (rec.owner != rmi.procid()) {
        if (scheduler_ptr->get_specific(lvid, uf) == sched_status::NEW_TASK) {
          rmi.remote_call(rec.owner, &engine_type::schedule_from_remote, rec.gvid, uf);
        }
      }
    }

    void rpc_begin_locking(vertex_id_type sched_vid) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Mirror Begin Locking: "
                           << sched_vid << std::endl;
      // immediately begin issuing the lock requests
      vertex_id_type sched_lvid = graph.local_vid(sched_vid);
      // set the vertex state
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate_locks[sched_lvid].lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      if (vstate[sched_lvid].state == NONE) {
        vstate[sched_lvid].state = LOCKING;
        add_internal_task(sched_lvid);
      }
      else if (vstate[sched_lvid].state == MIRROR_SCATTERING) {
        vstate[sched_lvid].state = MIRROR_SCATTERING_AND_NEXT_LOCKING;
      }
/*      else {
        ASSERT_TRUE(vstate[sched_lvid].state == NONE ||
                    vstate[sched_lvid].state == MIRROR_SCATTERING);
      }*/
      vstate_locks[sched_lvid].unlock();
    }

    // If I receive the call I am a mirror of this vid
    void rpc_begin_gathering(vertex_id_type sched_vid,
                             const update_functor_type& task) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Mirror Begin Gathering: "
                           << sched_vid << std::endl;
//      ASSERT_NE(graph.get_vertex_record(sched_vid).owner, rmi.procid());
      // immediately begin issuing the lock requests
      vertex_id_type sched_lvid = graph.local_vid(sched_vid);
      // set the vertex state
      BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
      vstate_locks[sched_lvid].lock();
      END_TRACEPOINT(disteng_waiting_for_vstate_locks);
//      ASSERT_EQ(vstate[sched_lvid].state, LOCKING);

      vstate[sched_lvid].state = MIRROR_GATHERING;
      vstate[sched_lvid].current = task;
      vstate_locks[sched_lvid].unlock();
      // lets go
      add_internal_task(sched_lvid);
    }

    /**
     * Task was added to the vstate. Now to begin scheduling the gathers
     */
    void master_broadcast_locking(lvid_type sched_lvid) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Broadcast Gathering: "
                           << graph.global_vid(sched_lvid) << std::endl;
//      ASSERT_I_AM_OWNER(sched_lvid);
      vertex_id_type sched_vid = graph.global_vid(sched_lvid);
      const typename graph_type::vertex_record& vrec =
        graph.l_get_vertex_record(sched_lvid);
      const unsigned char prevkey =
        rmi.dc().set_sequentialization_key(sched_vid % 254 + 1);
      rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                      &engine_type::rpc_begin_locking, sched_vid);
      rmi.dc().set_sequentialization_key(prevkey);
    }

    /**
     * Task was added to the vstate. Now to begin scheduling the gathers
     */
    void master_broadcast_gathering(lvid_type sched_lvid,
                                    const update_functor_type& task) {
      BEGIN_TRACEPOINT(disteng_init_gathering);
      logstream(LOG_DEBUG) << rmi.procid() << ": Broadcast Gathering: "
                           << graph.global_vid(sched_lvid) << std::endl;
//      ASSERT_I_AM_OWNER(sched_lvid);
      vertex_id_type sched_vid = graph.global_vid(sched_lvid);
      const typename graph_type::vertex_record& vrec =
        graph.l_get_vertex_record(sched_lvid);
      const unsigned char prevkey =
        rmi.dc().set_sequentialization_key(sched_vid % 254 + 1);
      rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                      &engine_type::rpc_begin_gathering, sched_vid, task);
      rmi.dc().set_sequentialization_key(prevkey);
      END_TRACEPOINT(disteng_init_gathering);
      add_internal_task(sched_lvid);
    }

    void rpc_begin_scattering(vertex_id_type vid, update_functor_type task,
                              const vertex_data_type &central_vdata) {
      vertex_id_type lvid = graph.local_vid(vid);
      vstate_locks[lvid].lock();
//      ASSERT_I_AM_NOT_OWNER(lvid);
//      ASSERT_EQ(vstate[lvid].state, MIRROR_SCATTERING);
      graph.get_local_graph().vertex_data(lvid) = central_vdata;
      vstate[lvid].current = task;
      vstate_locks[lvid].unlock();
      add_internal_task(lvid);
    }

    /**
     * Task was added to the vstate. Now to begin scheduling the gathers
     */
    void master_broadcast_scattering(lvid_type sched_lvid,
                                     const update_functor_type& task,
                                     const vertex_data_type &central_vdata) {
      BEGIN_TRACEPOINT(disteng_init_scattering);
      logstream(LOG_DEBUG) << rmi.procid() << ": Broadcast Scattering: "
                           << graph.global_vid(sched_lvid) << std::endl;
//      ASSERT_I_AM_OWNER(sched_lvid);
      vertex_id_type sched_vid = graph.global_vid(sched_lvid);
      const typename graph_type::vertex_record& vrec =
        graph.l_get_vertex_record(sched_lvid);
      const unsigned char prevkey =
        rmi.dc().set_sequentialization_key(sched_vid % 254 + 1);
      rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                      &engine_type::rpc_begin_scattering,
                      sched_vid, task, central_vdata);
      rmi.dc().set_sequentialization_key(prevkey);
      END_TRACEPOINT(disteng_init_scattering);
    }

    template <bool prelocked>
    void eval_sched_task(const lvid_type sched_lvid,
                         const update_functor_type& task) {
      BEGIN_TRACEPOINT(disteng_eval_sched_task);
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule Task: "
                           << graph.global_vid(sched_lvid) << std::endl;
      // If I am not the owner just forward the task to the other
      // scheduler and return
      const typename graph_type::vertex_record& rec = graph.l_get_vertex_record(sched_lvid);
      const procid_t owner = rec.owner;
      bool acquirelock = false;
      if (owner != rmi.procid()) {
        const vertex_id_type vid = rec.gvid;
        rmi.remote_call(owner, &engine_type::schedule_from_remote, vid, task);
        return;
      }
//      ASSERT_I_AM_OWNER(sched_lvid);
      // this is in local VIDs
      issued_tasks.inc();
      if (prelocked == false) {
        BEGIN_TRACEPOINT(disteng_waiting_for_vstate_locks);
        vstate_locks[sched_lvid].lock();
        END_TRACEPOINT(disteng_waiting_for_vstate_locks);
      }
      if (vstate[sched_lvid].state == NONE) {
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, UPDATE_EVENT, 1);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, WORK_ISSUED_EVENT, rec.num_in_edges + rec.num_out_edges);

        // we start gather right here.
        // set up the state
        vstate[sched_lvid].state = LOCKING;
        vstate[sched_lvid].hasnext = false;
        vstate[sched_lvid].current = task;
        vstate[sched_lvid].apply_count_down =
          graph.l_get_vertex_record(sched_lvid).num_mirrors() + 1;
        acquirelock = true;
        pending_edges.inc(rec.num_in_edges + rec.num_out_edges);
        // we are going to broadcast after unlock
      } else if (vstate[sched_lvid].state == LOCKING) {
         blocked_issues.inc();
         vstate[sched_lvid].current += task;
         joined_tasks.inc();
      } else {
        blocked_issues.inc();
        if (vstate[sched_lvid].hasnext) {
          scheduler_ptr->place(sched_lvid, task);
          joined_tasks.inc();
        } else {
          vstate[sched_lvid].hasnext = true;
          scheduler_ptr->place(sched_lvid, task);
        }
      }
      if (prelocked == false) vstate_locks[sched_lvid].unlock();
      END_TRACEPOINT(disteng_eval_sched_task);
      if (acquirelock) {
        BEGIN_TRACEPOINT(disteng_chandy_misra);
        cmlocks->make_philosopher_hungry_per_replica(sched_lvid);
        END_TRACEPOINT(disteng_chandy_misra);
        master_broadcast_locking(sched_lvid);
      }
    }

    void thread_start(size_t threadid) {
      bool has_internal_task = false;
      bool has_sched_task = false;
      std::deque<vertex_id_type> internal_lvid;
      vertex_id_type sched_lvid;
      update_functor_type task;
//      size_t ctr = 0;
      while(1) {
/*        ++ctr;
        if (max_clean_forks != (size_t)(-1) && ctr % 10000 == 0) {
          std::cout << cmlocks->num_clean_forks() << "/" << max_clean_forks << "\n";
        }*/
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
          eval_sched_task<false>(sched_lvid, task);
        }
        /*
         * We failed to obtain a task, try to quit
         */
        else if (!try_to_quit(threadid,
                              has_internal_task, internal_lvid,
                              has_sched_task, sched_lvid, task)) {
          if (has_internal_task) {
            while(!internal_lvid.empty()) {
              eval_internal_task(internal_lvid.front());
              internal_lvid.pop_front();
            }
          } else if (has_sched_task) {
            eval_sched_task<false>(sched_lvid, task);
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
      logstream(LOG_INFO) << "Spawning " << ncpus << " threads" << std::endl;
      ASSERT_TRUE(scheduler_ptr != NULL);
      // start the scheduler
      scheduler_ptr->start();
      started = true;
      threads_alive.value = ncpus;
      graph.get_local_graph().compute_coloring();
      
      context_range =  string_to_consistency_model(opts.get_scope_type());

      ASSERT_TRUE(context_range == EDGE_CONSISTENCY ||
                  context_range == VERTEX_CONSISTENCY);
      if (context_range == EDGE_CONSISTENCY) {
        std::cout << "Consistency Model is: Edge" << std::endl;
      }
      else {
        std::cout << "Consistency Model is: Vertex" << std::endl;
      }
      for (size_t i = 0;i < ncpus; ++i) {
        thrlocal[i].init_reorder_buffer(max_thread_reorder_buffer,
                                        schedule_preference_generator(&graph.get_local_graph(),
                                                                       cmlocks));
      }
      rmi.barrier();
      for (size_t i = 0; i < ncpus; ++i) {
        thrgroup.launch(boost::bind(&engine_type::thread_start, this, i));
      }
      thrgroup.join();
      size_t ctasks = completed_tasks.value;
      rmi.all_reduce(ctasks);
      completed_tasks.value = ctasks;

      ctasks = issued_tasks.value;
      rmi.all_reduce(ctasks);
      issued_tasks.value = ctasks;

      ctasks = blocked_issues.value;
      rmi.all_reduce(ctasks);
      blocked_issues.value = ctasks;

      ctasks = joined_tasks.value;
      ctasks += scheduler_ptr->num_joins();
      rmi.all_reduce(ctasks);
      joined_tasks.value = ctasks;

      if (rmi.procid() == 0) {
        std::cout << "Completed Tasks: " << completed_tasks.value << std::endl;
        std::cout << "Issued Tasks: " << issued_tasks.value << std::endl;
        std::cout << "Blocked Issues: " << blocked_issues.value << std::endl;
        std::cout << "------------------" << std::endl;
        std::cout << "Joined Tasks: " << joined_tasks.value << std::endl;
      }
      // test if all schedulers are empty
      for (size_t i = 0;i < ncpus; ++i) {
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
} // namespace

#include <graphlab/macros_undef.hpp>

#endif // GRAPHLAB_DISTRIBUTED_ENGINE_HPP

