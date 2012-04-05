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

#include <graphlab/scheduler/vertex_functor_set.hpp>
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
    typedef context<distributed_synchronous_engine>  context_type;

    
    
    consistency_model default_consistency;
    
    
  private:
    dc_dist_object<distributed_synchronous_engine<Graph, UpdateFunctor> > rmi;

    //! The local engine options
    graphlab_options opts; 
    //! the distributed graph
    graph_type& graph;
    //! local threads object
    thread_pool threads;
    size_t ncpus;
    
    size_t max_iterations;
    
    //! Engine state
    bool started;
    bool has_schedule_entries;
    
    std::vector<spinlock> vlocks;
    vertex_functor_set<update_functor_type>* scheduleset;
    vertex_functor_set<update_functor_type>* workingset;
    
    // if scheduleset[i] has a task 
    // perthread_scheduleset_bits[i % ncpus] will have bit i / ncpus set
    std::vector<dense_bitset> perthread_scheduleset_bits;
    std::vector<dense_bitset> perthread_workingset_bits; 
    
    atomic<uint64_t> completed_tasks;
    
    PERMANENT_DECLARE_DIST_EVENT_LOG(eventlog);

    enum {
      SCHEDULE_EVENT = 0,
      UPDATE_EVENT = 1,
      USER_OP_EVENT = 2,
      ENGINE_START_EVENT = 3,
      ENGINE_STOP_EVENT = 4,     
    };
    
    enum {
      BUFFER_LIMIT = 100
    };

    inline size_t lvid_to_thread(vertex_id_type lvid) const {
      return lvid % ncpus;
    }
    
    inline size_t lvid_to_bit(vertex_id_type lvid) const {
      return lvid / ncpus;
    }
    
    inline size_t bit_and_thread_to_lvid(size_t bit, size_t threadid) const {
      return bit * ncpus + threadid;
    }
    
  public:
    distributed_synchronous_engine(distributed_control &dc, graph_type& graph,
                              size_t ncpus) : 
      rmi(dc, this), graph(graph), threads(ncpus), ncpus(ncpus), max_iterations(-1),
      started(false), has_schedule_entries(false) {

#ifdef USE_EVENT_LOG
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 3000, 
                                dist_event_log::RATE_BAR);
#else
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 3000, 
                                dist_event_log::LOG_FILE);
#endif
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, SCHEDULE_EVENT, "Schedule");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, UPDATE_EVENT, "Updates");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, USER_OP_EVENT, "User Ops");
      PERMANENT_ADD_IMMEDIATE_DIST_EVENT_TYPE(eventlog, ENGINE_START_EVENT, "Engine Start");
      PERMANENT_ADD_IMMEDIATE_DIST_EVENT_TYPE(eventlog, ENGINE_STOP_EVENT, "Engine Stop");
      
      rmi.barrier();
    }

    /**
     * Unique to the distributed engine. This must be called prior
     * to any schedule() call. 
     */
    void initialize() {
      graph.finalize();
      vlocks.resize(graph.num_local_vertices());
      perthread_scheduleset_bits.resize(ncpus);
      perthread_workingset_bits.resize(ncpus);
      for (size_t i = 0;i < ncpus; ++i) {
        perthread_scheduleset_bits[i].resize(graph.num_local_vertices() / ncpus + 1);
        perthread_workingset_bits[i].resize(graph.num_local_vertices() / ncpus + 1);
        perthread_scheduleset_bits[i].clear();
        perthread_workingset_bits[i].clear();
      }
      scheduleset = new vertex_functor_set<update_functor_type>();
      workingset = new vertex_functor_set<update_functor_type>();
      scheduleset->resize(graph.num_local_vertices());
      workingset->resize(graph.num_local_vertices());
      logstream(LOG_INFO) 
        << rmi.procid() << ": Initializing..." << std::endl;
      rmi.barrier();
    }
    
    ~distributed_synchronous_engine() {
      delete workingset;
      delete scheduleset;
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

    void workingset_merge_local(vertex_id_type local_vid ,
                     const update_functor_type& update_functor) {
      workingset->merge(local_vid, update_functor);
      perthread_workingset_bits[lvid_to_thread(local_vid)].set_bit(lvid_to_bit(local_vid));
    }

    void workingset_merge(vertex_id_type vid,
                          const update_functor_type& update_functor) {
      workingset_merge_local(graph.local_vid(vid), update_functor);
    } // end of schedule

    void workingset_add_local(vertex_id_type local_vid ,
                     const update_functor_type& update_functor) {
      workingset->add(local_vid, update_functor);
      perthread_workingset_bits[lvid_to_thread(local_vid)].set_bit(lvid_to_bit(local_vid));
    }

    void workingset_add(vertex_id_type vid,
                          const update_functor_type& update_functor) {
      workingset_add_local(graph.local_vid(vid), update_functor);
    } // end of schedule


    void schedule_local(vertex_id_type local_vid ,
                        const update_functor_type& update_functor) {
      scheduleset->add(local_vid, update_functor);
      has_schedule_entries = true;
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_EVENT, 1);
      perthread_scheduleset_bits[lvid_to_thread(local_vid)].set_bit(lvid_to_bit(local_vid));
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
        if (graph.l_get_vertex_record(lvid).owner == rmi.procid()) {
          schedule_local(lvid, update_functor);
        }
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
      if(opts.engine_args.get_option("max_iterations", max_iterations)) {
        std::cout << "Max Iterations: " << max_iterations << std::endl;
      }
    } 

    /** \brief get the current engine options. */
    const graphlab_options& get_options() { return opts; }

    void do_apply(lvid_type lvid, update_functor_type& ufun) { 
      const vertex_id_type vid = graph.global_vid(lvid);
      context_type context(this, &graph, vid, VERTEX_CONSISTENCY);
      ufun.apply(context);
      completed_tasks.inc();
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, UPDATE_EVENT, 1);
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, 1);
    }
 

    void do_init_gather(lvid_type lvid, update_functor_type& ufun) {
      const vertex_id_type vid = graph.global_vid(lvid);
      context_type context(this, &graph, vid, ufun.gather_consistency());
      ufun.init_gather(context);
    }
    
    void do_gather(lvid_type lvid, update_functor_type& ufun) { // Do gather
      const vertex_id_type vid = graph.global_vid(lvid);
      context_type context(this, &graph, vid, ufun.gather_consistency());
      if(ufun.gather_edges() == graphlab::IN_EDGES || 
         ufun.gather_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_in_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          ufun.gather(context, edges[i]);
        }
      }
      if(ufun.gather_edges() == graphlab::OUT_EDGES ||
         ufun.gather_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_out_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          ufun.gather(context, edges[i]);
        }
      }
    } // end of do_gather
    
    void do_scatter(lvid_type lvid, update_functor_type& ufun) {
      const vertex_id_type vid = graph.global_vid(lvid);
      context_type context(this, &graph, vid, ufun.scatter_consistency());
      if(ufun.scatter_edges() == graphlab::IN_EDGES || 
         ufun.scatter_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_in_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          ufun.scatter(context, edges[i]);
        }
      }
      if(ufun.scatter_edges() == graphlab::OUT_EDGES ||
         ufun.scatter_edges() == graphlab::ALL_EDGES) {
        const local_edge_list_type edges = graph.l_out_edges(lvid);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, edges.size());
        for(size_t i = 0; i < edges.size(); ++i) {
          ufun.scatter(context, edges[i]);
        }
      }
    } // end of do scatter
    
    void workingset_add_batch(const std::vector<std::pair<vertex_id_type, update_functor_type> >& batch) {
      for (size_t i = 0;i < batch.size(); ++i) {
        workingset_add(batch[i].first, batch[i].second);
      }
    }


    void workingset_merge_batch(const std::vector<std::pair<vertex_id_type, update_functor_type> >& batch) {
      for (size_t i = 0;i < batch.size(); ++i) {
        workingset_merge(batch[i].first, batch[i].second);
      }
    }
    /**
     * Picks up all tasks on the working set for vertices I do 
     * not own and forward them to their owners
     */
    void transmit_schedule(size_t threadid) {
      update_functor_type uf;
      std::vector<std::vector<std::pair<vertex_id_type, update_functor_type> > > outbuffer;
      outbuffer.resize(rmi.numprocs());
      foreach(uint32_t b, perthread_workingset_bits[threadid]) {
        // convert back to local vid
        vertex_id_type local_vid = bit_and_thread_to_lvid(b, threadid);
        procid_t owner = graph.l_get_vertex_record(local_vid).owner;
        // if I am not the owner,  and there is a task, and I managed to 
        // extract it
        if (owner != rmi.procid() &&
            workingset->has_task(local_vid) && 
            workingset->test_and_get(local_vid, uf)) {
          // send it on to the right owner
          vertex_id_type global_vid = graph.global_vid(local_vid);
          outbuffer[owner].push_back(std::make_pair(global_vid, uf));
          if (outbuffer[owner].size() > BUFFER_LIMIT) {
            rmi.remote_call((procid_t)owner,
                        &distributed_synchronous_engine::workingset_add_batch, 
                        outbuffer[owner]);
            outbuffer[owner].clear();
          }
          
          // this operation may clear working set entries on vertices I do not own
          // note that this may only be called once for any given bit in this phase.
          // since I may not receive workingset_adds to vertices I do not own 
          perthread_workingset_bits[threadid].clear_bit(b);
        }
      }
      // forward to destination
      for (size_t i = 0; i < outbuffer.size(); ++i) {
        if (outbuffer[i].size() > 0) {
          rmi.remote_call((procid_t)i,
                          &distributed_synchronous_engine::workingset_add_batch, 
                          outbuffer[i]);
        }
      }
    }
    
      /**
     * Conceptual inverse of transmit_schedule.
     * Picks up all tasks on vertices I own, and send them to their 
     * mirrors.
     */
    void init_gather(size_t threadid) {
      update_functor_type* uf = NULL;
      std::vector<std::vector<std::pair<vertex_id_type, update_functor_type> > > outbuffer;
      outbuffer.resize(rmi.numprocs());
      // 
      foreach(uint32_t b, perthread_workingset_bits[threadid]) {
        // convert back to local vid
        
        vertex_id_type local_vid = bit_and_thread_to_lvid(b, threadid);
        const typename graph_type::vertex_record& vrec = graph.l_get_vertex_record(local_vid);
        // if I am the owner
        // and there are other mirrors
        // try to read a task
        if (vrec.owner == rmi.procid() && 
            workingset->get_reference_unsync(local_vid, uf)) {
          do_init_gather(local_vid, *uf);
        
          if (!vrec.mirrors().empty() && uf->gather_edges() != graphlab::NO_EDGES) {
            // loop it to the mirrors
            vertex_id_type global_vid = graph.global_vid(local_vid);

                    
/*            rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                            &distributed_synchronous_engine::workingset_add,
                            global_vid,
                            *uf);*/
            foreach(uint32_t m, vrec.mirrors()) {
              outbuffer[m].push_back(std::make_pair(global_vid, *uf));
              if (outbuffer[m].size() > BUFFER_LIMIT) {
                rmi.remote_call((procid_t)m,
                            &distributed_synchronous_engine::workingset_add_batch, 
                            outbuffer[m]);
                outbuffer[m].clear();
              }
            }
          }
        }
      }
      // forward to destination
      for (size_t i = 0; i < outbuffer.size(); ++i) {
        if (outbuffer[i].size() > 0) {
          rmi.remote_call((procid_t)i,
                          &distributed_synchronous_engine::workingset_add_batch, 
                          outbuffer[i]);
        }
      }
    }

    void perform_gather(size_t threadid) {
      update_functor_type uf;
      std::vector<std::vector<std::pair<vertex_id_type, update_functor_type> > > outbuffer;
      outbuffer.resize(rmi.numprocs());
      foreach(uint32_t b, perthread_workingset_bits[threadid]) {
        // convert back to local vid
        vertex_id_type local_vid = bit_and_thread_to_lvid(b, threadid);
        const typename graph_type::vertex_record& vrec = graph.l_get_vertex_record(local_vid);
        // If there is a task
        if (workingset->test_and_get(local_vid, uf)) {
          do_gather(local_vid, uf);
          // send it back to the owner if I am not the owner
          // otherwise just merge it back
          if (vrec.owner != rmi.procid()) {
            vertex_id_type global_vid = graph.global_vid(local_vid);
            outbuffer[vrec.owner].push_back(std::make_pair(global_vid, uf));
            if (outbuffer[vrec.owner].size() > BUFFER_LIMIT) {
              rmi.remote_call((procid_t)vrec.owner,
                          &distributed_synchronous_engine::workingset_merge_batch, 
                          outbuffer[vrec.owner]);
              outbuffer[vrec.owner].clear();
            }
            // this operation may clear working set entries on vertices I do not own
            // note that this may only be called once for any given bit in this phase.
            // since I may not receive workingset_adds to vertices I do not own 
            perthread_workingset_bits[threadid].clear_bit(b);
          }
          else {
            workingset_merge_local(local_vid, uf);
          }
        }
      }
      // forward to destination
      for (size_t i = 0; i < outbuffer.size(); ++i) {
        if (outbuffer[i].size() > 0) {
          rmi.remote_call((procid_t)i,
                          &distributed_synchronous_engine::workingset_merge_batch, 
                          outbuffer[i]);
        }
      }
    }
    
    struct apply_scatter_data {
      apply_scatter_data() { }
      apply_scatter_data(const vertex_data_type& vdata,
                         const update_functor_type& uf,
                         vertex_id_type global_vid):
                         vdata(vdata),uf(uf),global_vid(global_vid), has_uf(true) { }
      apply_scatter_data(vertex_id_type global_vid, 
                         const vertex_data_type& vdata):
                         vdata(vdata),global_vid(global_vid), has_uf(false) { }
      vertex_data_type vdata;
      update_functor_type uf;
      vertex_id_type global_vid;
      bool has_uf;
      
      void save(oarchive& oarc) const {
        if (has_uf) {
          oarc << has_uf << vdata << uf << global_vid;
        }
        else {
          oarc << has_uf << vdata << global_vid;
        }
      }
      void load(iarchive& iarc) {
        iarc >> has_uf;
        if (has_uf) {
          iarc >> vdata >> uf >> global_vid;
        }
        else {
          iarc >> vdata >> global_vid;
        }
      }
    };
    
    void update_vertex_data(vertex_id_type global_vid, 
                            vertex_data_type& vdata) {
      lvid_type lvid = graph.local_vid(global_vid);
      graph.get_local_graph().vertex_data(lvid) = vdata;
    }

    void workingset_add_with_data(vertex_id_type global_vid, 
                                  vertex_data_type& vdata,
                                  update_functor_type& uf) {
      lvid_type lvid = graph.local_vid(global_vid);
      graph.get_local_graph().vertex_data(lvid) = vdata;
      workingset_add_local(lvid, uf);
    }

    void workingset_add_batch_with_data(const std::vector<apply_scatter_data>& batch) {
      for (size_t i = 0;i < batch.size(); ++i) {
        lvid_type lvid = graph.local_vid(batch[i].global_vid);
        graph.get_local_graph().vertex_data(lvid) = batch[i].vdata;
        if (batch[i].has_uf)  workingset_add_local(lvid, batch[i].uf);
      }
    }
    
    void perform_apply_and_issue_scatter(size_t threadid) {
      update_functor_type uf;
      std::vector<std::vector<apply_scatter_data> > outbuffer;
      outbuffer.resize(rmi.numprocs());
      foreach(uint32_t b, perthread_workingset_bits[threadid]) {
        // convert back to local vid
        vertex_id_type local_vid = bit_and_thread_to_lvid(b, threadid);
        const typename graph_type::vertex_record& vrec = graph.l_get_vertex_record(local_vid);
        // If there is a task
        if (vrec.owner == rmi.procid() && 
            workingset->test_and_get(local_vid, uf)) {
          perthread_workingset_bits[threadid].clear_bit(b);
          do_apply(local_vid, uf);
          
          if (uf.scatter_edges() != graphlab::NO_EDGES) {
            // scatter!
            workingset_add_local(local_vid, uf);
/*            rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                &distributed_synchronous_engine::workingset_add_with_data,
                vrec.gvid,
                graph.get_local_graph().vertex_data(local_vid),
                uf);*/
            foreach(uint32_t m, vrec.mirrors()) {
              outbuffer[m].push_back(apply_scatter_data(graph.get_local_graph().vertex_data(local_vid),
                                                        uf, vrec.gvid));
              if (outbuffer[m].size() > BUFFER_LIMIT) {
                rmi.remote_call((procid_t)m,
                            &distributed_synchronous_engine::workingset_add_batch_with_data, 
                            outbuffer[m]);
                outbuffer[m].clear();
              }
            }
          }
          else {
            // only update data
            /*rmi.remote_call(vrec.mirrors().begin(), vrec.mirrors().end(),
                &distributed_synchronous_engine::update_vertex_data,
                vrec.gvid,
                graph.get_local_graph().vertex_data(local_vid));*/
            foreach(uint32_t m, vrec.mirrors()) {
              outbuffer[m].push_back(apply_scatter_data(vrec.gvid, 
                                                        graph.get_local_graph().vertex_data(local_vid)));
              if (outbuffer[m].size() > BUFFER_LIMIT) {
                rmi.remote_call((procid_t)m,
                            &distributed_synchronous_engine::workingset_add_batch_with_data, 
                            outbuffer[m]);
                outbuffer[m].clear();
              }
            }
          }
        }
      }
      // forward to destination
      for (size_t i = 0; i < outbuffer.size(); ++i) {
        if (outbuffer[i].size() > 0) {
          rmi.remote_call((procid_t)i,
                          &distributed_synchronous_engine::workingset_add_batch_with_data, 
                          outbuffer[i]);
        }
      }
    }
          


    void perform_scatter(size_t threadid) {
      update_functor_type uf;
      foreach(uint32_t b, perthread_workingset_bits[threadid]) {
        // convert back to local vid
        vertex_id_type local_vid = bit_and_thread_to_lvid(b, threadid);
        // If there is a task
        if (workingset->test_and_get(local_vid, uf)) {
          perthread_workingset_bits[threadid].clear_bit(b);
          do_scatter(local_vid, uf);
        }
      }
      ASSERT_TRUE(perthread_workingset_bits[threadid].empty());
    }
    
    /**
     * Transmits all schedules to remote machines.
     * Needs a full barrier for completion.
     */
    void parallel_transmit_schedule() {
      std::swap(scheduleset, workingset);
      std::swap(perthread_scheduleset_bits, perthread_workingset_bits);
      for (size_t i = 0;i < ncpus; ++i) {
        threads.launch(boost::bind(&distributed_synchronous_engine::transmit_schedule,
                                   this,
                                   i));
      }
      threads.join();
    }
    
    /**
     * Activates all scheduled tasks and initializes gathering.
     * Needs a full barrier for completion.
     */
    void parallel_init_gather() {
      for (size_t i = 0;i < ncpus; ++i) {
        threads.launch(boost::bind(&distributed_synchronous_engine::init_gather,
                                   this,
                                   i));
      }
      threads.join();
    }
    
    void parallel_gather() {
      // go through all the vertices
      for (size_t i = 0;i < ncpus; ++i) {
        threads.launch(boost::bind(&distributed_synchronous_engine::perform_gather,
                                   this,
                                   i));
      }
      threads.join();
    }


    void parallel_apply_and_issue_scatter() {
      // go through all the vertices
      for (size_t i = 0;i < ncpus; ++i) {
        threads.launch(boost::bind(&distributed_synchronous_engine::perform_apply_and_issue_scatter,
                                   this,
                                   i));
      }
      threads.join();
    }


    void parallel_scatter() {
      // go through all the vertices
      for (size_t i = 0;i < ncpus; ++i) {
        threads.launch(boost::bind(&distributed_synchronous_engine::perform_scatter,
                                   this,
                                   i));
      }
      threads.join();
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

      started = true;
      rmi.barrier();
      if (rmi.procid() == 0) {
        PERMANENT_IMMEDIATE_DIST_EVENT(eventlog, ENGINE_START_EVENT);
      }

      parallel_transmit_schedule();
      rmi.full_barrier();
      size_t iterationnumber = 1;
      for (size_t i = 0;i < max_iterations; ++i) {
        if (rmi.procid() == 0) std::cout << "Iteration " << iterationnumber << std::endl;
        ++iterationnumber;
        has_schedule_entries = false;
        parallel_init_gather();
        rmi.full_barrier();
        
        parallel_gather();
        rmi.full_barrier();
        
        parallel_apply_and_issue_scatter();
        rmi.full_barrier();
        
        parallel_scatter();
        rmi.all_reduce(has_schedule_entries);
        if (has_schedule_entries) {
          parallel_transmit_schedule();
        }
        rmi.full_barrier();
        if (has_schedule_entries == false) break;
      }
      if (rmi.procid() == 0) {
        PERMANENT_IMMEDIATE_DIST_EVENT(eventlog, ENGINE_STOP_EVENT);
      }

      size_t ctasks = completed_tasks.value;
      rmi.all_reduce(ctasks);
      completed_tasks.value = ctasks;
      if (rmi.procid() == 0) {
        std::cout << "Completed Tasks: " << completed_tasks.value << std::endl;
      }
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

