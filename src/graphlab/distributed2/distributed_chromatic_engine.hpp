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


#ifndef DISTRIBUTED_CHROMATIC_ENGINE_HPP
#define DISTRIBUTED_CHROMATIC_ENGINE_HPP

#include <functional>
#include <algorithm>
#include <ext/functional> // for select1st
#include <boost/bind.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>

#include <graphlab/engine/iengine.hpp>
#include <graphlab/logger/logger.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/distributed2/graph/distributed_graph_context.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {



/**
All processes must receive the same options at the same time.
i.e. if set_cpu_affinities is called, all processes mus call it at the same time.
This is true for all set_* functions.
*/
template<typename Graph, typename UpdateFunctor>
class distributed_chromatic_engine : public iengine<Graph> {
 public:
  typedef iengine<Graph> iengine_base;
  typedef typename iengine_base::update_task_type update_task_type;
  typedef typename iengine_base::update_function_type update_function_type;
  typedef typename iengine_base::termination_function_type termination_function_type;
  typedef typename iengine_base::iscope_type iscope_type;
  
  typedef typename iengine_base::sync_function_type sync_function_type;
  typedef typename iengine_base::merge_function_type merge_function_type;


  typedef typename distributed_chromatic_engine<Graph, UpdateFunctor> engine_type;

 private:
  // the local rmi instance
  dc_dist_object<engine_type> rmi;
  
  // the graph we are processing
  Graph &graph;
  
  /** Number of cpus to use */
  size_t ncpus; 

  /** Use processor affinities */
  bool use_cpu_affinity;

  /** Use schedule yielding when waiting on the scheduler*/
  bool use_sched_yield;
 
  
  /** Track the number of updates */
  std::vector<size_t> update_counts;

  atomic<size_t> numsyncs;

  /** terminators */ 
  std::vector<termination_function_type> term_functions;

  /** The timeout time in millis */
  size_t timeout_millis;
  timer ti;
  
  /// Used to identify when the engine is stopped through stop
  bool force_stop;
  
  /** The total number of tasks that should be executed */
  size_t task_budget;
  
  size_t randomize_schedule;
  
  /** If dynamic scheduling is used, the number of scheduled tasks */
  atomic<size_t> num_pending_tasks;
  
  
  /** The cause of the last termination condition */
  exec_status termination_reason;

  scope_range::scope_range_enum default_scope_range;
 
  std::vector<std::vector<vertex_id_t> > color_block; // set of localvids in each color
  
  vertex_functor_set<engine_type> vfun_set;
  size_t max_iterations;
  double barrier_time;
  size_t num_dist_barriers_called;
  
  // other optimizations
  bool const_nbr_vertices, const_edges;

  // Global Values ----------------------------------------------------------
  struct global_record { 
    std::vector<spinlock> locks;
    graphlab::any values;
  }; // end of global_record
  typedef std::map<std::string, global_record> global_map_type;
  global_map_type global_records;

  // Global Aggregates ------------------------------------------------------
  struct isync {
    size_t interval;
    size_t next;
    vertex_id_type begin_vid, end_vid;
    isync():interval(0),next(0) { }
    virtual ~isync() { }
    virtual void per_thread_accum(dgraph_context<engine_type> &context, size_t threadid) = 0;
    virtual void join_thread_accum() = 0;
    virtual any get_accum() = 0;
    virtual void join_other_accum(any& other) = 0;
    virtual void finalize(dgraph_context<engine_type>& globalcontext) = 0;
  }; // end of isync
  
  template<typename Accum >
  struct sync : public isync {
    typedef Accum       accumulator_type;
    using isync::begin_vid;
    using isync::end_vid;
    const accumulator_type zero;
    std::vector<acccumulator_type> intermediate_accumulator;
    accumulator_type shared_accumulator;
    sync(const accumulator_type& zero, size_t ncpus) : zero(zero), 
                                         intermediate_accumulator(ncpus, zero),
                                         shared_accumulator(zero){ }
                                         
    void per_thread_accum(dgraph_context<engine_type> &context, size_t threadid) {
      intermediate_accumulator[threadid](context);
    }
    void join_thread_accum() {
      shared_accumulator = zero;
      for (size_t i = 0;i < intermediate_accumulator.size(); ++i) {
        shared_accumulator += intermediate_accumulator[i];
        intermediate_accumulator[i] = zero;
      }
    }
    any get_accum() {
      return shared_accumulator;
    }
    virtual void join_other_accum(any& other) {
      shared_accumulator += other.as<Accum>();
    }
    
    virtual void finalize(dgraph_context<engine_type>& globalcontext) {
      shared_accumulator.finalize(globalcontext);
    }
    
  }; // end of sync

  mutex sync_master_lock;
  std::map<std::string, isync*> sync_map;
  
  metrics engine_metrics;  

 public:
  distributed_chromatic_engine(distributed_control &dc,
                                    Graph& graph,
                                    size_t ncpus = 1):
                            rmi(dc, this),
                            graph(graph),
                            callback(this),
                            ncpus( std::max(ncpus, size_t(1)) ),
                            use_cpu_affinity(false),
                            use_sched_yield(true),
                            update_counts(std::max(ncpus, size_t(1)), 0),
                            timeout_millis(0),
                            force_stop(false),
                            task_budget(0),
                            randomize_schedule(0),
                            termination_reason(EXEC_UNSET),
                            vfun_set(graph.owned_vertices().size()),
                            max_iterations(0),
                            barrier_time(0.0),
                            const_nbr_vertices(true),
                            const_edges(false),
                            engine_metrics("engine"),
                            thread_color_barrier(ncpus){ 
    rmi.barrier();
  }
  
  ~distributed_chromatic_engine() {
    rmi.barrier();
  }
  
  
  //! Get the number of cpus
  size_t get_ncpus() const { return ncpus; }

  //! set sched yield
  void set_sched_yield(bool value) {
    use_sched_yield = value;
    rmi.barrier();
  }

  void set_cpu_affinities(bool value) {
    use_cpu_affinity = value;
    rmi.barrier();
  }


  /**
   * Set the default scope range.  The scope ranges are defined in
   * iscope.hpp
   */
  void set_default_scope(scope_range::scope_range_enum default_scope_range_) {
    default_scope_range = default_scope_range_;
    rmi.barrier();
  }
  
  using iengine<Graph>::exec_status_as_string;   
  

  /**
   * Return the reason why the engine last terminated
   */
  exec_status last_exec_status() const {
    return termination_reason;
  }

  /**
   * This function computes the last update count by adding all the
   * update counts of the individual threads.  This is an underestimate
   * if the engine is currently running.
   */
  size_t thisproc_update_counts() const {
    size_t sum = 0;
    for(size_t i = 0; i < update_counts.size(); ++i)
      sum += update_counts[i];
    return sum;
  } 


  /**
   * Returns the total number of updates executed
   */
  size_t last_update_count() const {
    return total_update_count;
  } // end of last_update_count


  /**
   * Add a terminator to the engine.
   * Must be called by all machines simultaneously.
   */
  void add_terminator(termination_function_type term) {
    term_functions.push_back(term);
    rmi.barrier();
  }


  /**
   * Clear all terminators from the engine
   * Must be called by all machines simultaneously.
   */
   void clear_terminators() {
    term_functions.clear();
    rmi.barrier();
  }



  /**
   * Set a timeout. Disabled if set to 0.
   * Must be called by all machines simultaneously.
   */
  void set_timeout(size_t timeout_seconds = 0) {
    timeout_millis = timeout_seconds * 1000;
    rmi.barrier();
  }
  
  /**
   * Sets a Task budget - max number of tasks to allow.
   * Disabled if set to 0.
   * Must be called by all machines simultaneously.
   */
  void set_task_budget(size_t max_tasks) {
    task_budget = max_tasks;
    rmi.barrier();
  }
  

private:
  /** brief Adds an update task to local vertices 
   * with a particular priority
   */
  void schedule_one(vertex_id_type vid,
                const update_functor_type& update_functor) {
    vertex_id_t localvid = graph.globalvid_to_localvid(vid);
    if (vfun_set.add(localvid, fun)) num_pending_tasks.inc();
  }
public:

  /**
   * \brief Adds an update task with a particular priority.
   * add_task on any vertex can be called by any machine.
   * The call is asynchronous and may not be completed until
   * a full_barrier is issued.
   */
  void schedule(vertex_id_type vid,
                const update_functor_type& update_functor) {
    if (graph.is_owned(vid)) {
      schedule_one(vid, update_functor);
    }
    else {
      rmi.remote_call(graph.globalvid_to_owner(vid),
                      &engine_type::schedule_one,
                      vid,
                      update_functor);
    }
  }

  void schedule_in_neighbors(vertex_id_type vid,
                            const update_functor_type& update_functor) {
    if (graph.is_owned(vid)) {
      std::vector<vertex_id_type> inv = graph.in_vertices(vid);
      for (size_t i = 0;i < inv.size(); ++i) {
        schedule_one(inv[i], update_functor);
      }
    }
    else {
      rmi.remote_call(graph.globalvid_to_owner(vid),
                      &engine_type::schedule_in_neighbors,
                      vid,
                      update_functor);
    }
  }

  void schedule_out_neighbors(vertex_id_type vid,
                              const update_functor_type& update_functor) {
    if (graph.is_owned(vid)) {
      std::vector<vertex_id_type> outv = graph.out_vertices(vid);
      for (size_t i = 0;i < outv.size(); ++i) {
        schedule_one(outv[i], update_functor);
      }
    }
    else {
      rmi.remote_call(graph.globalvid_to_owner(vid),
                      &engine_type::schedule_out_neighbors,
                      vid,
                      update_functor);
    }
  }

  //! \brief Adds an update task with a particular priority.
  void schedule(const std::vector<vertex_id_type>& vid,
                const update_functor_type& update_functor) {
    // not the most effective way to do it
    for (size_t i = 0;i < vid.size(); ++i) {
      schedule(vid[i], update_functor);
    }
  }



 /**
   * \brief Creates a collection of tasks on all the vertices in the graph,
   * with the same update function and priority
   * Must be called by all machines simultaneously
   */
  void schedule_all(const update_functor_type& update_functor) {
    // local vertex IDs are all consecutive
    for (size_t i = 0;i < graph.owned_vertices().size(); ++i) {
      vertex_id_t localvid = (vertex_id_t)(i);
      if (vfun_set.add(localvid, fun)) num_pending_tasks.inc();
    }
    rmi.barrier();
  }



  template<typename T>
  void add_global(const std::string& key, const T& value, size_t size = 1) {
    global_record& record = global_records[key];
    // Set the initial value (this can change the type)
    typedef std::vector<T> vector_type;
    record.values = vector_type(size, value);
    record.locks.resize(size);
  } //end of set_global


  template<typename T>
  void set_global(const std::string& key, const T& value, size_t index = 0) {
    typename global_map_type::iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;
      return;
    }
    global_record& record = iter->second;
    typedef std::vector<T> vector_type;
    
    // graphlab::any& any_ref = record.values;
    // vector_type& values = any_ref.as<vector_type>();
    vector_type& values = record.values.template as<vector_type>();
   
    ASSERT_EQ(values.size(), record.locks.size());
    ASSERT_LT(index, values.size());
    record.locks[index].lock();
    values[index] = value; 
    // broadcast
    record.locks[index].unlock();
    for(procid p = 0;p < rmi.numprocs(); ++p) {
      if (p != rmi.procid()) {
        rmi.remote_call(&distributed_chromatic_engine<Graph, UpdateFunctor>::set_global<T>,
                        key,
                        value,
                        index);
      }
    }
  } // end of set_global


  template<typename T>
  void get_global(const std::string& key, T& ret_value, size_t index = 0) const {
    typename global_map_type::const_iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;
      return;
    }
    const global_record& record = iter->second;
    typedef std::vector<T> vector_type;
    // const graphlab::any& any_ref = record.values;
    // const vector_type& values = any_ref.as<vector_type>();
    const vector_type& values = record.values.template as<vector_type>();
    ASSERT_EQ(values.size(), record.locks.size());
    ASSERT_LT(index, values.size());
    record.locks[index].lock();
    ret_value = values[index];
    record.locks[index].unlock();
  } //end of get_global


  template<typename Accum>
  void add_sync(const std::string& key,
           const Accum& zero,                 
           size_t sync_interval,
           bool use_barrier,
           vertex_id_type begin_vid,
           vertex_id_type end_vid) {
    isync*& sync_ptr = sync_map[key];
    // Clear the old sync and remove from scheduling queue
    if(sync_ptr != NULL) { delete sync_ptr; sync_ptr = NULL; }
    sync_queue.remove(key);
    ASSERT_TRUE(sync_ptr == NULL);
    // Attach a new sync type
    typedef sync<Accum> sync_type;
    sync_ptr = new sync_type(zero, ncpus);
    sync_ptr->interval    = sync_interval;
    sync_ptr->next = 0;
    sync_ptr->begin_vid   = begin_vid;
    sync_ptr->end_vid     = end_vid;
  }// end of add_sync



  
  void generate_color_blocks() {
    // construct for each color, the set of vertices as well as the 
    // number of replicas for that vertex.
    // the number of replicas - 1 is the amount of communication
    // we have to perform to synchronize modifications to that vertex


    std::vector<std::vector<std::pair<size_t, vertex_id_t> > > color_block_and_weight;
    const size_t num_colors(graph.recompute_num_colors());
    // the list of vertices for each color
    color_block_and_weight.resize(num_colors);
    
    foreach(vertex_id_t v, graph.owned_vertices()) {
      color_block_and_weight[graph.get_color(v)].push_back(
                                    std::make_pair(graph.globalvid_to_replicas(v).size(), 
                                              graph.globalvid_to_localvid(v)));
    }
    color_block.clear();
    color_block.resize(num_colors);
    if (randomize_schedule) {
      for (size_t i = 0; i < color_block_and_weight.size(); ++i) {
        random::shuffle(color_block_and_weight[i].begin(),
                            color_block_and_weight[i].end());
      }
    }
    else {
      // optimize ordering. Sort in descending order
      // put all those which need a lot of communication in the front
      // to give communication the maximum amount if time possible.
      for (size_t i = 0; i < color_block_and_weight.size(); ++i) {
        std::sort(color_block_and_weight[i].rbegin(),
                  color_block_and_weight[i].rend());
      }
    }

    // insert the sorted vertices into the final color_block
    for (size_t i = 0;i < color_block_and_weight.size(); ++i ) {  
      std::transform(color_block_and_weight[i].begin(),
                    color_block_and_weight[i].end(), 
                    std::back_inserter(color_block[i]),
                    __gnu_cxx::select2nd<std::pair<size_t, vertex_id_t> >());
    }
    
  }
  
  /************  Actual Execution Engine ****************/
 private:

  atomic<size_t> curidx;
  barrier thread_color_barrier;
 public: 
  
  struct termination_evaluation{
    size_t pending_tasks;
    size_t executed_tasks;
    bool terminator;
    bool timeout;
    bool force_stop;
    termination_evaluation(): pending_tasks(0),
                              executed_tasks(0),
                              terminator(false),
                              timeout(false),
                              force_stop(false) { }
                              
    void save(oarchive &oarc) const {
      oarc << pending_tasks
           << executed_tasks
           << terminator
           << timeout
           << force_stop;
    }
    
    void load(iarchive &iarc) {
      iarc >> pending_tasks
           >> executed_tasks
           >> terminator
           >> timeout
           >> force_stop;
    }
  };


  /**
   * Called whenever a vertex is executed.
   * Accumulates the available syncs
   */
  void eval_syncs(vertex_id_t curvertex, dgraph_context<engine_type>& context, size_t threadid) {
    typedef std::map<std::string, isync*>::value_type pair_type;
    // go through all the active sync tasks
    foreach(const pair_type& pair, sync_map) {        
      // if in range, sync!
      if (pair.second->next <= num_executed_tasks && 
          pair.second->begin-vid <= curvertex && curvertex < pair.second->end_vid) {
        pair.second->per_thread_accum(context, threadid);
      }
    }
  }

  /** Called at the end of the iteration. Called by all threads after a barrier*/
  void sync_end_iteration(size_t threadid) {
    if (threadid == 0) {
      // one thread of each machine participates in |active_sync_tasks| gathers
      foreach(const pair_type& pair, sync_map) {        
        if (pair.second->next <= num_executed_tasks) {
          pair.second->join_thread_accum();
        }
      }
      foreach(const pair_type& pair, sync_map) {        
        isync* task = pair.second;
        if (task->next <= num_executed_tasks) {
          std::vector<any> gathervals(rmi.numprocs());
          gathervals[rmi.procid()] = task->get_accum();
          rmi.gather(gathervals, 0);
  
          // now if I am target I need to do the final merge and apply
          if (rmi.procid() == 0) {
            for (size_t i = 1; i < gathervals.size(); ++i) {
              task->join_other_accum(gathervals[i]);
            }
            dgraph_context<engine_type> globalctx(&graph, this, 0);
            task->finalize(globalctx);
            numsyncs.inc();
          }
          
          task->next = num_executed_tasks + task->interval;
        }
      }
    }
    thread_color_barrier.wait();
  }

  /** Checks all machines for termination and sets the termination reason.
      Also returns the number of update tasks completed globally */
  size_t check_global_termination(bool check_dynamic_schedule) {
    std::vector<termination_evaluation> termination_test;
    termination_test.resize(rmi.numprocs());
    
    if (check_dynamic_schedule) {
      termination_test[rmi.procid()].pending_tasks = num_pending_tasks.value;
    }

    size_t numupdates = 0;
    for (size_t i = 0; i < update_counts.size(); ++i) numupdates += update_counts[i];
    termination_test[rmi.procid()].executed_tasks = numupdates;
  
    if (timeout_millis > 0 && ti.current_time_millis() > timeout_millis) {
      termination_test[rmi.procid()].timeout = true;
    }
    
    for (size_t i = rmi.procid(); i < term_functions.size(); i += rmi.numprocs()) {
      if (term_functions[i]()) {
        termination_test[rmi.procid()].terminator = true;
        break;
      }
    }
    termination_test[rmi.procid()].force_stop = force_stop;
    // gather all to 0.
    // machine 0 evaluates termiation
    rmi.gather(termination_test, 0);
    // used to globally evaluate termination
    termination_evaluation aggregate;
    if (rmi.procid() == 0) {
      for (size_t i = 0;i < termination_test.size(); ++i) {
        aggregate.pending_tasks += termination_test[i].pending_tasks;
        aggregate.executed_tasks += termination_test[i].executed_tasks;
        aggregate.terminator |= termination_test[i].terminator;
        aggregate.timeout |= termination_test[i].timeout;
        aggregate.force_stop |= termination_test[i].force_stop;
      }
      
      if (check_dynamic_schedule && aggregate.pending_tasks == 0) {
        termination_reason = EXEC_TASK_DEPLETION;
      }
      else if (task_budget > 0 && aggregate.executed_tasks >= task_budget) {
        termination_reason = EXEC_TASK_BUDGET_EXCEEDED;
      }
      else if (timeout_millis > 0 && aggregate.timeout) {
        termination_reason = EXEC_TIMEOUT;
      }
      else if (aggregate.terminator) {
        termination_reason = EXEC_TERM_FUNCTION;
      }
      else if (aggregate.force_stop) {
        termination_reason = EXEC_FORCED_ABORT;
      }
    }
    size_t treason = termination_reason;
    // note this is OK because only machine 0 will have the right value for
    // executed_tasks. And everyone is receiving from machine 0
    std::pair<size_t, size_t> reason_and_task(treason, aggregate.executed_tasks);
    rmi.broadcast(reason_and_task, rmi.procid() == 0);
    termination_reason = exec_status(reason_and_task.first);
    return reason_and_task.second;
  }
 
  void start_thread(size_t threadid) {
    // create the scope
    dgraph_context<engine_type> context(&graph, this, threadid);
    timer ti;

    // loop over iterations
    size_t iter = 0;
    bool usestatic = max_iterations > 0;
    UpdateFunctor fn;
    while(1) {
      // if max_iterations is defined, quit
      if (usestatic && iter >= max_iterations) {
        termination_reason = EXEC_TASK_DEPLETION;
        break;
      }
      bool hassynctasks = active_sync_tasks.size() > 0;
      // loop over colors    
      for (size_t c = 0;c < color_block.size(); ++c) {
        // internal loop over vertices in the color
        while(1) {
          // grab a vertex  
          size_t i = curidx.inc_ret_last();  
          // if index out of scope, we are done with this color. break
          if (i >= color_block[c].size()) break;
          // otherwise, get the local and globalvid
          vertex_id_t localvid = color_block[c][i];
          vertex_id_t globalvid = graph.localvid_to_globalvid(color_block[c][i]);
          
          
          if (usestatic || vfun_set.test_and_get(localvid, fn)) {
            if (!usestatic) num_pending_tasks.dec();
            context.init(globalvid, consistency_model::NULL_CONSISTENCY);
            // call the update functor
            if(fn.is_factorizable()) {
              if(ufun.gather_edges() & update_functor_type::OUT_EDGES) {  
                const edge_list_type edges = graph.out_edge_ids(vid);
                foreach(const edge_id_type eid, edges) {
                  ufun.gather(context, eid);
                }
              }
              if(ufun.gather_edges() & update_functor_type::IN_EDGES) {  
                const edge_list_type edges = graph.in_edge_ids(vid);
                foreach(const edge_id_type eid, edges) {
                  ufun.gather(context, eid);
                } 
              }
              ufun.apply(context);
              if(ufun.gather_edges() & update_functor_type::OUT_EDGES) {  
                const edge_list_type edges = graph.out_edge_ids(vid);
                foreach(const edge_id_type eid, edges) {
                  ufun.scatter(context, eid);
                }
              }
              if(ufun.gather_edges() & update_functor_type::IN_EDGES) {  
                const edge_list_type edges = graph.in_edge_ids(vid);
                foreach(const edge_id_type eid, edges) {
                  ufun.scatter(context, eid);
                } 
              }
            } else { 
              ufun(context);
            }
            if (hassynctasks) eval_syncs(globalvid, context, threadid);
            scope.commit_async_untracked();
            update_counts[threadid]++;
          }
          else {
            // ok this vertex is not scheduled. But if there are syncs
            // to run I will still need to get the scope
            context.init(globalvid, consistency_model::NULL_CONSISTENCY);
            if (hassynctasks) eval_syncs(globalvid, context, threadid);
            scope.commit_async_untracked();
          }
        }
        // wait for all threads to synchronize on this color.
        thread_color_barrier.wait();
        curidx.value = 0;
        // full barrier on the color
        // this will complete synchronization of all add tasks as well
        if (threadid == 0) {
          ti.start();
          graph.wait_for_all_async_syncs();
          // TODO! If synchronize() calls were made then this barrier is necessary
          // but the time needed to figure out if a synchronize call is required 
          // could be as long as the barrier itself
          if (const_nbr_vertices == false || const_edges == false)  rmi.dc().barrier();
          rmi.dc().full_barrier();
          num_dist_barriers_called++;
          //std::cout << rmi.procid() << ": Full Barrier at end of color" << std::endl;
          barrier_time += ti.current_time();
        }
        thread_color_barrier.wait();

      }


      sync_end_iteration(threadid);
      thread_color_barrier.wait();
      if (threadid == 0) {
        ti.start();
        //std::cout << rmi.procid() << ": End of all colors" << std::endl;
        size_t numtasksdone = check_global_termination(!usestatic);

        //std::cout << numtasksdone << " tasks done" << std::endl;
        compute_sync_schedule(numtasksdone);
        barrier_time += ti.current_time();
      }
      // all threads must wait for 0
      thread_color_barrier.wait();

      ++iter;
      if (termination_reason != EXEC_UNSET) {
        //std::cout << rmi.procid() << ": Termination Reason: " << termination_reason << std::endl;
        break;
      }
    }
  }
  
  void set_const_edges(bool const_edges_ = true) {
    const_edges = const_edges_;
  }

  void set_const_nbr_vertices(bool const_nbr_vertices_ = true) {
    const_nbr_vertices = const_nbr_vertices_;
  }

  
  /** Execute the engine */
  void start() {
    assert(update_function != NULL);
    
    if (default_scope_range == scope_range::FULL_CONSISTENCY) {
      const_nbr_vertices = false;
    }
    // generate colors then
    // wait for everyone to enter start    
    generate_color_blocks();
    init_syncs();
    termination_reason = EXEC_UNSET;
    barrier_time = 0.0;
    force_stop = false;
    numsyncs.value = 0;
    num_dist_barriers_called = 0;
    std::fill(update_counts.begin(), update_counts.end(), 0);
    // two full barrers to complete flush replies
    rmi.dc().full_barrier();
    rmi.dc().full_barrier();

    // reset indices
    curidx.value = 0;
    ti.start();
    // spawn threads
    thread_group thrgrp; 
    for (size_t i = 0;i < ncpus; ++i) {
      size_t aff = use_cpu_affinity ? i : -1;
      thrgrp.launch(boost::bind(
                            &distributed_chromatic_engine<Graph>::start_thread,
                            this, i), aff);
    }
    
    thrgrp.join();              
    rmi.barrier();
    

    
    // proc 0 gathers all update counts
    std::vector<size_t> procupdatecounts(rmi.numprocs(), 0);
    procupdatecounts[rmi.procid()] = thisproc_update_counts();
    rmi.gather(procupdatecounts, 0);
    
    std::vector<double> barrier_times(rmi.numprocs(), 0);
    barrier_times[rmi.procid()] = barrier_time;
    rmi.gather(barrier_times, 0);
    // get RMI statistics
    std::map<std::string, size_t> ret = rmi.gather_statistics();

    if (rmi.procid() == 0) {
      engine_metrics.add("runtime",
                        ti.current_time(), TIME);
      total_update_count = 0;
      for(size_t i = 0; i < procupdatecounts.size(); ++i) {
        engine_metrics.add_vector_entry("updatecount", i, procupdatecounts[i]);
        total_update_count +=  procupdatecounts[i];
      }
      total_barrier_time = 0;
      for(size_t i = 0; i < barrier_times.size(); ++i) {
        engine_metrics.add_vector_entry("barrier_time", i, barrier_times[i]);
        total_barrier_time += barrier_times[i];
      }

      engine_metrics.set("termination_reason", 
                        exec_status_as_string(termination_reason));
      engine_metrics.add("dist_barriers_issued",
                        num_dist_barriers_called, INTEGER);

      engine_metrics.set("num_vertices", graph.num_vertices(), INTEGER);
      engine_metrics.set("num_edges", graph.num_edges(), INTEGER);
      engine_metrics.add("num_syncs", numsyncs.value, INTEGER);
      engine_metrics.set("isdynamic", max_iterations == 0, INTEGER);
      engine_metrics.add("iterations", max_iterations, INTEGER);
      engine_metrics.set("total_calls_sent", ret["total_calls_sent"], INTEGER);
      engine_metrics.set("total_bytes_sent", ret["total_bytes_sent"], INTEGER);
      total_bytes_sent = ret["total_bytes_sent"];
    }
    
    
    
  }
  
  
    /** \brief Update the scheduler options.  */
  void set_options(const graphlab_options& opts) {
    opts.get_engine_args().get_int_option("max_iterations", max_iterations);
    opts.get_engine_args().get_int_option("randomize_schedule", randomize_schedule);
    rmi.barrier();
  }

  void set_randomize_schedule(bool randomize_schedule_) {
    randomize_schedule = randomize_schedule_;
    rmi.barrier();
  }

  
  
  static void print_options_help(std::ostream &out) {
    out << "max_iterations = [integer, default = 0]\n";
    out << "randomize_schedule = [integer, default = 0]\n";
  };


  void stop() {
    force_stop = true;
  }

  
  size_t total_update_count;
  
  size_t get_tasks_done() const {
    return  total_update_count;
  }
  
  double total_barrier_time;
  double get_barrier_time() const {
      return total_barrier_time;
  }
  
  long long int total_bytes_sent;
  long long int get_total_bytes_sent() {
     return total_bytes_sent;
  }

  metrics get_metrics() {
    return engine_metrics;
  }


  void reset_metrics() {
    engine_metrics.clear();
  }

};

} // namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // DISTRIBUTED_CHROMATIC_ENGINE_HPP


