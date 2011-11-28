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


#ifndef DISTRIBUTED2_CHROMATIC_ENGINE_HPP
#define DISTRIBUTED2_CHROMATIC_ENGINE_HPP
#include <functional>
#include <algorithm>
#include <ext/functional> // for select1st
#include <boost/bind.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/engine/execution_status.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/distributed2/graph/dgraph_context.hpp>
#include <graphlab/distributed2/distributed_shared_data.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {



/**
All processes must receive the same options at the same time.
i.e. if set_cpu_affinities is called, all processes mus call it at the same time.
This is true for all set_* functions.
*/
template <typename Graph, typename UpdateFunctor>
class distributed_chromatic_engine : public iengine<Graph, UpdateFunctor> {
 public:
  // Include parent types
  typedef iengine<Graph, UpdateFunctor> iengine_base;
  typedef distributed_chromatic_engine<Graph, UpdateFunctor> engine_type;
  typedef typename iengine_base::graph_type graph_type;
  typedef typename iengine_base::update_functor_type update_functor_type;
  
  typedef typename graph_type::vertex_data_type vertex_data_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_id_type   edge_id_type;
  typedef typename graph_type::edge_list_type edge_list_type;

  typedef typename iengine_base::icontext_type  icontext_type;
  typedef dgraph_context<distributed_chromatic_engine>  context_type;
  typedef distributed_shared_data<context_type> shared_data_type;
  
  typedef typename iengine_base::termination_function_type termination_function_type;


 private:
  // the local rmi instance
  dc_dist_object<distributed_chromatic_engine<Graph, UpdateFunctor> > rmi;
  
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
  
  shared_data_type shared_data;
  
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
  execution_status::status_enum termination_reason;

  consistency_model::model_enum default_scope_range;
 
  std::vector<std::vector<vertex_id_type> > color_block; // set of localvids in each color
  vertex_functor_set<UpdateFunctor> vfunset;
  vertex_functor_set<UpdateFunctor> vfun_cacheset;
  
  graphlab_options opts;
  
  size_t max_iterations;
  double barrier_time;
  size_t num_dist_barriers_called;
  
  // other optimizations
  bool const_nbr_vertices, const_edges;
  struct sync_task {
    std::string key;
    size_t sync_interval;
    size_t next_time;
    vertex_id_type rangelow;
    vertex_id_type rangehigh;
  };
  
  /// A list of all registered sync tasks
  std::vector<sync_task> sync_tasks;
  
  /// The list of tasks which are currently being evaluated
  std::vector<sync_task*> active_sync_tasks;

  
  
 public:
  distributed_chromatic_engine(distributed_control &dc,
                                    Graph& graph,
                                    size_t ncpus = 1):
                            rmi(dc, this),
                            graph(graph),
                            ncpus( std::max(ncpus, size_t(1)) ),
                            use_cpu_affinity(false),
                            use_sched_yield(true),
                            update_counts(std::max(ncpus, size_t(1)), 0),
                            shared_data(dc, ncpus),
                            timeout_millis(0),
                            force_stop(false),
                            task_budget(0),
                            randomize_schedule(0),
                            termination_reason(execution_status::UNSET),
                            vfunset(graph.owned_vertices().size()),
                            vfun_cacheset(graph.get_local_store().num_vertices()),
                            max_iterations(0),
                            barrier_time(0.0),
                            const_nbr_vertices(true),
                            const_edges(false),
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
  void set_default_scope(consistency_model::model_enum default_scope_range_) {
    default_scope_range = default_scope_range_;
    rmi.barrier();
  }
  
  

  /**
   * Return the reason why the engine last terminated
   */
  execution_status::status_enum last_exec_status() const {
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
  void add_termination_condition(termination_function_type term) {
    term_functions.push_back(term);
    rmi.barrier();
  }


  /**
   * Clear all terminators from the engine
   * Must be called by all machines simultaneously.
   */
   void clear_termination_conditions() {
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
  
  void schedule_from_context(const vertex_id_type vid, 
                             const update_functor_type& fun) {
    vfun_cacheset.add(graph.globalvid_to_localvid(vid), fun);
  }
  
  
  void synchronize_vfun_cache() {
    // commit local changes
#pragma omp parallel for
    for (ssize_t i = 0;i < (ssize_t)graph.local_vertices(); ++i) {
      update_functor_type uf;
      // once again, local VIDs are sequential and at the start
      bool hasf = vfun_cacheset.test_and_get((vertex_id_type)(i), uf);
      if (hasf && vfunset.add((vertex_id_type)(i), uf)) {
        num_pending_tasks.inc();
      }      
    }
    
    std::vector<
        std::vector<
            std::pair<vertex_id_type, update_functor_type> > > remote_schedules;

    remote_schedules.resize(rmi.numprocs());
    std::vector<mutex> schedulelock(rmi.numprocs());
    // commit remote changes
#pragma omp parallel for
    for (ssize_t i = graph.local_vertices();
         i < (ssize_t)graph.get_local_store().num_vertices(); ++i) {
      update_functor_type uf;
      bool hasf = vfun_cacheset.test_and_get((vertex_id_type)(i), uf);
      if (hasf) {
        procid_t owner = graph.localvid_to_owner((vertex_id_type)(i));
        schedulelock[owner].lock();
        remote_schedules[owner].push_back(std::make_pair(graph.localvid_to_globalvid((vertex_id_type)(i)), 
                                                          uf));
        schedulelock[owner].unlock();
      }
    }
    
    // transmit
    for (procid_t i = 0;i < rmi.numprocs(); ++i) {
      if (i != rmi.procid() && remote_schedules[i].size() > 0) {
        rmi.remote_call(i, 
            &distributed_chromatic_engine<Graph, UpdateFunctor>::
                                                          schedule_collection,
            remote_schedules[i]);
      }
    }
  }

  void schedule_impl(const vertex_id_type vid, 
                const update_functor_type& fun) {
    if (graph.is_owned(vid)) {
      if (vfunset.add(graph.globalvid_to_localvid(vid), fun)) {
        num_pending_tasks.inc();
      }
    }
    else {
      rmi.remote_call(graph.globalvid_to_owner(vid),
                &distributed_chromatic_engine<Graph, UpdateFunctor>::schedule_impl,
                vid, fun);
    }
  }


  /**
   * \brief Adds an update task with a particular priority.
   * add_task on any vertex can be called by any machine.
   * The call is asynchronous and may not be completed until
   * a full_barrier is issued.
   */
  void schedule(const vertex_id_type vid, 
                const update_functor_type& fun) {
    schedule_impl(vid, fun);
  }

  /**
   * \brief Creates a collection of tasks on all the vertices in
   * 'vertices', and all with the same update function and priority
   * This function is forwarded to the scheduler.
   */
  void schedule(const std::vector<vertex_id_type>& vertices,
                const update_functor_type& update_functor) {
    // not the most efficient way to do it...
    for (size_t i = 0;i < vertices.size(); ++i) {
      schedule(vertices[i], update_functor);
    }
  }

  void schedule_collection(const std::vector<
                                     std::pair<vertex_id_type, 
                                             update_functor_type> >& tasks) {
    for (size_t i = 0;i < tasks.size(); ++i) {
      schedule(tasks[i].first, tasks[i].second);
    }
  }
  
  /**
   * \brief Creates a collection of tasks on all the vertices in the graph,
   * with the same update function and priority
   * Must be called by all machines simultaneously
   */
  void schedule_all(const update_functor_type& update_functor) { 
    schedule(graph.owned_vertices(), update_functor);
    rmi.barrier();
  }

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
    shared_data.add_global(key, value, size);
  }

  /**
    * Define a global constant.
    */
  template< typename T >
  void add_global_const(const std::string& key, const T& value, 
                        size_t size = 1) {
    shared_data.add_global_const(key, value, size);
  }


  //! Change the value of a global entry
  template< typename T >
  void set_global(const std::string& key, const T& value, 
                  size_t index = 0) {
    shared_data.set_global(key, value, index);
    shared_data.synchronize_global(key, value, index);
  }

  //! Get a copy of the value of a global entry
  template< typename T >
  T get_global(const std::string& key, size_t index = 0) const {
    return shared_data.get_global<T>(key, index);
  }



  /**
    Registers a sync operation.
    Must be called by all machine simultaneously
  */
    template<typename Accum>
    void add_sync(const std::string& key,            
                  const Accum& zero,                 
                  size_t sync_interval,
                  bool use_barrier = false,
                  vertex_id_type begin_vid = 0,
                  vertex_id_type end_vid = 
                  std::numeric_limits<vertex_id_type>::max()) {
    sync_task st;
    st.key = key;
    st.sync_interval = sync_interval;
    st.next_time = 0;
    st.rangelow = begin_vid;
    st.rangehigh = end_vid;
    sync_tasks.push_back(st);
    shared_data.add_sync(key, zero);
    rmi.barrier();
  }

  
  void generate_color_blocks() {
    // construct for each color, the set of vertices as well as the 
    // number of replicas for that vertex.
    // the number of replicas - 1 is the amount of communication
    // we have to perform to synchronize modifications to that vertex


    std::vector<std::vector<std::pair<size_t, vertex_id_type> > > color_block_and_weight;
    const size_t num_colors(graph.recompute_num_colors());
    // the list of vertices for each color
    color_block_and_weight.resize(num_colors);
    
    foreach(vertex_id_type v, graph.owned_vertices()) {
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
                    __gnu_cxx::select2nd<std::pair<size_t, vertex_id_type> >());
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
   * Initialize the sync tasks. Called by start()
   */
  void init_syncs() {
     active_sync_tasks.clear();
     shared_data.reset_all_syncs();
     for (size_t i = 0;i < sync_tasks.size(); ++i) {
       // everyone runs at the start even if scheduling interval is 0
       active_sync_tasks.push_back(&(sync_tasks[i]));
     }
  }

  /**
   * Called whenever a vertex is executed.
   * Accumulates the available syncs
   */
  void eval_syncs(vertex_id_type curvertex, 
                  context_type& context, 
                  size_t threadid) {
     // go through all the active sync tasks
     foreach(sync_task* task, active_sync_tasks) {
       // if in range, sync!
       if (task->rangelow <= curvertex && curvertex <= task->rangehigh) {
         shared_data.accumulate(task->key, &context, threadid);
       }
     }
  }

  /** Called at the end of the iteration. Called by all threads after a barrier*/
  void sync_end_iteration(size_t threadid, context_type& globalcontext) {
    if (threadid == 0) {
      foreach(sync_task* task, active_sync_tasks) {
        shared_data.finalize(task->key, &globalcontext);
      }
      shared_data.wait_for_all_communication();
      shared_data.reset_all_syncs();
    }
    thread_color_barrier.wait();
  }

  /** clears the active sync tasks and figure out what syncs to run next.
      Called by one thread from each machine after sync_end_iteration */
  void compute_sync_schedule(size_t num_executed_tasks) {
    // update the next time variable
    for (size_t i = 0;i < active_sync_tasks.size(); ++i) {
      sync_tasks[i].next_time = num_executed_tasks + sync_tasks[i].sync_interval;
      //  if it is within the next iteration 
      if (sync_tasks[i].next_time < num_executed_tasks + graph.num_vertices()) sync_tasks[i].next_time = num_executed_tasks;
      // if sync interval of 0, this was the first iteration.
      // then I just set next time to infinity and it will never be run again
      if (sync_tasks[i].sync_interval == 0) {
        sync_tasks[i].next_time = size_t(-1);
      }
    }
    active_sync_tasks.clear();
    // figure out what to run next
    for (size_t i = 0;i < sync_tasks.size(); ++i) {
      if (sync_tasks[i].next_time <= num_executed_tasks) {
        active_sync_tasks.push_back(&(sync_tasks[i]));
      }
    }
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
        termination_reason = execution_status::TASK_DEPLETION;
      }
      else if (task_budget > 0 && aggregate.executed_tasks >= task_budget) {
        termination_reason = execution_status::TASK_BUDGET_EXCEEDED;
      }
      else if (timeout_millis > 0 && aggregate.timeout) {
        termination_reason = execution_status::TIMEOUT;
      }
      else if (aggregate.terminator) {
        termination_reason = execution_status::TERM_FUNCTION;
      }
      else if (aggregate.force_stop) {
        termination_reason = execution_status::FORCED_ABORT;
      }
    }
    size_t treason = termination_reason;
    // note this is OK because only machine 0 will have the right value for
    // executed_tasks. And everyone is receiving from machine 0
    std::pair<size_t, size_t> reason_and_task(treason, aggregate.executed_tasks);
    rmi.broadcast(reason_and_task, rmi.procid() == 0);
    termination_reason = (execution_status::status_enum)(reason_and_task.first);
    return reason_and_task.second;
  }
 
  void start_thread(size_t threadid) {
    // create the scope
    dgraph_context<engine_type> context(this, &graph, &shared_data);
    timer ti;

    // loop over iterations
    size_t iter = 0;
    bool usestatic = max_iterations > 0;
    while(1) {
      // if max_iterations is defined, quit
      if (usestatic && iter >= max_iterations) {
        termination_reason = execution_status::TASK_DEPLETION;
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
          vertex_id_type localvid = color_block[c][i];
          vertex_id_type globalvid = graph.localvid_to_globalvid(color_block[c][i]);
          
          update_functor_type functor_to_run;
          bool has_functor_to_run = false;
          if (usestatic) {
              has_functor_to_run = vfunset.read_value(localvid, functor_to_run);
          }
          else {
            has_functor_to_run = vfunset.test_and_get(localvid, functor_to_run);
          }
          
          if (has_functor_to_run) {
            if (!usestatic) num_pending_tasks.dec();
            // otherwise. run the vertex
            // create the scope
            context.init(globalvid);
            // run the update function
            functor_to_run(context);
            // check if there are tasks to run
            if (hassynctasks) eval_syncs(globalvid, context, threadid);
            context.commit_async_untracked();
            update_counts[threadid]++;
          }
          else {
            // ok this vertex is not scheduled. But if there are syncs
            // to run I will still need to get the scope
            if (hassynctasks) {
              context.init(globalvid);
              eval_syncs(globalvid, context, threadid);
              context.commit_async_untracked();
            }
          }
        }
        // wait for all threads to synchronize on this color.
        thread_color_barrier.wait();
        curidx.value = 0;
        // full barrier on the color
        // this will complete synchronization of all add tasks as well
        if (threadid == 0) {
          ti.start();
          synchronize_vfun_cache();          
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


      sync_end_iteration(threadid, context);
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
      if (termination_reason != execution_status::UNSET) {
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
    if (default_scope_range == consistency_model::FULL_CONSISTENCY) {
      const_nbr_vertices = false;
    }
    // generate colors then
    // wait for everyone to enter start    
    generate_color_blocks();
    vfun_cacheset.clear_unsync();
    init_syncs();
    termination_reason = execution_status::UNSET;
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
                            &distributed_chromatic_engine<Graph, UpdateFunctor>::start_thread,
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

//     if (rmi.procid() == 0) {
//       engine_metrics.add("runtime",
//                         ti.current_time(), TIME);
//       total_update_count = 0;
//       for(size_t i = 0; i < procupdatecounts.size(); ++i) {
//         engine_metrics.add_vector_entry("updatecount", i, procupdatecounts[i]);
//         total_update_count +=  procupdatecounts[i];
//       }
//       total_barrier_time = 0;
//       for(size_t i = 0; i < barrier_times.size(); ++i) {
//         engine_metrics.add_vector_entry("barrier_time", i, barrier_times[i]);
//         total_barrier_time += barrier_times[i];
//       }
// 
//       engine_metrics.set("termination_reason", 
//                         exec_status_as_string(termination_reason));
//       engine_metrics.add("dist_barriers_issued",
//                         num_dist_barriers_called, INTEGER);
// 
//       engine_metrics.set("num_vertices", graph.num_vertices(), INTEGER);
//       engine_metrics.set("num_edges", graph.num_edges(), INTEGER);
//       engine_metrics.add("num_syncs", numsyncs.value, INTEGER);
//       engine_metrics.set("isdynamic", max_iterations == 0, INTEGER);
//       engine_metrics.add("iterations", max_iterations, INTEGER);
//       engine_metrics.set("total_calls_sent", ret["total_calls_sent"], INTEGER);
//       engine_metrics.set("total_bytes_sent", ret["total_bytes_sent"], INTEGER);
//       total_bytes_sent = ret["total_bytes_sent"];
//     }
//     
//     
    
  }
  
  /**
   * Performs a sync immediately. This function requires that the shared
   * variable already be registered with the engine.
   * and that the engine is not currently running
   * All processes must call simultaneously
   */
  void sync_now(const std::string& key) {
    // TODO
  }
  
  /**
    * \brief Update the engine options.  
    *
    * Setting the engine options will cause all existing state,
    * including scheduled update functors, to be cleared.
    */
  void set_options(const graphlab_options& newopts) {
    opts = newopts;
    opts.engine_args.get_option("max_iterations", max_iterations);
    opts.engine_args.get_option("randomize_schedule", randomize_schedule);
    rmi.barrier();
  }

  //! \brief Get the current engine options for this engine
  const graphlab_options& get_options() { return opts; } 
  
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
};

} // namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // DISTRIBUTED_CHROMATIC_ENGINE_HPP

