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
  typedef Graph graph_type;
  typedef typename iengine_base::update_functor_type update_functor_type;
  
  typedef typename graph_type::vertex_data_type vertex_data_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_id_type   edge_id_type;
  typedef typename graph_type::edge_list_type edge_list_type;
  typedef typename graph_type::edge_type edge_type;  
  
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
  float start_time;
  timer ti;

  
  /// Used to identify when the engine is stopped through stop
  bool force_stop;
  
  /** The total number of tasks that should be executed */
  size_t task_budget;
  
  size_t randomize_schedule;
  
  /** If dynamic scheduling is used, the number of scheduled tasks */
  atomic<size_t> num_pending_tasks;
  atomic<size_t> num_cached_tasks;
  
  
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
  bool use_factorized;
  bool no_graph_synchronization;
  bool no_colors;
  size_t factor_threshold;
  
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

  std::vector<mutex> synclock;
  
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
                            start_time(0),
                            force_stop(false),
                            task_budget(0),
                            randomize_schedule(0),
                            termination_reason(execution_status::UNSET),
                            vfunset(graph.owned_vertices().size()),
                            vfun_cacheset(graph.get_local_store().num_vertices()),
                            max_iterations(0),
                            barrier_time(0.0),
                            use_factorized(false),
                            no_graph_synchronization(false),
                            no_colors(false),
                            factor_threshold(1000),
                            const_nbr_vertices(true),
                            const_edges(false),
                            synclock(ncpus),
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
  
  void schedule_from_context(size_t threadid,
                             const vertex_id_type vid, 
                             const update_functor_type& fun) {
    if (vfun_cacheset.add(graph.globalvid_to_localvid(vid), fun)) {
      num_cached_tasks.inc();
    }
    
  }
  
  
   

  
  void synchronize_vfun_cache(std::vector<vertex_id_type> &nextschedule) {
    std::vector<
        std::vector<
            std::pair<vertex_id_type, update_functor_type> > > remote_schedules;

    remote_schedules.resize(rmi.numprocs());
    std::vector<mutex> schedulelock(rmi.numprocs());
    
#pragma omp parallel for
    for (ssize_t idx = 0; idx < nextschedule.size(); ++idx) {
      // commit local changes
      vertex_id_type i = nextschedule[idx];
      
      update_functor_type uf;
      bool hasf = vfun_cacheset.test_and_get(i, uf);
      if (hasf) {
        num_cached_tasks.dec();
        if (i < graph.local_vertices()) {
          if (vfunset.add((vertex_id_type)(i), uf)) {
            num_pending_tasks.inc();
          }      
        }
        else {
          procid_t owner = graph.localvid_to_owner((vertex_id_type)(i));
          schedulelock[owner].lock();
          remote_schedules[owner].push_back(std::make_pair(graph.localvid_to_globalvid(i), 
                                                            uf));
          schedulelock[owner].unlock();
        }
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
    // color blocks include ghosts


    std::vector<std::vector<std::pair<size_t, vertex_id_type> > > color_block_and_weight;
    const size_t num_colors(no_colors ? 1 : graph.recompute_num_colors());
    if (no_colors) {
      color_block_and_weight.resize(num_colors);
      foreach(vertex_id_type v, graph.owned_vertices()) {
        color_block_and_weight[0].push_back(
                                      std::make_pair(graph.globalvid_to_replicas(v).size(), 
                                                graph.globalvid_to_localvid(v)));
      }

    foreach(vertex_id_type v, graph.ghost_vertices()) {
        color_block_and_weight[0].push_back(
                                      std::make_pair(std::numeric_limits<size_t>::max(), 
                                                graph.globalvid_to_localvid(v)));
      }
    }
    else {
      // the list of vertices for each color
      color_block_and_weight.resize(num_colors);
      foreach(vertex_id_type v, graph.owned_vertices()) {
        color_block_and_weight[graph.get_color(v)].push_back(
                                      std::make_pair(graph.globalvid_to_replicas(v).size(), 
                                                graph.globalvid_to_localvid(v)));
      }

      foreach(vertex_id_type v, graph.ghost_vertices()) {
        color_block_and_weight[graph.get_color(v)].push_back(
                                      std::make_pair(std::numeric_limits<size_t>::max(), 
                                                graph.globalvid_to_localvid(v)));
      }
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
    std::cout << color_block_and_weight.size() << " colors\n";
    // insert the sorted vertices into the final color_block
    for (size_t i = 0;i < color_block_and_weight.size(); ++i ) {  
      std::cout << color_block_and_weight[i].size()<< " ";
      std::transform(color_block_and_weight[i].begin(),
                    color_block_and_weight[i].end(), 
                    std::back_inserter(color_block[i]),
                    __gnu_cxx::select2nd<std::pair<size_t, vertex_id_type> >());
    }
    std::cout << "\n";
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
    synclock[threadid].lock();
     // go through all the active sync tasks
     foreach(sync_task* task, active_sync_tasks) {
       // if in range, sync!
       if (task->rangelow <= curvertex && curvertex <= task->rangehigh) {
         shared_data.accumulate(task->key, &context, threadid);
       }
     }
     synclock[threadid].unlock();
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
      termination_test[rmi.procid()].pending_tasks = num_pending_tasks.value + num_cached_tasks.value;
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
          if (color_block[c][i] >= graph.local_vertices()) continue;
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
            context.init(globalvid, threadid);
            // run the update function
            functor_to_run(context);
            // check if there are tasks to run
            if (hassynctasks) eval_syncs(globalvid, context, threadid);
            if (no_graph_synchronization == false) context.commit_async_untracked();
            update_counts[threadid]++;
          }
          else {
            // ok this vertex is not scheduled. But if there are syncs
            // to run I will still need to get the scope
            if (hassynctasks) {
              context.init(globalvid, threadid);
              eval_syncs(globalvid, context, threadid);
              if (no_graph_synchronization == false) context.commit_async_untracked();
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
          synchronize_vfun_cache(color_block[(c + 1) % color_block.size()]);

          if (no_graph_synchronization == false) {
            graph.wait_for_all_async_syncs();
          }
          rmi.dc().barrier();
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
  
/////////////// Begin Factorized Implementation ////////////////////////
  
  struct gather_receive{
    update_functor_type accumulated_functor;
    mutex lock;
    size_t ctr;
  };
  
  boost::unordered_map<vertex_id_type, gather_receive> localvid_to_gather;

  
  dense_bitset seperator_set;
  size_t pending_factorized_updates;
  mutex pending_factorized_updates_lock;
  conditional pending_factorized_updates_cond;
  
  
  void compute_threshold_seperator_set() {
    std::vector<std::vector<vertex_id_type> > allseperatorsets(rmi.numprocs());
    
    // the seperator set is all bonudary vertices with degree > threshold
    for (vertex_id_type localvid = 0; 
        localvid < graph.local_vertices(); 
        ++localvid) {
      vertex_id_type globalvid = graph.localvid_to_globalvid(localvid);
      if (graph.on_boundary(globalvid)) {
        size_t numin = graph.get_local_store().num_in_neighbors(localvid);
        size_t numout = graph.get_local_store().num_out_neighbors(localvid);
        if (numin + numout >= factor_threshold && random::bernoulli()) {
          allseperatorsets[rmi.procid()].push_back(globalvid);
        }
      }
    }
    
    // not exactly the most efficient way to do this. really, we only
    // need to collect all the vertex IDs that will land on this machine
    rmi.all_gather(allseperatorsets);
    
    seperator_set.resize(graph.get_local_store().num_vertices());

    seperator_set.clear();
    
    // build the local seperator set 
    for (size_t i = 0;i < allseperatorsets.size(); ++i) {
      for (size_t j = 0; j < allseperatorsets[i].size(); ++j) {
        if (graph.vertex_is_local(allseperatorsets[i][j])) {
          vertex_id_type localvid = graph.globalvid_to_localvid(allseperatorsets[i][j]);
          seperator_set.set_bit(localvid);
        }
      }
    }
    if (rmi.procid() == 0) {
      // sum and print the number of seperator vertices
      size_t numsep = 0;
      for (size_t i = 0;i < allseperatorsets.size(); ++i) numsep += allseperatorsets[i].size();
      
      logstream(LOG_INFO) << "Seperator Set Size: " << numsep << std::endl;
    }
  }
  /* // used to debug seperators for the demo app
  void draw_seperator_set() {
    // try to draw the seperator set
    size_t v = 0;
    for (size_t i = 0;i < 20; ++i) {
      for (size_t j = 0;j < 20; ++j) {
        if (graph.vertex_is_local(v)) {
          if (graph.is_owned(v)) {
            if (seperator_set.get(graph.globalvid_to_localvid(v))) {
              textcolor(stdout, 1, 1);  //RED
              std::cout << rmi.procid() << " ";
              reset_color(stdout);
            }
            else {
              std::cout << rmi.procid() << " ";
            }
          }
          else {
            if (seperator_set.get(graph.globalvid_to_localvid(v))) {
              textcolor(stdout, 1, 1);  //RED
              std::cout << "*" << " ";
              reset_color(stdout);
            }
            else {
              std::cout << "*" << " ";
            }
          }
        }
        else {
          std::cout << "  ";
        }
        ++v;
      }
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }*/
  /**
   * 
   * Factorized update functors are somewhat awkward because the graph
   * is really an edge seperated graph. So I need to coerce the system
   * to behave like vertex seperators.
   * 
   * 1: A subset of vertices (vertices with degree > factor_threshold)
   *    are picked as vertex seperators.
   * 2: These subset of vertices use factorized updates, everything else 
   *    use regular updates.
   * 3: All synchronization rules hold just like regular updates with the exception that
   *    - if a LOCAL vertex v is
   *       a) adjacent to ghost vertices
   *       b) all adjacent ghost vertices are vertex seperators
   *    - then, no synchronization is necessary
   * 
   */
  void infer_vertex_seperator_set() {
    if (factor_threshold > 0) {
      compute_threshold_seperator_set();
    }
    else {
      graph.greedy_vertex_seperator(seperator_set);
    }
/*    rmi.barrier();
    for (procid_t p = 0;p < rmi.numprocs(); ++p) {
      if (p == rmi.procid()) draw_seperator_set();
      rmi.barrier();
    }
    if (rmi.procid() == 0) getchar();
    rmi.barrier(); */
    logstream(LOG_INFO) << "Rebuilding Replica Set" << std::endl;
    graph.rebuild_replica_set_assuming_vertex_seperator(seperator_set);
  }

  /**
   * proceeds to plan a factorized schedule for a subset of the vertices
   */
  size_t plan_factorized_schedule(std::vector<vertex_id_type> &next_localvset) {
    std::vector<
        std::vector<
            std::pair<vertex_id_type, update_functor_type> > > remote_schedules;

    remote_schedules.resize(rmi.numprocs());

    // look through the localvset and remotely 
    // schedule those which are vertex seperators
    size_t num_owned_factorized_updates = 0;
    foreach(vertex_id_type localvid, next_localvset) {
      if (localvid < graph.local_vertices() && seperator_set.get(localvid)) {
        update_functor_type uf;
        bool hasf = vfunset.read_value(localvid, uf);
        if (hasf) {
          ++num_owned_factorized_updates;
          vertex_id_type globalvid = graph.localvid_to_globalvid(localvid);
/*            std::cout << rmi.procid() << ": Planning : " << globalvid 
                    << ". color = " << graph.color(globalvid) << " " 
                      << graph.localvid_to_replicas(localvid).popcount() << " nodes" << std::endl;*/
          foreach(procid_t rep, graph.localvid_to_replicas(localvid)) {
            if (rep != rmi.procid()) {
              remote_schedules[rep].push_back(std::make_pair(globalvid, uf));
            }
          }
        }
      }
    }
    
    // transmit
    for (procid_t i = 0;i < rmi.numprocs(); ++i) {
      if (i != rmi.procid() && remote_schedules[i].size() > 0) {
        rmi.remote_call(i, 
            &distributed_chromatic_engine<Graph, UpdateFunctor>::
                                                          factorized_schedule_collection,
            remote_schedules[i]);
      }
    }
    return num_owned_factorized_updates;
  }
  
  
  void factorized_schedule_collection(const std::vector<
                                     std::pair<vertex_id_type, 
                                             update_functor_type> >& tasks) {
    for (size_t i = 0;i < tasks.size(); ++i) {
/*      std::cout << rmi.procid() << ": Adding fact. task: " 
                << tasks[i].first << " color = " << graph.color(tasks[i].first) << std::endl;*/
      
      if (vfunset.add(graph.globalvid_to_localvid(tasks[i].first), tasks[i].second)) {
        num_pending_tasks.inc();
      }
    }
  }

  
  
  
  
  
  
  
  
  
  void reset_gather_receive_structures() {
    foreach(uint32_t i, seperator_set) {
      if (i < graph.local_vertices()) {
        localvid_to_gather[i].ctr = graph.localvid_to_replicas(i).popcount();
        localvid_to_gather[i].accumulated_functor = update_functor_type();
      }
    }
  }
  
  void factorized_gather_receive(vertex_id_type globalvid, const update_functor_type &ufun) {
    vertex_id_type localvid = graph.globalvid_to_localvid(globalvid);
    ASSERT_LT(localvid, graph.local_vertices());
    typename boost::unordered_map<vertex_id_type, gather_receive>::iterator iter = 
                          localvid_to_gather.find(localvid);
    iter->second.lock.lock();
    iter->second.accumulated_functor.merge(ufun);
    iter->second.ctr--;
    bool done = (iter->second.ctr == 0);
//    std::cout << rmi.procid() << " gather receive on " << globalvid << 
//                ", " << iter->second.ctr << " remaining" << std::endl;
    iter->second.lock.unlock();
    if (done) {
      apply_and_broadcast_scatter(localvid);
    }
  }
  
  void apply_and_broadcast_scatter(vertex_id_type localvid) {
    typename boost::unordered_map<vertex_id_type, gather_receive>::iterator iter = 
                      localvid_to_gather.find(localvid);
    dgraph_context<engine_type> context(this, &graph, &shared_data);
    vertex_id_type globalvid = graph.localvid_to_globalvid(localvid);
    context.init(globalvid, localvid % ncpus);
    iter->second.accumulated_functor.apply(context);
    if (!active_sync_tasks.empty()) eval_syncs(globalvid, context, thread::thread_id() % ncpus);
    // scatter
    unsigned char oldkey = rmi.dc().set_sequentialization_key(globalvid % 255 + 1);
    if (context.own_modified) {
      graph.push_owned_vertex_to_replicas(globalvid, true, true);
    }
    foreach(procid_t p, graph.localvid_to_replicas(localvid)) {
      rmi.remote_call(p, 
                      &distributed_chromatic_engine<Graph,UpdateFunctor>::factorized_scatter,
                      globalvid,
                      iter->second.accumulated_functor);
    }
    rmi.dc().set_sequentialization_key(oldkey);
    pending_factorized_updates_lock.lock();
    pending_factorized_updates--;
//    std::cout << rmi.procid() << " scatter issued on " 
//      << globalvid << " " << pending_factorized_updates << " remaining" << std::endl;
    if (pending_factorized_updates == 0) {
      pending_factorized_updates_cond.broadcast();
    }
    pending_factorized_updates_lock.unlock();    
  }
  
  void factorized_scatter(vertex_id_type globalvid,
                          update_functor_type& ufun) {
    dgraph_context<engine_type> context(this, &graph, &shared_data);
    context.init(globalvid, globalvid % ncpus);
    if(ufun.gather_edges() == update_functor_type::IN_EDGES ||
       ufun.gather_edges() == update_functor_type::ALL_EDGES) {
      foreach(const edge_type edge, graph.in_edges_local_only(globalvid)) {
        if (graph.is_owned(edge.source())) ufun.scatter(context, edge);
      }
    }

    if(ufun.gather_edges() == update_functor_type::OUT_EDGES ||
       ufun.gather_edges() == update_functor_type::ALL_EDGES) {
      foreach(const edge_type edge, graph.out_edges_local_only(globalvid)) {
        if (graph.is_owned(edge.target())) ufun.scatter(context, edge);
      }
    }
  }
  
  void evaluate_factorized_update_functor(vertex_id_type localvid, 
                                          update_functor_type& ufun) {
    // Gather phase -----------------------------------------------------------
    vertex_id_type globalvid = graph.localvid_to_globalvid(localvid);
//    std::cout << rmi.procid() << ": Gather on vid " << globalvid << std::endl;
    dgraph_context<engine_type> context(this, &graph, &shared_data);
    context.init(globalvid, localvid % ncpus);
    ufun.init_gather(context);
    if(ufun.gather_edges() == update_functor_type::IN_EDGES ||
       ufun.gather_edges() == update_functor_type::ALL_EDGES) {
      foreach(const edge_type edge, graph.in_edges_local_only(globalvid)) {
        if (graph.is_owned(edge.source())) ufun.gather(context, edge);
      }
    }

    if(ufun.gather_edges() == update_functor_type::OUT_EDGES ||
       ufun.gather_edges() == update_functor_type::ALL_EDGES) {
      foreach(const edge_type edge, graph.out_edges_local_only(globalvid)) {
        if (graph.is_owned(edge.target())) ufun.gather(context, edge);
      }
    }
    
    // send to collected ufun to owner
    rmi.remote_call(graph.localvid_to_owner(localvid),
                    &distributed_chromatic_engine<Graph,UpdateFunctor>::factorized_gather_receive,
                    globalvid,
                    ufun);
  } // end of evaluate_update_functor

  
  
  void factorized_start_thread(size_t threadid) {
    // factorized start start begins with all tasks in the cache set
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
        if (threadid == 0) {
          reset_gather_receive_structures();
          synchronize_vfun_cache(color_block[c]);
          rmi.full_barrier();
          pending_factorized_updates = plan_factorized_schedule(color_block[c]);
          rmi.full_barrier();
        }
        thread_color_barrier.wait();
        
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
            // use the regular update if it is not a seperator set vertex
            if (localvid < graph.local_vertices() && seperator_set.get(localvid) == false) {
              // otherwise. run the vertex
              // create the scope
              context.init(globalvid, threadid);
              // run the update function
              functor_to_run(context);
              // check if there are tasks to run
              if (hassynctasks) eval_syncs(globalvid, context, threadid);
              if (no_graph_synchronization == false) context.commit_async_untracked();
              update_counts[threadid]++;
            }
            else {
              context.init(globalvid, threadid);
              evaluate_factorized_update_functor(localvid, 
                                                functor_to_run);
            }
          }
          else {
            // ok this vertex is not scheduled. But if there are syncs
            // to run I will still need to get the scope
            if (hassynctasks && localvid < graph.local_vertices()) {
              context.init(globalvid, threadid);
              eval_syncs(globalvid, context, threadid);
              if (no_graph_synchronization == false) context.commit_async_untracked();
            }
          }
        }
        pending_factorized_updates_lock.lock();
        while (pending_factorized_updates != 0) {
          pending_factorized_updates_cond.wait(pending_factorized_updates_lock);
        }
        pending_factorized_updates_lock.unlock();    

        //if (threadid == 0 && rmi.procid() == 0) getchar();
        // wait for all threads to synchronize on this color.
        thread_color_barrier.wait();
        curidx.value = 0;
        // full barrier on the color
        // this will complete synchronization of all add tasks as well
        if (threadid == 0) {
          ti.start();
          if (no_graph_synchronization == false) graph.wait_for_all_async_syncs();
          // TODO! If synchronize() calls were made then this barrier is necessary
          // but the time needed to figure out if a synchronize call is required 
          // could be as long as the barrier itself
          rmi.dc().barrier();
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

  float elapsed_time() const { return lowres_time_seconds() - start_time; }

  
  /** Execute the engine */
  void start() {
    start_time = lowres_time_seconds();
    
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
    
    // if factorized create funset extension.
    if (use_factorized) {
      vfunset.resize(graph.get_local_store().num_vertices());
      vfun_cacheset = vfunset;
      vfunset.clear_unsync();
      infer_vertex_seperator_set();
      num_cached_tasks.value = num_pending_tasks.value;
      num_pending_tasks.value = 0;
    }
    else {
      num_cached_tasks.value = 0;
    }
    // two full barrers to complete flush replies
    rmi.dc().full_barrier();
    rmi.dc().full_barrier();

    // reset indices
    curidx.value = 0;
    ti.start();
    // spawn threads
    thread_group thrgrp; 
    if (use_factorized) {
      for (size_t i = 0;i < ncpus; ++i) {
        size_t aff = use_cpu_affinity ? i : -1;
        thrgrp.launch(boost::bind(
                              &distributed_chromatic_engine<Graph, UpdateFunctor>::factorized_start_thread,
                              this, i), aff);
      }
    }
    else {
      for (size_t i = 0;i < ncpus; ++i) {
        size_t aff = use_cpu_affinity ? i : -1;
        thrgrp.launch(boost::bind(
                              &distributed_chromatic_engine<Graph, UpdateFunctor>::start_thread,
                              this, i), aff);
      }
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
      total_update_count = 0;
      for(size_t i = 0; i < procupdatecounts.size(); ++i) {
        total_update_count +=  procupdatecounts[i];
      }
      total_barrier_time = 0;
      for(size_t i = 0; i < barrier_times.size(); ++i) {
        total_barrier_time += barrier_times[i];
      }

      total_bytes_sent = ret["total_bytes_sent"];
    }
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
    opts.engine_args.get_option("use_factorized", use_factorized);
    opts.engine_args.get_option("no_graph_synchronization", no_graph_synchronization);
    opts.engine_args.get_option("factor_threshold", factor_threshold);
    opts.engine_args.get_option("no_colors", no_colors);
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
    out << "use_factorized = [integer, default = 0]\n";
    out << "factor_threshold = [integer, default = 1000]. \n";
    out << "no_graph_synchronization = [boolean, default = 0]. \n";
    out << "no_colors = [boolean, default = 0]. \n";
    out << "   Vertices with higher degree than factor_threshold will not use\n";
    out << "   distributed factorized updates\n";
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

