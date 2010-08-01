/**
  * \file Clustered priority scheduler.
  **/
#ifndef GRAPHLAB_CLUSTERED_PRIORITY_SCHEDULER_HPP
#define GRAPHLAB_CLUSTERED_PRIORITY_SCHEDULER_HPP

#include <cmath>
#include <cassert>
#include <queue>
#include <vector>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
#include <graphlab/schedulers/support/binary_vertex_task_set.hpp>
#include <graphlab/util/shared_termination.hpp>
#include <graphlab/util/dense_bitset.hpp>

// Not used
// #include <graphlab/extern/bitmagic/bm.h>

#include <graphlab/macros_def.hpp>

namespace graphlab {

  template<typename Graph>
  class clustered_priority_scheduler : 
    public ischeduler<Graph> {
  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type callback_type;
    typedef typename base::monitor_type monitor_type;

  private:
    using base::monitor;
 
  private:
    
    
    Graph& g;

    /** Remember the number of vertices in the graph */
    size_t num_vertices;  
    
    /** The queue over vertices */
    mutable_queue<size_t, double> task_queue;

    /** The lock on the priority queue */
    spinlock queuelock; 
    
    std::vector<std::vector<uint32_t> > id2vertex;
    std::vector<uint32_t> vertex2id;

    /** The vertex task set which maintains the tasks at each
        vertex */
    binary_vertex_task_set<Graph> task_set;
  
    /** The callbacks pre-created for each cpuid */
    std::vector<direct_callback<Graph> > callbacks;
    
    /** Used to assess termination */
    shared_termination terminator;
    
    size_t ncpus;
    bool aborted;
    
    struct cpustate {
      int clusterid; // the cluster ID this cpu is currently running
      size_t cluster_offset;  // the next "unread" offset in the
                              // cluster;
      uint64_t update_function_list;  // current list of update
                                      // functions to be executed on
                                      // the vertex v
      vertex_id_t vertex;
      std::vector<update_task_type> buffered_task_creations;
      std::vector<double> buffered_priority_updates;
      dense_bitset hasbuffer;
    };

    std::vector<cpustate> cpu_state;

  public:
    
    clustered_priority_scheduler(iengine_type* engine,
                                 Graph& _g, 
                                 size_t ncpus_) :
      g(_g),
      num_vertices(_g.num_vertices()),
      task_set(_g.num_vertices()),
      callbacks(ncpus_, direct_callback<Graph>(this, engine) ),
      terminator(ncpus_),
      ncpus(ncpus_) { 
      cpu_state.resize(ncpus);
      aborted = false;
      for (size_t i = 0; i < ncpus; ++i) {
        cpu_state[i].clusterid = -1;
        cpu_state[i].cluster_offset = 0;
        cpu_state[i].update_function_list = 0;
        cpu_state[i].vertex = 0;
        cpu_state[i].buffered_task_creations.reserve(4096);
      }
    }
    

    ~clustered_priority_scheduler() { }
    
    callback_type& get_callback(size_t cpuid) {
      return callbacks[cpuid];
    }
    
    void cluster(size_t verticespercluster,
                 partition_method::partition_method_enum method) {
      size_t numclusters = g.num_vertices() / verticespercluster;
      if (numclusters == 0) {
        numclusters = 1;
      }
      // cluster the vertices
      timer ti;
      ti.start();
      g.partition(method,numclusters, vertex2id);
      
      logger(LOG_INFO, "Partition took %f seconds", ti.current_time());

      id2vertex.resize(numclusters);
      for (size_t i = 0 ;i < g.num_vertices(); ++i) {
        id2vertex[vertex2id[i]].push_back(i);
      }
      for (size_t i = 0;i < numclusters; ++i) {
        task_queue.push(i,0.0);
      }
      
      for (size_t i = 0; i < ncpus; ++i) {
        cpu_state[i].buffered_priority_updates.resize(numclusters);
        cpu_state[i].hasbuffer.resize(numclusters);
        cpu_state[i].hasbuffer.clear();
      }
    }
    
    /** Get the next element in the queue */
    sched_status::status_enum get_next_task_from_list(size_t cpuid, 
                                         update_task_type &ret_task) {
      // Loop until we either can't build a splash or we find an
      // element
      cpustate &curstate = cpu_state[cpuid];
      while(true) {
        if (aborted) return sched_status::WAITING;
        // check for work to do in the updatefunction list
        if (curstate.update_function_list != 0) {
          // look for a bit to pop
          size_t bitpos = __builtin_ctzl(curstate.update_function_list);
          const size_t mask(size_t(1) << size_t(bitpos)); 
          curstate.update_function_list = 
            curstate.update_function_list  & (~mask);
          // return the function
          ret_task = 
            update_task_type(curstate.vertex, task_set.updatefuncs[bitpos]);
          if (monitor != NULL) 
            monitor->scheduler_task_scheduled(ret_task, 0.0);
          return sched_status::NEWTASK;
        } else if (curstate.clusterid >= 0 && 
                   curstate.cluster_offset < 
                   id2vertex[curstate.clusterid].size()) {
          // if update_function_list is empty
          // try to get something from the cluster
         
          // recall clusteroffset is the"next vertex"
          for(; curstate.cluster_offset < 
                id2vertex[curstate.clusterid].size();
              curstate.cluster_offset++ ) {
            // try to pop the list of update functions
            curstate.vertex = 
              id2vertex[curstate.clusterid][curstate.cluster_offset];
            curstate.update_function_list = 
              task_set.pop_all_tasks(curstate.vertex); 
            // if there are update functions, increment the vertex
            // list and break
            if (curstate.update_function_list != 0) {
              ++curstate.cluster_offset;
              break;
            }
          }

        } else {         // no update, nothing in the cluster.
          // before we get the next list, we need to check this one in
          curstate.clusterid = -1;
          completed_task(cpuid, update_task_type(0, NULL));
          get_next_list(cpuid);
          // if I can't get next list
          if (curstate.clusterid < 0) {
            return sched_status::WAITING;
          }
        }
      }
    } // end of get next task
    
    
    sched_status::status_enum get_next_task(size_t cpuid, update_task_type &ret_task) {
      // While the scheduler is active
      while(true) {
        // Try and get next task for splash
        sched_status::status_enum ret = get_next_task_from_list(cpuid, ret_task);
        // If we are not waiting then just return
        if (ret != sched_status::WAITING) {
          return ret;        
        } else {
          // Otherwise enter the shared terminator code
          terminator.begin_sleep_critical_section(cpuid);
          // Try once more to get the next task
          ret = get_next_task_from_list(cpuid, ret_task);
          // If we are waiting then 
          if (ret != sched_status::WAITING) {
            // If we are either complete or succeeded then cancel the
            // critical section and return
            terminator.cancel_sleep_critical_section(cpuid);
            return ret;
          } else {
            // Otherwise end sleep waiting for either completion
            // (true) or some new work to become available.
            if (terminator.end_sleep_critical_section(cpuid)) {
              return sched_status::COMPLETE;
            }
          } 
        }
      } // End of while loop
      // We reach this point only if we are no longer active
      return sched_status::COMPLETE;
    }
    
    void start() {
      // need to dump CPU0
      completed_task(0, update_task_type(0, NULL));
      terminator.reset();
    }
    
    void get_next_list(size_t cpuid) {
      
      cpustate &curstate = cpu_state[cpuid];
      queuelock.lock();
      // get the top, demote it to the bottom.
      std::pair<size_t, double> t = task_queue.top();
      if (t.second > 0) {
        //demote it to the bottom
        //full the states
        task_queue.update(t.first, 0);
        curstate.clusterid = t.first;
        curstate.cluster_offset = 0;
        curstate.update_function_list = 0;
        curstate.vertex = 0;
        // std::cout << "cpuid " << cpuid << " got cluster " <<
        // t.first << " with priority " << t.second << "\n";
      }
      else {
        curstate.clusterid = -1;
      }
      queuelock.unlock();
    }
    
    void add_task(update_task_type task, double priority) {
      size_t tid = thread::thread_id();
      // buffer the change in priority
      size_t clusterid = vertex2id[task.vertex()];
      cpu_state[tid].buffered_priority_updates[clusterid ] = 
        std::max(cpu_state[tid].buffered_priority_updates[clusterid], 
                 priority);
      // set the bit so we know which cluster has changed
      cpu_state[tid].hasbuffer.set_bit(clusterid);
      cpu_state[tid].buffered_task_creations.push_back(task);
    } // end of add_task
    

    void add_tasks(const std::vector<vertex_id_t> &vertices,
                   update_function_type func,
                   double priority) {
      foreach(vertex_id_t vertex, vertices) {
        add_task(update_task_type(vertex, func), priority);
      }
    } // end of add tasks 
    
    
    void add_task_to_all(update_function_type func, double priority) {
      for (vertex_id_t vertex = 0; vertex < num_vertices; ++vertex){
        add_task(update_task_type(vertex, func), priority);
      }
    } // end of add tasks to all
    
    void update_state(size_t cpuid,
                      const std::vector<vertex_id_t> &updated_vertices,
                      const std::vector<edge_id_t>& updatededges) { 
    }

    void scoped_modifications(size_t cpuid, vertex_id_t rootvertex,
                              const std::vector<edge_id_t>& updatededges){}

    void completed_task(size_t cpuid, const update_task_type &task) {
      // if I am at the end of the current list, dump the buffer
      cpustate &curstate = cpu_state[cpuid];
      if (curstate.clusterid == -1 || 
          (curstate.cluster_offset == id2vertex[curstate.clusterid].size())) {
        bool insertionsmade = 
          curstate.buffered_task_creations.size() > 0;
        // insert the tasks
        for (size_t i = 0; i < curstate.buffered_task_creations.size(); ++i) {
          task_set.add(curstate.buffered_task_creations[i]);
        }
        curstate.buffered_task_creations.clear();
  
        uint32_t b = 0;
        // update the the priority queue
        if (curstate.hasbuffer.first_bit(b) == false) return;
        queuelock.lock();
        do {
          task_queue.insert_max(b, curstate.buffered_priority_updates[b]);
          curstate.buffered_priority_updates[b] = 0;
          curstate.hasbuffer.clear_bit(b);
        } while(curstate.hasbuffer.next_bit(b));
        queuelock.unlock();
        if (insertionsmade) {
          terminator.new_job();
        }
      }
    }
    
    void abort() { aborted=true; }
    
    void restart() { 
      aborted=false;
      for (size_t i = 0; i < ncpus; ++i) {
        cpu_state[i].clusterid = -1;
        cpu_state[i].cluster_offset = 0;
        cpu_state[i].update_function_list = 0;
        cpu_state[i].vertex = 0;
        cpu_state[i].buffered_task_creations.clear();

        for (size_t j = 0; 
             j < cpu_state[i].buffered_priority_updates.size(); ++j) {
          cpu_state[i].buffered_priority_updates[j] = 0;;
        }
        cpu_state[i].hasbuffer.clear();
      }
      completed_task(0, update_task_type(0,NULL)); 
    }
    
  }; // end of priority_queue class


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
