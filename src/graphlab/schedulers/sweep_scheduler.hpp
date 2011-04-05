/**
 * This class defines a very simple scheduler that loops vertices that
 * are "dirty". Each cpu loops vertices of which id%num_cpus==cpuid. Also called
 * "partitioned" scheduler.
 **/

#ifndef GRAPHLAB_SWEEP_SCHEDULER_HPP
#define GRAPHLAB_SWEEP_SCHEDULER_HPP

#include <queue>
#include <cmath>
#include <cassert>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/synchronized_queue.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
#include <graphlab/schedulers/support/vertex_task_set.hpp>
#include <graphlab/util/shared_termination.hpp>

#include <graphlab/parallel/atomic.hpp>


#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** \ingroup group_schedulers
    */
  template<typename Graph>
  class sweep_scheduler : 
    public ischeduler<Graph> {
  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type callback_type;
    typedef typename base::monitor_type monitor_type;
    typedef shared_termination terminator_type;
    
  private:
    using base::monitor;
   
    static const size_t LS_MAX_UPDATEFUNCTIONS = 3;
    
    size_t INTERNAL_TO_VID(size_t vid) { 
      return v_to_int[vid]; 
    }
    
    size_t VID_TO_INTERNAL(size_t intid) { 
      return int_to_v[intid]; 
    }

    size_t GETIDX(size_t funcid, size_t cpuid, size_t vid) {
      return funcid * vrcpu + cpuid * vreserve + vid;
    }

    size_t vreserve;
    size_t vrcpu;

    vertex_id_t* v_to_int;
    vertex_id_t* int_to_v;

  public:
    sweep_scheduler(iengine_type* engine,
                    Graph& g, 
                    size_t ncpus) : 
      terminator(ncpus) {
      this->g = &g;
      numvertices = g.local_vertices();
      num_cpus = ncpus;
      callbacks.resize(num_cpus);
      for(unsigned int i=0; i<num_cpus; i++) {
        callbacks[i] = new direct_callback<Graph>(this, engine);
      }
      // pad to guarantee not go over boundary
      vreserve = LS_MAX_UPDATEFUNCTIONS * (numvertices + 8 - numvertices%8); 
      vrcpu = vreserve * ncpus;
      dirty_vertices = (char*) calloc(vrcpu * 
                                      LS_MAX_UPDATEFUNCTIONS, sizeof(char));
      last_index = std::vector<size_t>(ncpus,0);
      for(unsigned int i=0; i<num_cpus; i++) {
        last_index[i] = i;
      }
      last_col = std::vector<uint8_t>(ncpus, 0);
      updatefuncs = (update_function_type *) malloc(LS_MAX_UPDATEFUNCTIONS * 
                                                    sizeof(update_function_type));
      num_of_updatefunctions = 0;
      v_to_int = NULL;
      int_to_v = NULL;
      permute_vertices = false;
    }
    
    
    
    ~sweep_scheduler() { 
      delete dirty_vertices;
      delete updatefuncs;
      free(v_to_int);
      free(int_to_v);
    }
    
    void start() { }
    
    /** Get next dirty vertex in queue. Each cpu checks vertices with
        modulo num_cpus = cpuid */
    sched_status::status_enum get_next_task(size_t cpuid,
                                 update_task_type &ret_task) {
      size_t lastidx = last_index[cpuid];
      
      if (lastidx >= numvertices) lastidx = cpuid;
            
      unsigned int nvp = numvertices + num_cpus;
      for(unsigned int i=0; i<nvp; i+=num_cpus) {
       
        vertex_id_t vid = lastidx;
        assert(vid%num_cpus == cpuid);
        assert(vid >= 0 && vid < numvertices);
        
        for(int col=last_col[cpuid]; col<num_of_updatefunctions; col++) {
          if (dirty_vertices[GETIDX(col,cpuid,vid)]) {
            dirty_vertices[GETIDX(col,cpuid,vid)] = 0;
            ret_task = update_task_type(INTERNAL_TO_VID(vid), updatefuncs[col]);
            if (col < num_of_updatefunctions-1) {
              // Maybe more functions left scheduled for this vertex
              last_index[cpuid] = lastidx;
              last_col[cpuid] = col+1;
            } else {
              // Was last update function, move to next vertex
              last_index[cpuid] = (lastidx+num_cpus);
              last_col[cpuid] = 0;
            }
            if (monitor != NULL) 
              monitor->scheduler_task_scheduled(ret_task, 0.0);
            
            return sched_status::NEWTASK;
          }
        }
        lastidx += num_cpus;
        if (lastidx >= numvertices) lastidx = cpuid;
        last_col[cpuid] = 0;
      }
      // printf("Cpu %d Task new %d done%d\n", cpuid, task_new.value,
      // task_done.value);
      return sched_status::EMPTY;
    }


    uint8_t get_update_func_id(update_function_type upf) {
      for(int i=0; i<num_of_updatefunctions; i++ ){
        if (updatefuncs[i] == upf) return i;
      }
      // Ok not found, have to make it 
      updflock.lock();
      // Check once more
      for(int i=0; i<num_of_updatefunctions; i++ ){
        if (updatefuncs[i] == upf) return i;
      }
      assert(num_of_updatefunctions < LS_MAX_UPDATEFUNCTIONS);
      updatefuncs[num_of_updatefunctions] = upf;
      int newid = num_of_updatefunctions;
      num_of_updatefunctions++;
      updflock.unlock();
      return newid;
    }
    
    
    /** Simple sets a flag to dirty **/
    void add_task(update_task_type task, double priority, 
                  int generated_by_cpuid) {
      if (v_to_int == NULL) init();  // Check if init was forgotten!
      vertex_id_t task_vid = VID_TO_INTERNAL(task.vertex());
      size_t targetcpuid = task_vid%num_cpus;
      size_t idx = GETIDX(get_update_func_id(task.function()), 
                          targetcpuid , task_vid);
      if (dirty_vertices[idx] == 0) {
        dirty_vertices[idx] = 1;
        terminator.new_job(targetcpuid);
        if (monitor != NULL) 
          monitor->scheduler_task_added(task, priority);
      }
      if (monitor != NULL) 
        monitor->scheduler_task_pruned(task);
    }
    
    
      
    bool is_task_scheduled(update_task_type task) {
      vertex_id_t task_vid = VID_TO_INTERNAL(task.vertex_id());
      size_t targetcpuid = task_vid%num_cpus;
      size_t idx = GETIDX(get_update_func_id(task.function()), 
                          targetcpuid , task_vid);
      return (dirty_vertices[idx] != 0);
    }

          
    
    callback_type& get_callback(size_t cpuid) {
      return *(callbacks[cpuid]);
    }
    
    void add_tasks(const std::vector<vertex_id_t> &vertices,
                   update_function_type func,
                   double priority)  {
      foreach(vertex_id_t vertex, vertices) {
        add_task(update_task_type(vertex, func), priority, -1);
      }
    }
    
    void add_task_to_all(update_function_type func, double priority) {
      for (vertex_id_t vertex = 0; vertex < numvertices; ++vertex){
        add_task(update_task_type(vertex, func), priority, -1);
      }
    }
    
    
    void add_task(update_task_type task, double priority) {
      add_task(task, priority, -1);
    }
  
    void init() {
      /* Mappings */
      v_to_int = (vertex_id_t*) calloc(numvertices, sizeof(vertex_id_t));
      int_to_v = (vertex_id_t*) calloc(numvertices, sizeof(vertex_id_t));
      if (permute_vertices) {
        random::shuffle(&(v_to_int[0]), &(v_to_int[numvertices]));
        for(vertex_id_t i=0; i<numvertices; i++) {
          int_to_v[v_to_int[i]] = i;
        }
      }
      else {
        for(vertex_id_t i=0; i<numvertices; i++) {
          v_to_int[i] = i;
          int_to_v[i] = i;
        }
      }
    }
    
    void completed_task(size_t cpuid, const update_task_type &task) { }


    terminator_type& get_terminator() {
      return terminator;
    };


    void set_options(const scheduler_options &opts) {
      std::string ordertype;
      if (opts.get_string_option("ordering", ordertype)) {
        if (ordertype == "permute") {
          permute_vertices = true;
        }
        else if (ordertype == "linear") {
          permute_vertices = false;
        }
        else {
          logstream(LOG_WARNING) << "Invalid Ordering Method" << ordertype
                                                               << std::endl;
        }
      }
    }
    static void print_options_help(std::ostream &out) {
      out << "ordering = [string: linear/permute, default=linear]\n";
    };

  private:
    
    char* dirty_vertices;
    std::vector<size_t> last_index;
    std::vector<uint8_t> last_col;
    uint8_t num_of_updatefunctions;
    update_function_type* updatefuncs;
    spinlock updflock;
    shared_termination terminator;
    size_t numvertices;
    size_t num_cpus;
    std::vector< direct_callback<Graph>* > callbacks;
    Graph* g;
    bool permute_vertices;
  }; 
  
  
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
