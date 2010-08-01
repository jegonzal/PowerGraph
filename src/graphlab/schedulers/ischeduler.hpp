#ifndef GRAPHLAB_ISCHEDULER_HPP
#define GRAPHLAB_ISCHEDULER_HPP

#include <vector>

#include <graphlab/engine/iengine.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/schedulers/icallback.hpp>


namespace graphlab {
        
  
  /**
   * This is an enumeration for the possible return values for
   * get_next_tasks
   */
  struct sched_status {
    enum status_enum {
      NEWTASK,      /// The get_next_tasks function returned a new task
                    /// to be executed
      WAITING,      /// The get_next_tasks function did not return a new
                    /// task, but the program should not terminate
      COMPLETE      /// The get_next_tasks function did not return a new
                    /// task, and the program should terminate
    };
  };
  
  struct scheduler_options {
    enum options_enum {
      UPDATE_FUNCTION,
      SPLASH_SIZE,
      MAX_ITERATIONS,
      START_VERTEX,
      SCHEDULING_FUNCTION,
      BARRIER,
      DISTRIBUTED_CONTROL
    };
  };


  /**
   * This describes the interface/concept for the scheduler. The
   * engine will be passed the scheduler type as a template argument,
   * so the scheduler must inherit and satisfy this interface
   * EXACTLY. Note that all functions (with the exception of the
   * constructor and destructor) must be thread-safe.
   */
  template<typename Graph>
  class ischeduler {
  public:

    typedef Graph graph_type;
    typedef update_task<Graph> update_task_type;
    typedef typename update_task_type::update_function_type 
    update_function_type;

    typedef iengine<Graph> iengine_type;
    typedef icallback<Graph> callback_type;
    typedef imonitor<Graph> monitor_type;
    
    /** 
     * Constructor: The scheduler must be provided with the graph, and the 
     * number of
     * cpus.  All initialization of the scheduler internal state must
     * be done here.
     *
     *    \note This constructor here does not actually do
     *     anything. It just exists to force the derived class
     *     constructors to look like this
     */
    //    ischeduler(iengine_type* engine, Graph& g, size_t ncpus) : monitor(NULL) { }
    ischeduler() : monitor(NULL) {}
    
    /// destructor
    virtual ~ischeduler() {};
        
    /** Called by engine before executing the schedule */
    virtual void start() {};

    /** Called when the engine stops */
    virtual void stop() {};

    
    /// Adds an update task with a particular priority
    virtual void add_task(update_task_type task, double priority) = 0;
    
    /** 
     * Creates a collection of tasks on all the vertices in
     * 'vertices', and all with the same update function and priority
     */
    virtual void add_tasks(const std::vector<vertex_id_t>& vertices, 
                           update_function_type func, double priority) = 0;
    
    /** 
     * Creates a collection of tasks on all the vertices in the graph,
     * with the same update function and priority
     */
    virtual void add_task_to_all(update_function_type func, 
                                 double priority) = 0;
    
    /**
     * This function returns a reference to the scheduling callback to
     * be used for a particular cpu. This callback will be passed to
     * update functions, and is the main interface which allow the
     * update functions to create new tasks.
     */
    virtual callback_type& get_callback(size_t cpuid) = 0;

    /**
     * This function is called by the engine to ask for new work to
     * do.  The update task to be executed is returned in ret_task.
     *
     *  \retval NEWTASK There is an update task in ret_task to be
     *   executed
     * 
     *  \retval WAITING ret_task is empty. But the engine should wait
     *  as execution is still not complete
     *
     *  \retval COMPLETE ret_task is empty and the engine should
     *  proceed to terminate
     */
    virtual sched_status::status_enum get_next_task(size_t cpuid, 
                                       update_task_type &ret_task) = 0;

    /**
     * This is called after a task has been executed
     */
    virtual void completed_task(size_t cpuid, 
                                const update_task_type &task) = 0;


    /** Installs a listener (done by the engine) */
    virtual void register_monitor(monitor_type* monitor_) { 
      monitor = monitor_;
    }        

    virtual void set_option(scheduler_options::options_enum opt, void* value) { };


  protected:
    monitor_type* monitor;

  };

}
#endif

