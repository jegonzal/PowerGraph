/* \file inengine.hpp
   \brief The file containing the iengine description
   
  This file contains the description of the engine interface.  All
  graphlab engines (single_threaded, multi_threaded, distributed, ...)
  should satisfy the functionality described below.
 */

#ifndef IENGINE_HPP
#define IENGINE_HPP

#include <graphlab/graph/graph.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/shared_data/ishared_data.hpp>
#include <graphlab/shared_data/ishared_data_manager.hpp>
namespace graphlab {
  
  /**
   * \brief the reasons for execution completion.
   *
   * Because there are several reasons why the graphlab engine might
   * terminate the exec_status value is returned from the start
   * function after completing execution. 
   *
   */
  enum exec_status {
    
    /** Execution completed successfully due to task depletion */
    EXEC_TASK_DEPLETION,

    /** 
     * Execution completed successfully due to termination
     * function.
     */
    EXEC_TERM_FUNCTION,

    //! The execution completed after timing out
    EXEC_TIMEOUT,

    /** The execution completed because the maximum number of tasks
        was exceeded */
    EXEC_TASK_BUDGET_EXCEEDED,

    /** the engine was stopped by calling force abort */
    EXEC_FORCED_ABORT
  };
  

  
  /**
     \brief The abstract interface of a GraphLab engine.
     
     The graphlab engine interface describes the core functionality
     provided by all graphlab engines.  The engine is templatized over
     the type of graph.
  */
  template<typename Graph>
  class iengine {
  public:

    //! The type of graph that the engine operates on
    typedef Graph graph_type;

    //! The type of update task
    typedef update_task<Graph> update_task_type;

    //! The type of update function
    typedef typename update_task_type::update_function_type 
    update_function_type;

    //! The type of scheduler
    typedef ischeduler<Graph> ischeduler_type;

    //! The type of monitor
    typedef imonitor<Graph> imonitor_type;

    //! The type of scope 
    typedef iscope<Graph> iscope_type;

    //! The type of ishared_data
    typedef ishared_data<Graph> ishared_data_type;

    //! The type of shared data manager
    typedef ishared_data_manager<Graph> ishared_data_manager_type;
    
    /**
     * The termination function is a function that reads the shared
     * data and returns true if the engine should terminate execution.
     * The termination function is called at fixed millisecond
     * intervals and therefore the engine may continue to execute even
     * after a termination function evaluates to true.  Because
     * termination functions are executed frequently and cannot
     * directly contribut to the computation, they should return
     * quickly.
     *
     */
    typedef bool (*termination_function_type) (const ishared_data_type* manager);
    

    //! Virtual destructor required for inheritance 
    virtual ~iengine() {};

    //! get the number of cpus
    virtual size_t get_ncpus() const = 0;


    /**
     * \brief Set the shared data manager
     *
     * If a shared data is to be available to update functions called
     * by this engine then it must be set here.  The shared data can
     * be set to NULL in which case a null pointer will be passed to
     * the update functions.
     */
    virtual void set_shared_data_manager(ishared_data_manager_type* manager) = 0;


    /**
     * \brief Set the default scope range.
     *
     * The default scope range determines the locking extent of an
     * update function. 
     *
     * \param default_scope_range can take on any of the values
     * described in scope_range::scope_range_enum.
     *
     * \todo make sure this documentation properly references the
     * enum.
     * 
     */
    virtual void set_default_scope(scope_range::scope_range_enum default_scope_range) = 0;
    
    /**
     * \brief Start the engine execution.
     *
     * This <b>blocking</b> function starts the engine and does not
     * return until either one of the termination conditions evaluate
     * true or the scheduler has no tasks remaining.
     *
     */
    virtual void start() = 0;


    /**
     * \brief Force engine to terminate immediately.
     *
     * This function is used to stop the engine execution by forcing
     * immediate termination.  Any existing update tasks will finish
     * but no new update tasks will be started and the call to start()
     * will return.
     */
    virtual void stop() = 0;

    
    /**
     * \brief Describe the reason for termination.
     *
     * Return the reason for the last termination.
     *
     */
    virtual exec_status last_exec_status() const = 0;


    
    /**
     * \brief Get the number of updates executed by the engine.
     *
     * This function returns the numbe of updates executed by the last
     * run of this engine.
     * 
     * \return the total number of updates
     *
     */
    virtual size_t last_update_count() const = 0;

        
    /**
     * \brief Register a monitor with an engine. 
     *
     * A monitor tracks the execution of an engine can be useful when
     * debugging. 
     */
    virtual void register_monitor(imonitor_type* listener) = 0;
    
  /**
    * Adds an update task with a particular priority.
    * This function is forwarded to the scheduler.
    */
  virtual void add_task(update_task_type task, double priority) = 0;

  /**
    * Creates a collection of tasks on all the vertices in
    * 'vertices', and all with the same update function and priority
    * This function is forwarded to the scheduler.
    */
  virtual void add_tasks(const std::vector<vertex_id_t>& vertices,
                          update_function_type func, double priority) = 0;

  /**
    * Creates a collection of tasks on all the vertices in the graph,
    * with the same update function and priority
    * This function is forwarded to the scheduler.
    */
  virtual void add_task_to_all(update_function_type func,
                                double priority) = 0;
    /**
     * \brief associate a termination function with this engine.
     *
     * An engine can typically have many termination functions
     * associated with it. A termination function is a function which
     * takes a constant reference to the shared data and returns a
     * boolean which is true if the engine should terminate execution.
     *
     */
    virtual void add_terminator(termination_function_type term) = 0;

    //!  remove all associated termination functions
    virtual void clear_terminators() = 0;
    

    virtual void enable_sched_yield(bool value) { }

    virtual void enable_cpu_affinities(bool value) { }

    
    
    /**
     *  Timeout. Default - no timeout. The timeout is the total
     *  ammount of time in seconds that the engine may run before
     *  exeuction is automatically terminated.
     *
     * \todo Should we continue to support this function. 
     */
    virtual void set_timeout(size_t timeout_secs) = 0;
    
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
     *
     * \todo Should we continue to support this function?
     */
    virtual void set_task_budget(size_t max_tasks) = 0;


    
    virtual scheduler_options& sched_options() = 0;

    virtual const scheduler_options& sched_options() const = 0;
  };

}

#endif

