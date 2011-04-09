/* \file iengine.hpp
   \brief The file containing the iengine description
   
   This file contains the description of the engine interface.  All
   graphlab engines (single_threaded, multi_threaded, distributed, ...)
   should satisfy the functionality described below.
*/

#ifndef GRAPHLAB_IENGINE_HPP
#define GRAPHLAB_IENGINE_HPP

#include <graphlab/graph/graph.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/shared_data/ishared_data.hpp>
#include <graphlab/shared_data/glshared.hpp>
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

    EXEC_UNSET,  /** The default termination reason */

    EXEC_TASK_DEPLETION, /**<Execution completed successfully due to
                            task depletion */

    EXEC_TERM_FUNCTION,  /**< Execution completed successfully due to
                            termination function. */

    EXEC_TIMEOUT,       /**< The execution completed after timing
                           out */

    EXEC_TASK_BUDGET_EXCEEDED, /**< The execution completed because
                                  the maximum number of tasks was
                                  exceeded */

    EXEC_FORCED_ABORT,     /**< the engine was stopped by calling force
                             abort */
                             
    EXEC_EXCEPTION        /**< the engine was stopped by an exception */
  };
  

  
  /**
     \brief The abstract interface of a GraphLab engine.
     The graphlab engine interface describes the core functionality
     provided by all graphlab engines.  The engine is templatized over
     the type of graph.
     
     The GraphLab engines are a core element of the GraphLab
     framework.  The engines are responsible for applying a the update
     tasks and sync operations to a graph and shared data using the
     scheduler to determine the update schedule. This class provides a
     generic interface to interact with engines written to execute on
     different platforms.
     
     While users are free to directly instantiate the engine of their
     choice we highly recommend the use of the \ref core data
     structure to manage the creation of engines. Alternatively, users
     can use the \ref engine_factory static functions to create
     engines directly from configuration strings.
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

    typedef void(*sync_function_type)(iscope_type& scope,
                                      any& accumulator);

    typedef void(*merge_function_type)(any& merge_dest,
                                       const any& merge_src);

    
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
     * \deprecated Use set_sync and glshared
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
     * update function. See \ref Scopes for details.
     *
     * \param default_scope_range can take on any of the values
     * described in \ref scope_range
     *
     */
    virtual void set_default_scope(scope_range::scope_range_enum default_scope_range) = 0;
    
    /**
     * \brief Start the engine execution.
     *
     * This \b blocking function starts the engine and does not
     * return until either one of the termination conditions evaluate
     * true or the scheduler has no tasks remaining.
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
     */
    virtual exec_status last_exec_status() const = 0;


    
    /**
     * \brief Get the number of updates executed by the engine.
     *
     * This function returns the numbe of updates executed by the last
     * run of this engine.
     * 
     * \return the total number of updates
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
     * \brief Adds an update task with a particular priority.
     * This function is forwarded to the scheduler.
     */
    virtual void add_task(update_task_type task, double priority) = 0;

    /**
     * \brief Add an update function to a particular vertex.
     */
    virtual void add_vtask(vertex_id_t vid, 
                          update_function_type fun, 
                          double priority = 1.0) {
      add_task(update_task_type(vid, fun),  priority);
    }

    /**
     * \brief Creates a collection of tasks on all the vertices in
     * 'vertices', and all with the same update function and priority
     * This function is forwarded to the scheduler.
     */
    virtual void add_tasks(const std::vector<vertex_id_t>& vertices,
                           update_function_type func, double priority) = 0;

    /**
     * \brief Creates a collection of tasks on all the vertices in the graph,
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
     * A termination function has the following type:
     * \code
     * bool term_fun(const ishared_data_type* shared_data)
     * \endcode
     */
    virtual void add_terminator(termination_function_type term) = 0;

    //!  remove all associated termination functions
    virtual void clear_terminators() = 0;
    

    /**
     * Set whether sched yield should be used when waiting on new
     * jobs
     */
    virtual void set_sched_yield(bool value) = 0;

    /**
     * Set whether cpu affinities should be used.
     */
    virtual void set_cpu_affinities(bool value) = 0;

    
    
    /**
     *  \brief The timeout is the total
     *  ammount of time in seconds that the engine may run before
     *  exeuction is automatically terminated.
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
     */
    virtual void set_task_budget(size_t max_tasks) = 0;


    /** \brief Update the scheduler options.  */
    virtual void set_scheduler_options(const scheduler_options& opts) = 0;



    /**
     * Registers a sync with the engine.
     * The sync will be performed every "interval" updates,
     * and will perform a reduction over all vertices from rangelow
     * to rangehigh inclusive.
     * The merge function may be NULL, in which it will not be used.
     *
     * \param shared The shared variable to synchronize
     * \param sync The reduction function
     * \param apply The final apply function which writes to the shared value
     * \param zero The initial zero value passed to the reduction
     * \param sync_interval Frequency at which the sync is initiated.
     *                      Corresponds approximately to the number of
     *                     update function calls before the sync is reevaluated.
     *                     If 0, the sync will only be evaluated once
     *                     at engine start,  and will never be evaluated again.
     *                     Defaults to 0.
     * \param merge Combined intermediate reduction value. defaults to NULL.
     *              in which case, it will not be used.
     * \param rangelow he lower range of vertex id to start syncing.
     *                 The range is inclusive. i.e. vertex with id 'rangelow'
     *                 and vertex with id 'rangehigh' will be included.
     *                 Defaults to 0.
     * \param rangehigh The upper range of vertex id to stop syncing.
     *                  The range is inclusive. i.e. vertex with id 'rangelow'
     *                  and vertex with id 'rangehigh' will be included.
     *                  Defaults to infinity.
     */
    virtual void set_sync(glshared_base& shared,
                          sync_function_type sync,
                          glshared_base::apply_function_type apply,
                          const any& zero,
                          size_t sync_interval = 0,
                          merge_function_type merge = NULL,
                          size_t rangelow = 0,
                          size_t rangehigh = -1) { }

    /**
     * Performs a sync immediately. This function requires that the shared
     * variable already be registered with the engine.
     */
    virtual void sync_now(glshared_base& shared) { };
    
    // Convenience function.
    static std::string exec_status_as_string(exec_status es) {
      switch(es) {
      case EXEC_UNSET: return "engine not run!";
      case EXEC_FORCED_ABORT: return "forced abort";
      case EXEC_TASK_BUDGET_EXCEEDED: return "budget exceed";
      case EXEC_TERM_FUNCTION: return "termination function";
      case EXEC_TASK_DEPLETION: return "task depletion (natural)";
      case EXEC_TIMEOUT: return "timeout";
      case EXEC_EXCEPTION: return "exception";
      };
      return "unknown";
    }

  };

}

#endif

