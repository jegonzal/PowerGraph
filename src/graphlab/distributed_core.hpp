/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_DISTRIBUTED_CORE_HPP
#define GRAPHLAB_DISTRIBUTED_CORE_HPP

#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/engine_options.hpp>
#include <graphlab/distributed2/distributed2_includes.hpp>
#include <graphlab/util/command_line_options.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>

#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/graph/graph.hpp>



#include <graphlab/metrics/metrics.hpp>
#include <graphlab/metrics/reporters/null_reporter.hpp>
#include <graphlab/metrics/reporters/basic_reporter.hpp>
#include <graphlab/metrics/reporters/file_reporter.hpp>
#include <graphlab/metrics/reporters/html_reporter.hpp>



#include <graphlab/macros_def.hpp>
namespace graphlab {

  // Predecleration 
  template<typename Graph> struct types;
  


  /**
     \brief A GraphLab core is the base (or core) data structure in GraphLab.
     
     This is like \ref graphlab::core but for the distributed setting.

     The core is templatized over the VertexType and EdgeType however
     by using the ref types typedef, one can simply create a core by
     doing the following:
   
     \code
     gl::distributed_core glcore;
     \endcode
   
     The core contains the 
   
     \li Data Graph: which represents the structured data dependencies.
     \li Engine: The computational structure which contains the
     scheduling and execution statistics for the GraphLab program. The
     core provides pass-through calls for many engine functions.
        
     The core also manages the engine and scheduler construction
     parameters.
     
    The distributed core is more limited as compared to the 
    shared memory \ref graphlab::core version. In particular, engine construction
    must be executed manually through build_engine() and the
    engine options / scheduler options cannot be modified after engine construction.
    
    Also, some functions must be called by all machines simultaneously, 
    while others are "parallel" allowing any machine to call the function 
    seperately. This behavior is documented in each function. The user must
    take care to obey this requirement or it may result in unexpected behavior.
  */
  template <typename VertexType, typename EdgeType>
  class distributed_core {
  public:
    typedef graphlab::types<graphlab::distributed_graph<VertexType, EdgeType> > types;

  public:
    /** default constructor. Graph is constructed using the atom index.
     * All machines must construct simultaneously.
    */
    distributed_core(distributed_control &dc, std::string atomindex) :
      dc(dc),
      mgraph(dc, atomindex),
      mengine(NULL),
      coremetrics("distributed_core"), reporter(new null_reporter) { }
  private:
    //! Core is not copyable
    distributed_core(const distributed_core& other);
    //! Core is not copyable
    distributed_core& operator=(const distributed_core& other);

  public:

    /**
     * Destructor. 
     * All machines must call simultaneously.
     */
    ~distributed_core() { 
      if (meopts.get_metrics_type() != "none") {        
        // Write options to metrics
        fill_metrics();
        report_metrics();
      }
      delete mengine;
      delete reporter;
    } 
       
    /** Get a modifiable reference to the graph associated with this core
     * This function is parallel.
     */
    typename types::graph& graph() { return mgraph; }

    /** Get a constant reference to the graph associated with this core
     * This function is parallel.
     */
    const typename types::graph& graph() const { return mgraph; }

    /**
     * \brief Set the type of scheduler.
     * The engine must not be constructed yet.
     * All machines must call simultaneously.
     */
    void set_scheduler_type(const std::string& scheduler_type) {
      ASSERT_EQ(mengine, NULL);
      bool success = meopts.set_scheduler_type(scheduler_type);
      ASSERT_TRUE(success);
    }

    /**
     * \brief Set the scope consistency model used in this engine.
     *
     * The engine must not be constructed yet. 
     * All machines must call simultaneously.
     * The available scopes are:
     * 
     *  \li \b "full" This ensures full data consistency within the scope
     *  \li \b "edge" This ensures data consistency with just the
     *     vertex and edges
     *  \li \b "vertex" This ensures that a vertex cannot be updated
     *     by two processors simultaneously
     *
     * See \ref Scopes for details
     */
    void set_scope_type(const std::string& scope_type) {
      ASSERT_EQ(mengine, NULL);
      bool success = meopts.set_scope_type(scope_type);
      ASSERT_TRUE(success);
    }


    /**
     * \brief Set the engine type.
     *
     * The engine must not be constructed yet. 
     * All machines must call simultaneously.
     *
     *  \li \b "dist_locking" Distributed engine with consistency ensured 
     *                        through locking
     *  \li \b "dist_chromatic" Distributed engien with consistency ensured
     *                          through coloring
     */
    void set_engine_type(const std::string& engine_type) {
      ASSERT_EQ(mengine, NULL);
      bool success = meopts.set_engine_type(engine_type);
      ASSERT_TRUE(success);
    }
    
    /**
     * \brief Sets the output format of any recorded metrics
     *  This function is parallel.
     * 
     *  \li \b "none" No reporting
     *  \li \b "basic" Outputs to screen
     *  \li \b "file" Outputs to a text file graphlab_metrics.txt
     *  \li \b "html" Outputs to a html file graphlab_metrics.html
     */
    void set_metrics_type(const std::string& metrics_type) {
      bool metrics_set_success = meopts.set_metrics_type(metrics_type);
      ASSERT_TRUE(metrics_set_success);
      
      delete reporter;
      if (meopts.get_metrics_type() == "file") {
        reporter = new file_reporter("graphlab_metrics.txt");
      } else if (meopts.get_metrics_type() == "html") {
        reporter = new  html_reporter("graphlab_metrics.html");
      } else if (meopts.get_metrics_type() == "basic") {
        reporter = new basic_reporter;
      } else {
        reporter = new null_reporter;
      }
    }

    
    /**
     * \brief Set the number of cpus that the engine will use.
     *
     * The engine must not be constructed yet. 
     * All machines must call simultaneously.
     *
     */
    void set_ncpus(size_t ncpus) {
      ASSERT_EQ(mengine, NULL);
      meopts.set_ncpus(ncpus);
    }


    /**
     * Get a reference to the active engine.  
     * build_engine() must be called prior to this.
     * This function is parallel.
     */
    typename types::iengine& engine() {
      ASSERT_NE(mengine, NULL);
      return *mengine; 
    }




    /**
     * \brief Constructs the engine using the current defined options
     * Once an engine is constructed, options cannot be modified
     * All machines must call simultaneously.
     */
    bool build_engine() {
      ASSERT_EQ(mengine, NULL);
      // create the engine
      mengine = distributed_engine_factory::new_engine(dc, meopts, mgraph);
      if(mengine == NULL) return false;
      return true;
    }

    /**
     * \brief Set the engine options by passing in an engine options object.
     * The engine must not be constructed yet. 
     * All machines must call simultaneously.
     */
    void set_engine_options(const engine_options& opts) {
      ASSERT_EQ(mengine, NULL);
      meopts = opts;
      
      delete reporter;
      if (meopts.get_metrics_type() == "file") {
        reporter = new file_reporter("graphlab_metrics.txt");
      } else if (meopts.get_metrics_type() == "html") {
        reporter = new  html_reporter("graphlab_metrics.html");
      } else if (meopts.get_metrics_type() == "basic") {
        reporter = new basic_reporter;
      } else {
        reporter = new null_reporter;
      }
    }

    /**
     * \brief Gets the reporter
     * This function is parallel.
     */
    imetrics_reporter& get_reporter() {
      return *reporter;
    }

    /**
     * \brief Returns the engine options
     * This function is parallel
     */
    const engine_options& get_engine_options() const { 
      return meopts;
    }

    /**
     * \brief Returns a modifiable reference to the scheduler options
     * 
     * This function is parallel <b> but> any modifications to the options must be 
     * made the same way across all machines.
     */
    scheduler_options& sched_options() {
      return meopts.get_scheduler_options();
    }

    /**
     * \brief Returns a constant reference to the scheduler options
     * This function is parallel
     */
    const scheduler_options& sched_options() const{
      return meopts.get_scheduler_options();
    }


    /**
     * \brief Set the engine options by simply parsing the command line
     * arguments. 
     * The engine must not be constructed yet. 
     * All machines must call simultaneously.
     */
    bool parse_engine_options(int argc, char **argv) {
      ASSERT_EQ(mengine, NULL);
      command_line_options clopts;
      bool success = clopts.parse(argc, argv);
      ASSERT_TRUE(success);
      return set_engine_options(clopts);
    }


    /**
     * \brief Run the engine until a termination condition is reached or
     * there are no more tasks remaining to execute. This function
     * will call build_engine() internally if the engine has not yet been
     * constructed.
     * All machines must call simultaneously.
     */
    double start() {
      if (mengine == NULL) {
        bool success = build_engine();
        ASSERT_TRUE(success);
        ASSERT_NE(mengine, NULL);
      }
      // merge in options from command line and other manually set options
      mengine->set_scheduler_options( meopts.get_scheduler_options() );
      graphlab::timer ti;
      ti.start();
      mengine->start();
      return ti.current_time();
    }
  

    /**
     * \brief Add a single update function to a single vertex.
     * This function is parallel. Engine must have been constructed
     * using build_engine() prior to calling this function.
     */
    void add_task(vertex_id_t vertex,
                  typename types::update_function func,
                  double priority) {
      typename types::update_task task(vertex, func);
      add_task(task, priority);
    }


    /**
     * \brief Add a single task with a fixed priority.
     * This function is parallel. Engine must have been constructed
     * using build_engine() prior to calling this function.
     */
    void add_task(typename types::update_task task, double priority) {
      engine().add_task(task, priority);
    }

    /**
     * \brief Add the update function to all the veritces in the provided
     * vector with the given priority.
     * This function is parallel. Engine must have been constructed
     * using build_engine() prior to calling this function.
     */
    void add_tasks(const std::vector<vertex_id_t>& vertices, 
                   typename types::update_function func, double priority) {
      engine().add_tasks(vertices, func, priority);
    }


    /**
     * \brief Add the given function to all vertices using the given priority
     * This function is parallel. Engine must have been constructed
     * using build_engine() prior to calling this function.
     */
    void add_task_to_all(typename types::update_function func, 
                         double priority) {
      engine().add_task_to_all(func, priority);
    }
    
    /**
     * \brief Get the number of updates executed by the engine
     * This function is parallel. Engine must have been constructed
     * using build_engine() prior to calling this function.
     */
    size_t last_update_count() {
      ASSERT_NE(mengine, NULL);
      return mengine->last_update_count();
    }
    
    /**
     * \brief Fills the metrics with the engine options.
     * This function is parallel.
     */
    void fill_metrics() {
      coremetrics.set("ncpus", meopts.get_ncpus());
      coremetrics.set("engine", meopts.get_engine_type());
      coremetrics.set("scope", meopts.get_scope_type());
      coremetrics.set("scheduler", meopts.get_scheduler_type());
      coremetrics.set("affinities", meopts.get_cpu_affinities() ? "true" : "false");
      coremetrics.set("schedyield", meopts.get_sched_yield() ? "true" : "false");
      coremetrics.set("compile_flags", meopts.get_compile_flags());
    }

    /**
     * \brief Clears all recorded metrics. This function is parallel.
     */  
    void reset_metrics() {
      coremetrics.clear();
      if (mengine) engine().reset_metrics();
    }
      
    /**
       \brief Outputs the recorded metrics. This function is parallel.
    */
    void report_metrics() {
      coremetrics.report(get_reporter());
      engine().report_metrics(get_reporter());
    }
    
    /**
     * \brief Registers a sync with the engine.
     *
     * Registers a sync with the engine. All machines must call simultaneously.
     * 
     * The sync will be performed approximately every "interval" updates,
     * and will perform a reduction over all vertices from rangelow
     * to rangehigh inclusive.
     * The merge function may be NULL, in which it will not be used.
     * However, it is highly recommended to provide a merge function since
     * this allow the sync operation to be parallelized.
     *
     * The sync operation is guaranteed to be strictly sequentially consistent
     * with all other execution.
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
     * \param merge Combined intermediate reduction value. Required.
     * \param rangelow he lower range of vertex id to start syncing.
     *                 The range is inclusive. i.e. vertex with id 'rangelow'
     *                 and vertex with id 'rangehigh' will be included.
     *                 Defaults to 0.
     * \param rangehigh The upper range of vertex id to stop syncing.
     *                  The range is inclusive. i.e. vertex with id 'rangelow'
     *                  and vertex with id 'rangehigh' will be included.
     *                  Defaults to infinity.
     */
    void set_sync(glshared_base& shared,
                  typename types::iengine::sync_function_type sync,
                  glshared_base::apply_function_type apply,
                  const any& zero,
                  size_t sync_interval ,
                  typename types::iengine::merge_function_type merge ,
                  vertex_id_t rangelow = 0,
                  vertex_id_t rangehigh = -1) { 
      engine().set_sync(shared, sync, apply, zero, 
                        sync_interval, merge, rangelow, rangehigh);
      
    }
    

    /**
     * Performs a sync immediately. This function requires that the shared
     * variable already be registered with the engine.
     * Not implemented.
     */
    void sync_now(glshared_base& shared) ;
  private:


    distributed_control& dc;
    // graph and data objects
    typename types::graph mgraph;
    engine_options meopts;
    typename types::iengine *mengine;
    
    metrics coremetrics;
    imetrics_reporter* reporter;
  };

}
#include <graphlab/macros_undef.hpp>
#endif

