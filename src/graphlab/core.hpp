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


#ifndef GRAPHLAB_CORE_HPP
#define GRAPHLAB_CORE_HPP

#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/engine_options.hpp>
#include <graphlab/engine/engine_factory.hpp>

#include <graphlab/util/command_line_options.hpp>

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
     
     Because many GraphLab programs will consists of a graph and an
     engine we have created a single data-structure, called a core,
     which manages all the pieces of GraphLab including engine and
     scheduler construction parameters.

     The core is templatized over the VertexType and EdgeType however
     by using the ref types typedef, one can simply create a core by
     doing the following:
   
     \code
     gl::core glcore;
     \endcode
   
cd     The core contains the 
   
     \li Data Graph: which represents the structured data dependencies.
     \li Engine: The computational structure which contains the
     scheduling and execution statistics for the GraphLab program. The
     core provides pass-through calls for many engine functions.
        
     The core also manages the engine and scheduler construction
     parameters.
   
     The core will invisibly recreate the engine each time engine
     options are modified. This will mean that this internal behavior of
     the core should be pretty much "transparent" for the typical use
     case where engine options and scheduler options are defined before
     tasks are added to the scheduler.
   
     Otherwise, modifications to the engine options will result in the
     clearing of all scheduler tasks.
  */
  template <typename VertexType, typename EdgeType>
  class core {
  public:
    typedef graphlab::types<graphlab::graph<VertexType, EdgeType> > types;
    typedef typename types::vertex_id vertex_id_type;
    typedef typename types::edge_id edge_id_type;
    typedef typename types::edge_list edge_list_type;

  public:
    /// default constructor does nothing
    core() : 
      mengine(NULL),
      engine_has_been_modified(false), 
      coremetrics("core"), reporter(new null_reporter) { }
  private:
    //! Core is not copyable
    core(const core& other);
    //! Core is not copyable
    core& operator=(const core& other);

  public:


    ~core() { 
      if (meopts.get_metrics_type() != "none") {        
        // Write options to metrics
        fill_metrics();
        report_metrics();
      }
      destroy_engine(); 
      delete reporter;
    } 
       
    /// \brief Get a modifiable reference to the graph associated with this core
    typename types::graph& graph() { return mgraph; }

    /// \brief Get a constant reference to the graph associated with this core
    const typename types::graph& graph() const { return mgraph; }

    /**
     * \brief Set the type of scheduler.
     *
     * This will destroy the current engine and any tasks currently
     * associated with the scheduler.  See \ref Schedulers for the
     * list of supported schedulers.
     */
    void set_scheduler_type(const std::string& scheduler_type) {
      check_engine_modification();
      bool success = meopts.set_scheduler_type(scheduler_type);
      ASSERT_TRUE(success);
      destroy_engine();
    }

    /**
     * \brief Set the scope consistency model used in this engine.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler.  The available scopes are:
     * 
     *  \li \b "full" This ensures full data consistency within the scope
     *  \li \b "edge" This ensures data consistency with just the
     *     vertex and edges
     *  \li \b "vertex" This ensures that a vertex cannot be updated
     *     by two processors simultaneously
     *  \li \b "none" This eliminates all locking 
     *
     * See \ref Scopes for details
     */
    void set_scope_type(const std::string& scope_type) {
      check_engine_modification();
      bool success = meopts.set_scope_type(scope_type);
      ASSERT_TRUE(success);
      destroy_engine();
    }


    /**
     * \brief Set the engine type.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler. 
     *
     *  \li \b "async" This is the regular multithreaded engine
     *  \li \b "async_sim" This is a single threaded engine. But it can be 
     *                     be started with multiple "simulated threads".
     *                     The simulation is low-fidelity however, and should
     *                     be used with caution.
     */
    void set_engine_type(const std::string& engine_type) {
      check_engine_modification();
      bool success = meopts.set_engine_type(engine_type);
      ASSERT_TRUE(success);
      destroy_engine();
    }
    
    /**
     * \brief Sets the output format of any recorded metrics
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
       \brief Destroys a created engine (if any).
    */
    void reset() {
      engine_has_been_modified = false;
      destroy_engine();
    }
    
    /**
     * \brief Set the number of cpus that the engine will use.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler. 
     *
     */
    void set_ncpus(size_t ncpus) {
      check_engine_modification();
      meopts.set_ncpus(ncpus);
      destroy_engine();
    }


    /**
     * Get a reference to the active engine.  If no engine exists one is
     * created.
     */
    typename types::iengine& engine() { 
      bool engine_build_success = auto_build_engine();
      ASSERT_TRUE(engine_build_success);
      return *mengine; 
    }




    /**
     * \brief Destroys and reconstructs the current engine,
     * reprocessing the engine arguments.  
     */
    bool rebuild_engine() {
      destroy_engine();
      ASSERT_EQ(mengine, NULL);
      return auto_build_engine();
    }

    /**
     * \brief Set the engine options by passing in an engine options object.
     */
    void set_engine_options(const engine_options& opts) {
      check_engine_modification();
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

    imetrics_reporter& get_reporter() {
      return *reporter;
    }

    /**
     * \brief Returns the engine options
     */
    const engine_options& get_engine_options() const { 
      return meopts;
    }

    /**
     * \brief Returns a modifiable reference to the scheduler options
     */
    scheduler_options& sched_options() {
      return meopts.get_scheduler_options();
    }

    /**
     * \brief Returns a constant reference to the scheduler options
     */
    const scheduler_options& sched_options() const{
      return meopts.get_scheduler_options();
    }


    /**
     * \brief Set the engine options by simply parsing the command line
     * arguments. 
     */
    bool parse_engine_options(int argc, char **argv) {
      check_engine_modification();
      command_line_options clopts;
      bool success = clopts.parse(argc, argv);
      ASSERT_TRUE(success);
      return set_engine_options(clopts);
    }


    /**
     * \brief Run the engine until a termination condition is reached or
     * there are no more tasks remaining to execute.
     */
    double start() {
      bool success = auto_build_engine();
      ASSERT_TRUE(success);
      ASSERT_NE(mengine, NULL);
      // merge in options from command line and other manually set options
      mengine->set_scheduler_options( meopts.get_scheduler_options() );
      graphlab::timer ti;
      ti.start();
      mengine->start();
      return ti.current_time();
    }
  

    /**
     * \brief Add a single update function to a single vertex.
     */
    void add_task(vertex_id_type vertex,
                  typename types::update_function func,
                  double priority) {
      engine_has_been_modified = true;
      typename types::update_task task(vertex, func);
      add_task(task, priority);
    }


    /**
     * \brief Add a single task with a fixed priority.
     */
    void add_task(typename types::update_task task, double priority) {
      engine_has_been_modified = true;
      engine().add_task(task, priority);
    }

    /**
     * \brief Add the update function to all the veritces in the provided
     * vector with the given priority.
     */
    void add_tasks(const std::vector<vertex_id_type>& vertices, 
                   typename types::update_function func, double priority) {
      engine_has_been_modified = true;
      engine().add_tasks(vertices, func, priority);
    }


    /**
     * \brief Add the given function to all vertices using the given priority
     */
    void add_task_to_all(typename types::update_function func, 
                         double priority) {
      engine_has_been_modified = true;
      engine().add_task_to_all(func, priority);
    }
    
    /**
     * \brief Get the number of updates executed by the engine
     */
    size_t last_update_count() {
      if(mengine == NULL) return 0;
      else return mengine->last_update_count();
    }
    
    void fill_metrics() {
      coremetrics.set("ncpus", meopts.get_ncpus());
      coremetrics.set("engine", meopts.get_engine_type());
      coremetrics.set("scope", meopts.get_scope_type());
      coremetrics.set("scheduler", meopts.get_scheduler_type());
      coremetrics.set("affinities", meopts.get_cpu_affinities() ? "true" : "false");
      coremetrics.set("schedyield", meopts.get_sched_yield() ? "true" : "false");
      coremetrics.set("compile_flags", meopts.get_compile_flags());
    }

    void reset_metrics() {
      coremetrics.clear();
      engine().reset_metrics();
    }
      
    /**
       \brief Outputs the recorded metrics
    */
    void report_metrics() {
      coremetrics.report(get_reporter());
      engine().report_metrics(get_reporter());
    }
    
    /**
     * \brief Registers a sync with the engine.
     *
     * Registers a sync with the engine.
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
    void set_sync(glshared_base& shared,
                  typename types::iengine::sync_function_type sync,
                  glshared_base::apply_function_type apply,
                  const any& zero,
                  size_t sync_interval = 0,
                  typename types::iengine::merge_function_type merge = NULL,
                  vertex_id_type rangelow = 0,
                  vertex_id_type rangehigh = -1) { 
      engine_has_been_modified = true;
      engine().set_sync(shared, sync, apply, zero, 
                        sync_interval, merge, rangelow, rangehigh);
      
    }
    

    /**
     * Performs a sync immediately. This function requires that the shared
     * variable already be registered with the engine.
     */
    void sync_now(glshared_base& shared) { 
      engine().sync_now(shared);
    };
  private:

    /**
     * Build the engine if it has not already been built.
     */
    bool auto_build_engine() {
      if(mengine == NULL) {
        // create the engine
        mengine = engine_factory::new_engine(meopts, mgraph);
        if(mengine == NULL) return false;
        mengine->set_engine_options(meopts.get_engine_options());
      }
      // scheduler options is one parameter that is allowed
      // to change without rebuilding the engine
      return true;
    }

    /**
     * Destroy the engine if one exists.
     */
    void destroy_engine() {
      if(mengine != NULL) {
        delete mengine;
        mengine = NULL;
      }
    }



    /** Save the core to a file */
    void save(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      ASSERT_TRUE(fout.good());
      oarchive oarc(fout);
      oarc << *this;
      fout.close();
    } // end of save
    
    /** Save the core to an archive */
    void save(oarchive& arc) const {
      arc << mgraph
          << meopts;
    } // end of save


    /** Load the core from a file. */
    void load(const std::string& filename) {
      std::ifstream fin(filename.c_str());
      ASSERT_TRUE(fin.good());
      iarchive iarc(fin);
      iarc >> *this;
      fin.close();
    } // end of load


    /** Load the core from an archive. */
    void load(iarchive& arc) {
      arc >> mgraph
          >> meopts;
    } // end of load


    void check_engine_modification() {
      ASSERT_MSG(engine_has_been_modified == false, 
                 "Modifications to the engine/scheduler parameters are not"
                 "allowed once tasks have been inserted into the engine.");
    }
    
    // graph and data objects
    typename types::graph mgraph;
    engine_options meopts;
    typename types::iengine *mengine;
    /** For error tracking. Once engine has been modified, any scheduler/
     * engine parameter modifications will reset the modifications
     */
    bool engine_has_been_modified;
    metrics coremetrics;

    imetrics_reporter* reporter;
  };

}
#include <graphlab/macros_undef.hpp>
#endif

