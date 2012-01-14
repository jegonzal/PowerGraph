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


#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/options/command_line_options.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/engine/shared_memory_engine.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {

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
   
     The core contains the 
   
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
  template <typename Graph, typename UpdateFunctor>
  class core {
  public:
    typedef Graph graph_type;
    typedef UpdateFunctor update_functor_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef typename graph_type::edge_id_type   edge_id_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef shared_memory_engine<graph_type, update_functor_type>
    engine_type;

  private:
    
    // graph and data objects
    graph_type mgraph;
    engine_type mengine;



  public:
    /// default constructor does nothing
    core() :  mengine(mgraph) { } 
      // coremetrics("core"), reporter(new null_reporter) { }
  private:
    //! Core is not copyable
    core(const core& other);
    //! Core is not copyable
    core& operator=(const core& other);



  public:


    // ~core() { 
    //   // if (opts.get_metrics_type() != "none") {        
    //   //   // Write options to metrics
    //   //   fill_metrics();
    //   //   report_metrics();
    //   // }
    //   //      delete reporter;
    // } 
       
    /// \brief Get a modifiable reference to the graph associated with this core
    graph_type& graph() { return mgraph; }

    /// \brief Get a constant reference to the graph associated with this core
    const graph_type& graph() const { return mgraph; }

    /**
     * \brief Set the type of scheduler.
     *
     * This will destroy the current engine and any tasks currently
     * associated with the scheduler.  See \ref Schedulers for the
     * list of supported schedulers.
     */
    void set_scheduler_type(const std::string& scheduler_str) {
      graphlab_options opts = mengine.get_options();
      opts.set_scheduler_type(scheduler_str);
      mengine.set_options(opts);
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
    void set_scope_type(const std::string& scope_str) {
      graphlab_options opts = mengine.get_options();
      opts.set_scope_type(scope_str);
      mengine.set_options(opts);
    }


    // /**
    //  * \brief Set the engine type.
    //  *
    //  * This will destroy the current engine and any tasks associated
    //  * with the current scheduler. 
    //  *
    //  *  \li \b "async" This is the regular multithreaded engine
    //  *  \li \b "async_sim" This is a single threaded engine. But it can be 
    //  *                     be started with multiple "simulated threads".
    //  *                     The simulation is low-fidelity however, and should
    //  *                     be used with caution.
    //  */
    // void set_engine_type(const std::string& engine_type) {
    //   check_engine_modification();
    //   bool success = opts.set_engine_type(engine_type);
    //   ASSERT_TRUE(success);
    //   destroy_engine();
    // }

    
    // /**
    //  * \brief Sets the output format of any recorded metrics
    //  *  \li \b "none" No reporting
    //  *  \li \b "basic" Outputs to screen
    //  *  \li \b "file" Outputs to a text file graphlab_metrics.txt
    //  *  \li \b "html" Outputs to a html file graphlab_metrics.html
    //  */
    // void set_metrics_type(const std::string& metrics_type) {
    //   bool metrics_set_success = opts.set_metrics_type(metrics_type);
    //   ASSERT_TRUE(metrics_set_success);
      
    //   delete reporter;
    //   if (opts.get_metrics_type() == "file") {
    //     reporter = new file_reporter("graphlab_metrics.txt");
    //   } else if (opts.get_metrics_type() == "html") {
    //     reporter = new  html_reporter("graphlab_metrics.html");
    //   } else if (opts.get_metrics_type() == "basic") {
    //     reporter = new basic_reporter;
    //   } else {
    //     reporter = new null_reporter;
    //   }
    // }

    /**
       \brief Destroys a created engine (if any).
    */
    void reset() {  mengine.reset(); }
    
    /**
     * \brief Set the number of cpus that the engine will use.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler. 
     *
     */
    void set_ncpus(size_t ncpus) {
      graphlab_options opts = mengine.get_options();
      opts.set_ncpus(ncpus);
      mengine.set_options(opts);
    }


    /**
     * Get a reference to the active engine.  If no engine exists one is
     * created.
     */
    engine_type& engine() { return mengine; }




    /**
     * \brief Set the engine options by passing in an engine options object.
     */
    void set_options(const graphlab_options& opts) {
      mengine.set_options(opts);
    }

    // imetrics_reporter& get_reporter() { return *reporter; }

    /**
     * \brief Returns the engine options
     */
    const graphlab_options& get_options() const { 
      return mengine.get_options();
    }

    /**
     * \brief Set the engine options by simply parsing the command line
     * arguments. 
     */
    void parse_options(int argc, char **argv) {
      command_line_options clopts;
      bool success = clopts.parse(argc, argv);
      ASSERT_TRUE(success);
      return set_options(clopts);
    }


    /**
     * \brief Run the engine until a termination condition is reached or
     * there are no more tasks remaining to execute.
     */
    double start() {
      graphlab::timer ti;
      ti.start();
      mengine.start();
      return ti.current_time();
    }
  

    /**
     * \brief Add a single update function to a single vertex.
     */
    void schedule(vertex_id_type vid, const update_functor_type& fun) {
      mengine.schedule(vid, fun);
    }

    /**
     * \brief Add an update function to a vector of vertices
     */
    void schedule(const std::vector<vertex_id_type>& vid,
                  const update_functor_type& fun) {
      mengine.schedule(vid, fun);
    }

    
    

    /**
     * \brief Add the given function to all vertices using the given priority
     */
    void schedule_all(const update_functor_type& fun) {
      mengine.schedule_all(fun);
    }
    
    /**
     * \brief Get the number of updates executed by the engine
     */
    size_t last_update_count() {
      return mengine.last_update_count();
    }
    
    // void fill_metrics() {
    //   coremetrics.set("ncpus", opts.get_ncpus());
    //   coremetrics.set("engine", opts.get_engine_type());
    //   coremetrics.set("scope", opts.get_scope_type());
    //   coremetrics.set("scheduler", opts.get_scheduler_type());
    //   coremetrics.set("affinities", opts.get_cpu_affinities() ? "true" : "false");
    //   coremetrics.set("schedyield", opts.get_sched_yield() ? "true" : "false");
    //   coremetrics.set("compile_flags", opts.get_compile_flags());
    // }

    // void reset_metrics() {
    //   coremetrics.clear();
    //   engine().reset_metrics();
    // }
      
    // /**
    //    \brief Outputs the recorded metrics
    // */
    // void report_metrics() {
    //   coremetrics.report(get_reporter());
    //   engine().report_metrics(get_reporter());
    // }
    

    //! Add a global entry 
    template< typename T >
    void add_global(const std::string& key, const T& value, size_t size = 1) {
      engine().add_global(key, value, size); 
    }

    //! Add a global entry 
    template< typename T >
    void add_global_const(const std::string& key, const T& value, size_t size = 1) {
      engine().add_global_const(key, value, size); 
    }

    //! Change the value of a global entry
    template< typename T >
    void set_global(const std::string& key, const T& value, size_t index = 0) {
      engine().set_global(key, value, index);
    }

    //! Get a copy of the value of a global entry
    template< typename T >
    T get_global(const std::string& key, size_t index = 0) {
      return engine().get_global<T>(key, index);
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
    template<typename Accum>
    void add_sync(const std::string& key,           
                  const Accum& zero,                 
                  size_t sync_interval,
                  bool use_barrier = false,
                  vertex_id_type begin_vid = 0,
                  vertex_id_type end_vid = 
                  std::numeric_limits<vertex_id_type>::max()) {
      engine().add_sync(key, zero, sync_interval,
                        use_barrier, begin_vid, end_vid);
    }    

    /**
     * Performs a sync immediately. This function requires that the shared
     * variable already be registered with the engine.
     */
    void sync_now(const std::string& key) { 
      engine().sync_now(key);
    };





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
      arc << mgraph << mengine.get_options();
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
      graphlab_options opts;
      arc >> mgraph >> opts;
      mengine.set_options(opts);
    } // end of load


    // metrics coremetrics;
    // imetrics_reporter* reporter;
  };

}
#include <graphlab/macros_undef.hpp>
#endif

