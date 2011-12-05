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


#ifndef GRAPHLAB_DISTRIBUTED_DELTA_CORE_HPP
#define GRAPHLAB_DISTRIBUTED_DELTA_CORE_HPP

#include <graphlab/engine/iengine.hpp>
#include <graphlab/distributed2/distributed2_includes.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/graph/graph.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
     \brief A GraphLab core is the base (or core) data structure in GraphLab.
     
     This is like \ref graphlab::core but for the distributed setting.

     The core is templatized over the VertexType and EdgeType however
     by using the ref types typedef, one can simply create a core by
     doing the following:
   
     \code
     gl::distributed_delta_core glcore;
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
  template <typename Graph, typename UpdateFunctor>
  class distributed_delta_core {
  public:
    typedef Graph graph_type;
    typedef UpdateFunctor update_functor_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef typename graph_type::edge_id_type   edge_id_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef distributed_delta_engine<graph_type, update_functor_type>
          engine_type;

    
  public:
    /** default constructor. Graph is constructed using the atom index.
     * All machines must construct simultaneously.
    */
    distributed_delta_core(distributed_control &dc, std::string atomindex,
                     disk_graph_atom_type::atom_type atomtype = disk_graph_atom_type::WRITE_ONLY_ATOM) :
      dc(dc),
      mgraph(dc, atomindex, false, false, atomtype),
      mengine(NULL){ }
  private:
    //! Core is not copyable
    distributed_delta_core(const distributed_delta_core& other);
    //! Core is not copyable
    distributed_delta_core& operator=(const distributed_delta_core& other);

  public:

    /**
     * Destructor. 
     * All machines must call simultaneously.
     */
    ~distributed_delta_core() { 
    } 
       
    /** Get a modifiable reference to the graph associated with this core
     * This function is parallel.
     */
    graph_type& graph() { return mgraph; }

    /** Get a constant reference to the graph associated with this core
     * This function is parallel.
     */
    const graph_type& graph() const { return mgraph; }



    /**
     * \brief Set the type of scheduler.
     *
     * This will destroy the current engine and any tasks currently
     * associated with the scheduler.  See \ref Schedulers for the
     * list of supported schedulers.
     */
    void set_scheduler_type(const std::string& scheduler_str) {
      opts.set_scheduler_type(scheduler_str);
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
      opts.set_scope_type(scope_str);
    }

    /**
       \brief Destroys a created engine (if any).
    */
    void reset() {  engine().reset(); }
    
    /**
     * \brief Set the number of cpus that the engine will use.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler. 
     *
     */
    void set_ncpus(size_t ncpus) {
      opts.set_ncpus(ncpus);
    }


    void build_engine() {
      mengine = new engine_type(dc, mgraph, opts.get_ncpus());
      mengine->set_options(opts);
    }

    /**
     * Get a reference to the active engine.  If no engine exists one is
     * created.
     */
    engine_type& engine() { 
      ASSERT_TRUE(mengine != NULL);
      return *mengine; 
    }




    /**
     * \brief Set the engine options by passing in an engine options object.
     */
    void set_options(const graphlab_options& newopts) {
      opts = newopts;
    }

    // imetrics_reporter& get_reporter() { return *reporter; }

    /**
     * \brief Returns the engine options
     */
    const graphlab_options& get_options() const { 
      return opts;
    }

    /**
     * \brief Set the engine options by simply parsing the command line
     * arguments. 
     */
    bool parse_options(int argc, char **argv) {
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
      engine().start();
      return ti.current_time();
    }
  

    /**
     * \brief Add a single update function to a single vertex.
     */
    void schedule(vertex_id_type vid, const update_functor_type& fun) {
      engine().schedule(vid, fun);
    }

    /**
     * \brief Add an update function to a vector of vertices
     */
    void schedule(const std::vector<vertex_id_type>& vid,
                  const update_functor_type& fun) {
      engine().schedule(vid, fun);
    }

    
    

    /**
     * \brief Add the given function to all vertices using the given priority
     */
    void schedule_all(const update_functor_type& fun) {
      engine().schedule_all(fun);
    }
    
    /**
     * \brief Get the number of updates executed by the engine
     */
    size_t last_update_count() {
      return engine().last_update_count();
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


  private:


    distributed_control& dc;
    // graph and data objects
    graph_type mgraph;
    engine_type *mengine;
    graphlab_options opts;
  };

}
#include <graphlab/macros_undef.hpp>
#endif

