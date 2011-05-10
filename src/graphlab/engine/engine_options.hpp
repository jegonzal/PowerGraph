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

#ifndef GRAPHLAB_ENGINE_OPTIONS_HPP
#define GRAPHLAB_ENGINE_OPTIONS_HPP

#include <boost/program_options.hpp>
#include <graphlab/engine/iengine.hpp>
// #include <graphlab/engine/engine_factory.hpp>
#include <graphlab/schedulers/scheduler_options.hpp>
 
namespace graphlab {


  /**
   * The engine options class is really a simple struct that contains
   * the basic options needed to create an engine.  These options
   * include:
   
   <ul>

   <li> size_t ncpus: The number of cpus (threads) to use for this
   engine. </li>

   <li> std::string engine_type: The type of engine to use.  Currently
   we support {async, async_sim, synchronous}. </li>

   <li> std::string scope_type: The type of locking protocol (scope)
   to use. Currently we support {none, vertex, edge, full}. </li>

   <li> std::string scheduler_type: The type of scheduler to user.
   Currently we support a wide range of schedulers: {synchronous,
   fifo, priority, sampling, splash,  sweep, multiqueue_fifo,
   multiqueue_priority,  set, clustered_priority, round_robin,
   chromatic} </li>

   <li> size_t splash_size: The size parameter for the splash
   scheduler. </li>
   </ul>
   */
  class engine_options {
  public:
    //! The number of cpus
    size_t ncpus;
    //! The type of engine {async, async_sim, synchronous}
    std::string engine_type;
    //! The type of scope
    std::string scope_type;
    //! The type of scheduler to use
    std::string scheduler_type;

    //! Metrics type
    std::string metrics_type;


    //! The options associated with the scheduler
    scheduler_options scheduler_opts;

    //! The compiler flags
    std::string compile_flags;

    
    //! Use CPU affinities
    bool enable_cpu_affinities;

    bool enable_sched_yield;
 
    bool distributed_options;
    
    engine_options() :
      ncpus(2),
      engine_type("async"),
      scope_type("edge"),
      scheduler_type("fifo"),
      metrics_type("basic"),
      enable_cpu_affinities(false),
      enable_sched_yield(true),
      distributed_options(false){
      // Grab all the compiler flags 
/*#ifdef COMPILEFLAGS
#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)
      compile_flags = QUOTEME(COMPILEFLAGS);
#undef QUOTEME
#undef QUOTEME_
#endif*/
    } // end of constructor

    //! Use distributed options instead of shared memory options
    void use_distributed_options() {
      engine_type = "dist_chromatic";
      scheduler_type = "sweep";
      distributed_options = true;
    }

    //! Set the cpu affinities value (true = enabled)
    void set_cpu_affinities(bool enabled) {
      enable_cpu_affinities = enabled;
    }

    //! Get the cpu affinities value (true = enabled)
    bool get_cpu_affinities() const { 
      return enable_cpu_affinities;
    }

    //! Turn on schedule yielding  (true = enabled)
    void set_sched_yield(bool enabled) {
      enable_sched_yield = enabled;
    }

    //! Get schedule yielding value (true = enabled)
    bool get_sched_yield() const {
      return enable_sched_yield;
    }


    //! Set the number of cpus
    void set_ncpus(size_t n) { ncpus = n; }

    //! Get the number of cpus
    size_t get_ncpus() const { return ncpus; }

    //! Set the engine type
    bool set_engine_type(const std::string& etype) {
      engine_type = etype;
      return true; // TODO: ADD CHECKING
    }

    //! Get the engine type
    const std::string& get_engine_type() const {
      return engine_type;
    }

    //! Set the scope type
    bool set_scope_type(const std::string& stype) {
      scope_type = stype;
      return true; // TODO: ADD CHECKING
    }

    //! Get the scope type
    const std::string& get_scope_type() const {
      return scope_type;
    }


    //! Set the scope type
    bool set_scheduler_type(const std::string& stype) {
      scheduler_type = scheduler_opts.parse_scheduler_string(stype);
      return true; // TODO: ADD CHECKING
    }

    //! Get the scope type
    const std::string& get_scheduler_type() const {
      return scheduler_type;
    }


    //! Set the metrics type
    bool set_metrics_type(const std::string& mtype) {
      metrics_type = mtype;
      return true; // TODO: ADD CHECKING
    }

    //! Get the metrics type
    const std::string& get_metrics_type() const {
      return metrics_type;
    }
    


    //! Get the scheduler options
    const scheduler_options& get_scheduler_options() const { 
      return scheduler_opts;
    }

    //! Get the scheduler options
    scheduler_options& get_scheduler_options() { 
      return scheduler_opts;
    }

    //! Get the compiler options (flags)
    const std::string& get_compile_flags() const {
      return compile_flags;
    }




    // /**
    //  * create an engine for the given graph using the engine options.
    //  * If the engine options are not set correctly this will return
    //  * NULL.
    //  */
    // template<typename Graph>
    // iengine<Graph>* create_engine(Graph& graph) {
    //   iengine<Graph>* eng =
    //     engine_factory::new_engine(*this, graph);
    //   assert(eng != NULL);
    //   return eng;
    // }

    /**
     * Display the current engine options
     */
    virtual void print() const {
      std::cout << "GraphLab Options -------------------\n" 
                << "ncpus:       " << ncpus << "\n"
                << "engine:      " << engine_type << "\n"
                << "scope:       " << scope_type  << "\n"
                << "scheduler:   " << scheduler_type << "\n"
                << "affinities:  " << enable_cpu_affinities << "\n"
                << "metrics:     " << metrics_type << "\n"
                << "schedyield:  " << enable_sched_yield  << std::endl;
      std::cout << "\n";
      std::cout << "Scheduler Options: \n";
      std::cout << scheduler_opts;
      std::cout << std::endl;
    }



    /**
     * Save the engine options to a serialized archive
     * TODO: does not save scheduler options
     */
    void save(oarchive& arc) const {
      arc << ncpus
          << engine_type
          << scope_type
          << scheduler_type
        //          << scheduler_opts
          << compile_flags
          << enable_cpu_affinities
          << enable_sched_yield;
    } // end of save


    /**
     * Load the engine options from a serialized archive
     * TODO: does not save scheduler options
     */
    void load(iarchive& arc) {
      arc >> ncpus
          >> engine_type
          >> scope_type
          >> scheduler_type
        //          << scheduler_opts
          >> compile_flags
          >> enable_cpu_affinities
          >> enable_sched_yield;
    } // end of load
  };


  
}
#endif
