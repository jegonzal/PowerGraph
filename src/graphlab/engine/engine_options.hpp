#ifndef ENGINE_OPTIONS_HPP
#define ENGINE_OPTIONS_HPP

#include <boost/program_options.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/engine_factory.hpp>
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
   colored} </li>

   <li> size_t splash_size: The size parameter for the splash
   scheduler. </li>
   </ul>
   */
  struct engine_options {
    //! The number of cpus
    size_t ncpus;
    //! The type of engine {async, async_sim, synchronous}
    std::string engine_type;
    //! The type of scope
    std::string scope_type;
    //! The type of scheduler to use
    std::string scheduler_type;
    //! The compiler flags
    std::string compile_flags;

    //! Use CPU affinities
    bool enable_cpu_affinities;
    bool enable_sched_yield;

    scheduler_options schedopts;
    
    engine_options() :
      ncpus(2),
      engine_type("async"),
      scope_type("edge"),
      scheduler_type("fifo"),
      enable_cpu_affinities(false),
      enable_sched_yield(true) {
      // Grab all the compiler flags 
#ifdef COMPILEFLAGS
#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)
      compile_flags = QUOTEME(COMPILEFLAGS);
#undef QUOTEME
#undef QUOTEME_
#endif
    } // end of constructor


    /**
     * create an engine for the given graph using the engine options.
     * If the engine options are not set correctly this will return
     * NULL.
     */
    template<typename Graph>
    iengine<Graph>* create_engine(Graph& graph) {
      iengine<Graph>* eng =
        engine_factory::new_engine(engine_type,
                                   scheduler_type,
                                   scope_type,
                                   graph,
                                   ncpus);
      assert(eng != NULL);
      eng->enable_sched_yield(enable_sched_yield);
      eng->enable_cpu_affinities(enable_cpu_affinities);
      eng->sched_options().merge_options(sched_options());
      return eng;
    }

    scheduler_options& sched_options() {
      return schedopts;
    }

    const scheduler_options& sched_options() const{
      return schedopts;
    }

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
                << "schedyield: " << enable_sched_yield  << std::endl;
      std::cout << "\n";
      std::cout << "Scheduler Options: \n";
      std::cout << sched_options();
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
          >> compile_flags
          >> enable_cpu_affinities
          >> enable_sched_yield;
    } // end of load
  };


  
}
#endif
