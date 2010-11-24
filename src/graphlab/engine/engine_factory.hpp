
#ifndef GRAPHLAB_ENGINE_FACTORY_HPP
#define GRAPHLAB_ENGINE_FACTORY_HPP

#include <string>
#include <sstream>
// The graph is needed to build an engine
#include <graphlab/graph/graph.hpp>

// The engines
#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/asynchronous_engine.hpp>



// Scopes
#include <graphlab/scope/general_scope_factory.hpp>


// Schedulers
#include <graphlab/schedulers/fifo_scheduler.hpp>
#include <graphlab/schedulers/priority_scheduler.hpp>
#include <graphlab/schedulers/sampling_scheduler.hpp>
#include <graphlab/schedulers/round_robin_scheduler.hpp>
#include <graphlab/schedulers/colored_scheduler.hpp>
#include <graphlab/schedulers/sweep_scheduler.hpp>
#include <graphlab/schedulers/splash_scheduler.hpp>
#include <graphlab/schedulers/multiqueue_fifo_scheduler.hpp>
#include <graphlab/schedulers/multiqueue_priority_scheduler.hpp>
#include <graphlab/schedulers/clustered_priority_scheduler.hpp>

namespace graphlab {
  
  
  /**
   * Class helper for constructing graphlab engines.
   **/
  
  namespace engine_factory {

    template<typename Graph, typename Scheduler, typename ScopeFactory>
    iengine<Graph>* new_engine(const std::string& engine,
                               Graph& _graph,
                               size_t ncpus) {
      if(engine == "async") {
        typedef asynchronous_engine<Graph, Scheduler, ScopeFactory> engine_type;
        return new engine_type(_graph, ncpus, engine_type::THREADED);
      } if(engine == "async_sim") {
        typedef asynchronous_engine<Graph, Scheduler, ScopeFactory> engine_type;
        return new engine_type(_graph, ncpus, engine_type::SIMULATED);
      } else {
        std::cout << "Invalid engine type: " << engine
                  << std::endl;
        return NULL;
      }
    } // end of new engine 

    
    template<typename Graph, typename Scheduler>
    iengine<Graph>* new_engine(const std::string& engine,
                               const std::string& scope_factory,
                               Graph& _graph,
                               size_t ncpus) {
      iengine<Graph>* ret = new_engine<Graph, Scheduler,
                          general_scope_factory<Graph> >(engine,
                                                         _graph,
                                                         ncpus);
      if(ret == NULL) { return NULL; }
      
      if(scope_factory == "unsync" || scope_factory == "vertex") {
        ret->set_default_scope(scope_range::VERTEX_CONSISTENCY);
      } else if(scope_factory == "locked" || scope_factory == "edge") {
        ret->set_default_scope(scope_range::EDGE_CONSISTENCY);
      } else if(scope_factory == "extlocked" || scope_factory == "full") {
        ret->set_default_scope(scope_range::FULL_CONSISTENCY);
      } else if(scope_factory == "null" || scope_factory == "none") {
        ret->set_default_scope(scope_range::NULL_CONSISTENCY);
      } else {
        std::cout << "Invalid scope type: " << scope_factory
                  << std::endl;
        return NULL;
      }
      
      if(scope_factory == "unsync" || scope_factory == "locked" ||
        scope_factory == "extlocked") {
        logstream(LOG_WARNING) << "The scope names \"unsync\", \"locked\""
                                  "and \"extlocked\" have been deprecated";
        logstream(LOG_WARNING) << "Please use the new names \"vertex\", "
                                  "\"edge\" and \"full\"" << std::endl;
      }
      return ret;
    }

  

    /**
     * Allocate an engine given the strings for the engine type, scope
     * factory, and scheduler.
     *
     * engine       = {async, async_sim, synchronous}
     * scope        = {none, vertex, edge, full}
     * scheduler    = {synchronous, fifo, priority, sampling,
     *                 sweep, multiqueue_fifo, multiqueue_priority,
     *                 clustered_priority({metis, bfs, random}, verts. per part)
     *                 round_robin, colored}
     * 
     * Note that the caller is responsible for freeing the
     * corresponding engine
     */
    template<typename Graph>
    iengine<Graph>* new_engine(const std::string& engine,
                               const std::string& scheduler_raw,
                               const std::string& scope_factory,
                               Graph& _graph,
                               size_t ncpus = 1) {
      // Break the scheduler string appart
      size_t first_paren = scheduler_raw.find_first_of('(');
      size_t last_paren = scheduler_raw.find_last_of(')');
      std::string scheduler = scheduler_raw.substr(0, first_paren);
      std::string arguments;
      // Fill in the arguments if such are possibe
      if(first_paren != std::string::npos &&
         last_paren != std::string::npos) {
        arguments = scheduler_raw.substr(first_paren + 1,
                                         last_paren - first_paren - 1 );
      }
      
      if(!arguments.empty()) {
        std::replace(arguments.begin(), arguments.end(), ',', ' ');
        std::replace(arguments.begin(), arguments.end(), ';', ' ');        
      }     
      std::stringstream arg_strm(arguments);
      
      iengine<Graph>* eng = NULL;
      if(scheduler == "fifo") {
        eng =  new_engine<Graph, fifo_scheduler<Graph> >(engine,
                                                  scope_factory,
                                                  _graph,
                                                  ncpus);
      } else if(scheduler == "priority") {
        eng = new_engine<Graph, priority_scheduler<Graph> >(engine,
                                                      scope_factory,
                                                      _graph,
                                                      ncpus);
      } else if(scheduler == "sampling") {
        eng = new_engine<Graph, sampling_scheduler<Graph> >(engine,
                                                      scope_factory,
                                                      _graph,
                                                      ncpus);
      } else if(scheduler == "sweep") {
        eng = new_engine<Graph, sweep_scheduler<Graph> >(engine,
                                                          scope_factory,
                                                          _graph,
                                                          ncpus);
      } else if(scheduler == "multiqueue_fifo") {
        eng = new_engine<Graph, multiqueue_fifo_scheduler<Graph> >(engine,
                                                             scope_factory,
                                                             _graph,
                                                             ncpus);     

      } else if(scheduler == "multiqueue_priority") {
        eng = new_engine<Graph, multiqueue_priority_scheduler<Graph> >(engine,
                                                                        scope_factory,
                                                                        _graph,
                                                                        ncpus);     
      } else if(scheduler == "clustered_priority") {
        eng = new_engine<Graph, clustered_priority_scheduler<Graph> >(engine,
                                                                  scope_factory,
                                                                  _graph,
                                                                  ncpus);
      } else if(scheduler == "round_robin") {
        eng = new_engine<Graph, round_robin_scheduler<Graph> >(engine,
                                                                scope_factory,
                                                                _graph,
                                                                ncpus);      
      } else if(scheduler == "splash") {
        eng = new_engine<Graph, splash_scheduler<Graph> >(engine,
                                                          scope_factory,
                                                          _graph,
                                                          ncpus);            
      } else if(scheduler == "colored") {
        eng = new_engine<Graph, colored_scheduler<Graph> >
                                                    (engine,
                                                    scope_factory,
                                                    _graph,
                                                    ncpus);
      } else {
        std::cout << "Invalid scheduler type: " << scheduler
                  << std::endl;
        return NULL;
      }
      if(eng != NULL) {
        if(!arguments.empty()) {
          eng->set_sched_option(arg_strm);
        }
      }
      return eng;
      
      
    } // end of new engine  

  }; // end of class engine factory
}; // End of namespace graphl




#endif
