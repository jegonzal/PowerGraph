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

#ifndef GRAPHLAB_ENGINE_FACTORY_HPP
#define GRAPHLAB_ENGINE_FACTORY_HPP

#include <string>
#include <sstream>
// The graph is needed to build an engine
#include <graphlab/graph/graph.hpp>

// The engines
#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/asynchronous_engine.hpp>
#include <graphlab/engine/engine_options.hpp>



// Scopes
#include <graphlab/scope/general_scope_factory.hpp>


// Schedulers
#include <graphlab/schedulers/scheduler_list.hpp>

#include <boost/preprocessor.hpp>


namespace graphlab {
  
  
  /**
   *  helper for constructing graphlab engines.
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
     * Allocate an engine given engine options.
     */
    template<typename Graph>
    iengine<Graph>* new_engine(const engine_options& eopts,
                               Graph& _graph) {
      iengine<Graph>* eng = NULL;

#define __GENERATE_NEW_ENGINE__(r_unused, data_unused, i,  elem)        \
      BOOST_PP_EXPR_IF(i, else)                                         \
        if (eopts.get_scheduler_type() == BOOST_PP_TUPLE_ELEM(3,0,elem)) { \
          eng = new_engine<Graph, BOOST_PP_TUPLE_ELEM(3,1,elem) <Graph> > \
            ( eopts.get_engine_type(),                                  \
              eopts.get_scope_type(),                                   \
              _graph,                                                   \
              eopts.get_ncpus() );                                      \
        }
      
      // generate the construction calls
      BOOST_PP_SEQ_FOR_EACH_I(__GENERATE_NEW_ENGINE__, _, __SCHEDULER_LIST__)
        /*
          if(scheduler == "fifo") {
          eng =  new_engine<Graph, fifo_scheduler<Graph> >(engine,
          scope_factory,
          _graph,
          ncpus);
          } else if ....*/
      else {
        std::cout << "Invalid scheduler type: " << eopts.get_scheduler_type()
                  << std::endl;
        return NULL;
      }

      if(eng != NULL) {
        // Should we merge instead?
        eng->set_scheduler_options( eopts.get_scheduler_options() );
      }
      return eng;
#undef __GENERATE_NEW_ENGINE__ 
      
    } // end of new engine  





    /**
     * Deprecated but used by some older apps
     *
     * Allocate an engine given the strings for the engine type, scope
     * factory, and scheduler.
     *
     * \param engine  Type of engine to construct. {async, async_sim, synchronous}
     * \param scope    Type of scope to use.  {none, vertex, edge, full}
     * \param scheduler Type of scheduler to use synchronous, fifo, priority, sampling,
     *                 sweep, multiqueue_fifo, multiqueue_priority,
     *                 clustered_priority({metis, bfs, random}, verts. per part)
     *                 round_robin, chromatic}
     * 
     * Note that the caller is responsible for freeing the
     * corresponding engine
     * \anchor gl_new_engine
     */
    template<typename Graph>
    iengine<Graph>* new_engine(const std::string& engine_type,
                               const std::string& scheduler_type,
                               const std::string& scope_type,
                               Graph& _graph,
                               size_t ncpus) {
      engine_options eopts;
      eopts.set_engine_type(engine_type);
      eopts.set_scheduler_type(scheduler_type);
      eopts.set_scope_type(scope_type);
      eopts.set_ncpus(ncpus);
      return new_engine(eopts, _graph);
    }



  }; // end of class engine factory
}; // End of namespace graphl




#endif

