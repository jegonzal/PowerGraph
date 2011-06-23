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


#ifndef GRAPHLAB_DISTRIBUTED_ENGINE_FACTORY_HPP
#define GRAPHLAB_DISTRIBUTED_ENGINE_FACTORY_HPP

#include <string>
#include <sstream>
#include <graphlab/rpc/dc.hpp>

// The graph is needed to build an engine
#include <graphlab/distributed2/graph/distributed_graph.hpp>

// The engines
#include <graphlab/engine/iengine.hpp>
#include <graphlab/distributed2/distributed_chromatic_engine.hpp>
#include <graphlab/distributed2/distributed_locking_engine.hpp>
#include <graphlab/engine/engine_options.hpp>



// Scopes
#include <graphlab/scope/general_scope_factory.hpp>


// Schedulers
#include <graphlab/distributed2/distributed_scheduler_list.hpp>

#include <boost/preprocessor.hpp>


namespace graphlab {
  
  
  /**
   *  helper for constructing graphlab engines.
   **/
  namespace distributed_engine_factory {

    template<typename Graph, typename Scheduler, typename ScopeFactory>
    iengine<Graph>* new_engine(distributed_control& dc,
                               const std::string& engine,
                               Graph& _graph,
                               size_t ncpus) {
      if(engine == "dist_locking") {
        typedef distributed_locking_engine<Graph, Scheduler> engine_type;
        return new engine_type(dc, _graph, ncpus);
      } if(engine == "dist_chromatic") {
        typedef distributed_chromatic_engine<Graph> engine_type;
        return new engine_type(dc, _graph, ncpus);
      } else {
        std::cout << "Invalid engine type: " << engine
                  << std::endl;
        return NULL;
      }
    } // end of new engine 

    
    template<typename Graph, typename Scheduler>
    iengine<Graph>* new_engine(distributed_control& dc,
                               const std::string& engine,
                               const std::string& scope_factory,
                               Graph& _graph,
                               size_t ncpus) {
      iengine<Graph>* ret = new_engine<Graph, Scheduler,
        general_scope_factory<Graph> >(dc, engine, _graph, ncpus);
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
    iengine<Graph>* new_engine(distributed_control& dc,
                               const engine_options& eopts,
                               Graph& _graph) {
      iengine<Graph>* eng = NULL;

#define __GENERATE_NEW_ENGINE__(r_unused, data_unused, i,  elem)        \
      BOOST_PP_EXPR_IF(i, else)                                         \
        if (eopts.get_scheduler_type() == BOOST_PP_TUPLE_ELEM(3,0,elem)) { \
          eng = new_engine<Graph, BOOST_PP_TUPLE_ELEM(3,1,elem) <Graph> > \
            ( dc,                                                       \
              eopts.get_engine_type(),                                  \
              eopts.get_scope_type(),                                   \
              _graph,                                                   \
              eopts.get_ncpus() );                                      \
        }
      
      // generate the construction calls
      BOOST_PP_SEQ_FOR_EACH_I(__GENERATE_NEW_ENGINE__, _, __DISTRIBUTED_SCHEDULER_LIST__)
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
    iengine<Graph>* new_engine(distributed_control &dc,
                               const std::string& engine_type,
                               const std::string& scheduler_type,
                               const std::string& scope_type,
                               Graph& _graph,
                               size_t ncpus) {
      engine_options eopts;
      eopts.set_engine_type(engine_type);
      eopts.set_scheduler_type(scheduler_type);
      eopts.set_scope_type(scope_type);
      eopts.set_ncpus(ncpus);
      return new_engine(dc, eopts, _graph);
    }



  }; // end of class engine factory
}; // End of namespace graphl




#endif

