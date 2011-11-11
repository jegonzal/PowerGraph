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


// The engines
#include <graphlab/engine/iengine.hpp>
#include <graphlab/distributed2/distributed_chromatic_engine.hpp>
//#include <graphlab/distributed2/distributed_locking_engine.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/options/command_line_options.hpp>


namespace graphlab {
  
  
  /**
   *  helper for constructing graphlab engines.
   **/
  namespace distributed_engine_factory {

    template<typename Graph>
    iengine<Graph>* new_engine(distributed_control& dc,
                               const std::string& engine,
                               Graph& _graph,
                               size_t ncpus) {
     /* if(engine == "dist_locking") {
        typedef distributed_locking_engine<Graph, Scheduler> engine_type;
        return new engine_type(dc, _graph, ncpus);
      }*/ 
      if(engine == "dist_chromatic") {
        typedef distributed_chromatic_engine<Graph> engine_type;
        return new engine_type(dc, _graph, ncpus);
      } else {
        std::cout << "Invalid engine type: " << engine
                  << std::endl;
        return NULL;
      }
    } // end of new engine 

  }; // end of class engine factory
}; // End of namespace graphl




#endif

