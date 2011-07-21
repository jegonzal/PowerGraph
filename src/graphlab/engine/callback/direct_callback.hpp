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


/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */



#ifndef GRAPHLAB_DIRECT_CALLBACK_HPP
#define GRAPHLAB_DIRECT_CALLBACK_HPP

#include <vector>
#include <algorithm>

#include <graphlab/graph/graph.hpp>
#include <graphlab/engine/icallback.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {

  /**
     This is a callback class that is passed to the update functions,
     and allow the update functions to add tasks back into the
     scheduler \see ischeduler_callback 
  */
  template<typename Engine>
  class direct_callback : 
    public icallback<typename Engine::graph_type, 
                     typename Engine::update_functor_type> {
  public:
    typedef Engine engine_type;
    typedef typename engine_type::graph_type graph_type;
    typedef typename engine_type::update_functor_type update_functor_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef typename engine_type::ischeduler_type ischeduler_type;

  private:
    engine_type* engine_ptr;
    //    ischeduler_type* scheduler_ptr;
  public:
  
    direct_callback(engine_type* engine_ptr = NULL) :
      engine_ptr(engine_ptr) {} 

    void schedule(const vertex_id_type& vertex, 
                  const update_functor_type& update_fun) {
      engine_ptr->schedule(vertex, update_fun);
    }

    void schedule_in_neighbors(const vertex_id_type& vertex, 
                               const update_functor_type& update_fun) {
      logstream(LOG_FATAL) << "Unsupported call!" << std::endl;
    }

    virtual void schedule_in_neighbors(const vertex_id_type& vertex, 
                                       const update_functor_type& update_fun) {
      logstream(LOG_FATAL) << "Unsupported call!" << std::endl;
    }
    
    /**
     * Calling this function will force the engine to abort
     * immediately
     */
    virtual void terminate() {
      engine_ptr->stop();
    };

  };
}

#include <graphlab/macros_undef.hpp>


#endif

