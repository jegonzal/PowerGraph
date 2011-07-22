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



#ifndef GRAPHLAB_ICALLBACK_HPP
#define GRAPHLAB_ICALLBACK_HPP


#include <graphlab/macros_def.hpp>
namespace graphlab {




  /**
   * This is the class that is passed to update functors and allows
   * them to communicate with the engine.
   */
  template<typename Graph, typename UpdateFunctor>
  class icallback {
  public:
    typedef Graph graph_type;
    typedef UpdateFunctor update_functor_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    
    icallback() {}
    virtual ~icallback() {}

    /**
     * Adds a task to execute the update function on the vertex with
     * the given priority.
     */
    virtual void schedule(const vertex_id_type& vertex, 
                          const update_functor_type& update_fun) = 0;    


    /**
     * Schedule an update on all the neighbors of a particular vertex
     */
    virtual void schedule_in_neighbors(const vertex_id_type& vertex, 
                                       const update_functor_type& update_fun) = 0;

    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    virtual void schedule_in_neighbors(const vertex_id_type& vertex, 
                                       const update_functor_type& update_fun) = 0;
                                                  

    /**
     * Calling this function will force the engine to abort
     * immediately
     */
    virtual void terminate() = 0;
  };

}; //end graphlab namespace

#include <graphlab/macros_undef.hpp>
#endif

