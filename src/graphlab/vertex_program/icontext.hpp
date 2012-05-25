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


#ifndef GRAPHLAB_ICONTEXT_HPP
#define GRAPHLAB_ICONTEXT_HPP

#include <set>
#include <vector>
#include <cassert>

#include <graphlab/macros_def.hpp>
namespace graphlab {



  /**
   * \brief The context object mediates the interaction between the
   * vertex program and the graphlab execution environment.
   *
   * Each of the vertex program methods is passed a reference to the
   * engine's context.  
   */
  template<typename VertexType,
           typename GatherType, 
           typename MessageType>
  class icontext {
  public:
    // Type members ===========================================================

    /** 
     * The opaque vertex object type 
     * TODO: add a reference back to the graph type
     */
    typedef VertexType vertex_type;   

    /**
     * The message type specified by the user-defined vertex-program.
     * TODO: add a reference back to vertex program type
     */
    typedef MessageType message_type;

    /**
     * The type returned by the gather operation.
     * TODO: add a reference back to vertex program type
     */
    typedef GatherType gather_type;

   
  public:        
    /** icontext destructor */
    virtual ~icontext() { }
    
    /**
     * Get the number of vertices in the graph.
     */
    virtual size_t num_vertices() const = 0;

    /**
     * Get the number of edges in the graph
     */
    virtual size_t num_edges() const = 0;

    /**
     * Get an estimate of the number of update functions executed up
     * to this point.
     */
    virtual size_t num_updates() const = 0;

    /**
     * Get the elapsed time in seconds
     */
    virtual float elapsed_seconds() const = 0;

    /**
     * Return the current interation number (if supported).
     */
    virtual size_t iteration() const = 0;

    /**
     * Force the engine to stop executing additional update functions.
     */
    virtual void terminate() = 0;

    /**
     * Send a message to a vertex.
     */
    virtual void send_message(const vertex_type& vertex, 
                              const message_type& message = message_type()) = 0;

    /**
     * Post a change to the cached sum for the vertex
     */
    virtual void post_delta(const vertex_type& vertex, 
                            const gather_type& delta) = 0;    

    /**
     * Invalidate the cached gather on the vertex.
     */
    virtual void clear_gather_cache(const vertex_type& vertex) = 0; 

                                                

  }; // end of icontexty
  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

