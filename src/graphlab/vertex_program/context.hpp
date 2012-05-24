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


#ifndef GRAPHLAB_CONTEXT_HPP
#define GRAPHLAB_CONTEXT_HPP

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
  template<typename Engine>
  class context : public icontext<typename Engine::vertex_program_type> {
  public:
    // Type members ===========================================================

    typedef Engine engine_type;

    /** the graph type used by this context */
    typedef typename engine_type::graph_type graph_type;
    typedef typename graph_type::lvid_type lvid_type;

    /** The type of the user-defined vertex program */
    typedef typename engine_type::vertex_program_type vertex_program_type;



    /** 
     * The opaque vertex object type 
     * TODO: add a reference back to the graph type
     */
    typedef typename vertex_program_type::vertex_type vertex_type;   

    /**
     * The message type specified by the user-defined vertex-program.
     * TODO: add a reference back to vertex program type
     */
    typedef typename vertex_program_type::message_type message_type;

    /**
     * The type returned by the gather operation.
     * TODO: add a reference back to vertex program type
     */
    typedef typename vertex_program_type::gather_type gather_type;


    typedef std::pair<vertex_type, message_type> message_pair_type;
    typedef std::vector<message_pair_type> message_buffer_type;

    typedef std::pair<vertex_type, gather_type> delta_pair_type;
    typedef std::vector<delta_pair_type> delta_buffer_type;


  private:
    engine_type& engine;
    graph_type& graph;
    

    message_buffer_type message_buffer;
    delta_buffer_type delta_buffer;
    

    
  public:        

    context(engine_type& engine, graph_type& graph) : 
      engine(engine), graph(graph) { }
    
    /**
     * Clear the message and delta buffers
     */
    void clear() {
      message_buffer.clear();
      delta_buffer.clear();
    } // end of clear

    message_buffer_type messages() { return message_buffer; }
    delta_buffer_type deltas() { return delta_buffer; }
    
    /**
     * Get the number of vertices in the graph.
     */
    size_t num_vertices() const { return graph.num_vertices(); }

    /**
     * Get the number of edges in the graph
     */
    size_t num_edges() const { return graph.num_edges(); }

    /**
     * Get an estimate of the number of update functions executed up
     * to this point.
     */
    size_t num_updates() const { return engine.num_updates(); }

    /**
     * Get the elapsed time in seconds
     */
    size_t elapsed_time() const { return engine.elapsed_time(); }

    /**
     * Return the current interation number (if supported).
     */
    size_t iteration() const { return engine.iteration(); }

    /**
     * Force the engine to stop executing additional update functions.
     */
    void terminate() { engine.stop(); }

    /**
     * Send a message to a vertex.
     */
    void send_message(const vertex_type& vertex, 
                      const message_type& message = message_type()) {
      message_buffer.push_back(message_pair_type(vertex, message));
    }

    /**
     * Post a change to the cached sum for the vertex
     */
    void post_delta(const vertex_type& vertex, 
                       const gather_type& delta) {
      delta_buffer.push_back(delta_pair_type(vertex, delta));
    }

    /**
     * Invalidate the cached gather on the vertex.
     */
    virtual void clear_gather(const vertex_type& vertex) { 
      
    }


                                                

  }; // end of context
  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

