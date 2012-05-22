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
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved. 
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



#ifndef GRAPHLAB_IVERTEX_PROGRAM_HPP
#define GRAPHLAB_IVERTEX_PROGRAM_HPP


#include <graphlab/context/icontext.hpp>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/graph.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  
  /**
   * A vertex program represents the primary user define computation
   * in graphlab.  A unique instance of the vertex program is run on
   * each vertex in the graph and can interact with neighboring vertex
   * programs through the gather and scatter functions as well as by
   * messaging.
   *
   *
   * All user define vertex programs must extend the ivertex_program
   * interface and implement the gather, apply, and scatter
   * functions. In addition, all vertex programs must provide the
   * following types:
   * 
   *   1) vertex_data_type: the type of the data stored on each vertex
   *   2) edge_data_type: the type of the edge data
   *   3) gather_type: the type used in the gather phase
   *   4) message_type: The type used for messaging
   *
   * Both the gather_type and message_type must be serializable (i.e.,
   * a primitive or implement load/save) and must support the
   * commutative associative += operation.
   *   
   * To enable the interface to be aware of the user defined gather
   * and message types the ivertex_program interface takes the
   * vertex_program as an argument.
   * 
   */
  template<typename VertexProgram> 
  class ivertex_program {    
  public:

    // User defined type members ==============================================
    /**
     * The vertex program interface takes as an argument the vertex
     * program.
     */
    typedef VertexProgram vertex_program_type;

    /**
     * The type of the vertex data which must be defined by the
     * vertex-program.
     */
    typedef typename vertex_program_type::vertex_data_type vertex_data_type;


    /**
     * The type of the edge data which must be defined by the
     * vertex-program.
     */
    typedef typename vertex_program_type::edge_data_type edge_data_type;

    /**
     * The gather type which must be provided by the vertex program.
     */
    typedef typename vertex_program_type::gather_type gather_type;

    /**
     * The message type which must be provided by the vertex_program
     */
    typedef typename vertex_program_type::message_type message_type;


    // Graph specific type members ============================================
    /**
     * The graph type associative with this vertex program.  The
     * vertex data type is the vertex_program_type and the
     * edge_data_type is the type specified by the vertex program as
     * the edge_data_type.
     */
    typedef graph<vertex_program_type, edge_data_type> graph_type;

    /**
     * The unique integer id used to reference vertices in the graph.
     * TODO: add doc link to definition in graph_type
     */
    typedef typename graph_type::vertex_id_type vertex_id_type;
    
    /**
     * The opaque vertex object type used to get vertex information.
     * TODO: add doc link to definition in graph_type
     */
    typedef typename graph_type::vertex_type vertex_type;
    
    /**
     * The opaque edge_object type used to access edge information.
     * TODO: add doc link to definition in graph_type
     */
    typedef typename graph_type::vertex_type edge_type;

    /**
     * The type used to define the direction of edges used in gather
     * and scatter.
     * TODO: add doc link to definition in graph_type
     */
    typedef graphlab::edge_dir_type edge_dir_type;

    // Additional Types =======================================================
    
    /**
     * The context type is used by the vertex program to communicate
     * with the engine and provides facilities for sending messages,
     * posting deltas, and accessing engine state.
     */
    typedef icontext<vertex_program_type> icontext_type;
   
    // Functions ==============================================================
    virtual ~ivertex_program() { }

    /**
     * The init function is called once for each vertex before the
     * start of the GraphLab program.  If the vertex program does not
     * implement this function then the default implementation (NOP)
     * is used.
     */
    virtual void init(icontext_type& context,
                      vertex_type& vertex) { /** NOP */ }

    /**
     * Recv message is called by the engine to receive a message to
     * this vertex program.  The vertex program can use this to
     * initialize any state before entering the gather phase.  If the
     * vertex program does not implement this function then the
     * default implementation (NOP) is used.
     */
    virtual void recv_message(icontext_type& context,
                              const vertex_type& vertex, 
                              const message_type& msg) { /** NOP */ }
    
    /**
     * Returns the set of edges on which to run the gather function.
     * The default edge direction is the in edges.
     */
    virtual edge_dir_type gather_edges(icontext_type& context,
                                       const vertex_type& vertex) const { 
      return IN_EDGES; 
    }

    /**
     * Gather is called on all gather_edges() in parallel and returns
     * the gather_type which are added to compute the final output of
     * the gather.
     */
    virtual gather_type gather(icontext_type& context, 
                               edge_type& edge) const {
      logstream(LOG_FATAL) << "Gather not implemented!" << std::endl;
    };

    /**
     * The apply function is called once the gather has completed and
     * must be implemented by all vertex programs. 
     */
    virtual void apply(icontext_type& context, 
                       vertex_type& vertex, 
                       const gather_type& total) = 0;

    /**
     * Returns the set of edges on which to run the scatter function.
     * The default edge direction is the out edges.
     */
    virtual edge_dir_type scatter_edges(icontext_type& context,
                                        const vertex_type& vertex) const { 
      return OUT_EDGES; 
    }

    /**
     * Scatter is called on all scatter_edges() in parallel after the
     * apply function has completed.  The scatter function can post
     * deltas.
     */
    virtual void scatter(icontext_type& context, edge_type& edge) const { 
      logstream(LOG_FATAL) << "Scatter not implemented!" << std::endl;
    };

  };  // end of ivertex_program
 
}; //end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
