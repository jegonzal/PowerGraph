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

#ifndef GRAPHLAB_IVERTEX_PROGRAM_HPP
#define GRAPHLAB_IVERTEX_PROGRAM_HPP


#include <graphlab/vertex_program/icontext.hpp>
#include <graphlab/util/empty.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/distributed_graph.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  
  /**
   * \brief The ivertex_program class defines the vertex program
   * interface that all vertex programs should extend and implement.
   *
   * Overview
   * ==================
   *
   * A vertex program represents the primary user define computation
   * in graphlab.  A unique instance of the vertex program is run on
   * each vertex in the graph and can interact with neighboring vertex
   * programs through the gather and scatter functions as well as by
   * messaging.  The vertex program's state is persistent throughout
   * the execution of the GraphLab program.
   *
   * Vertex programs express computation by implementing what we call
   * the *Gather-Apply-Scatter (GAS)* model which decomposes the
   * vertex program into a parallel gather phase, followed by an
   * atomic apply phase, and finally a parallel scatter phase.  This
   * decomposition allows us to execute a single vertex program on
   * several machines simultaneously and move computation to the data.
   *
   * \code
   * For vertex vtx:
   *   // During execution:
   *   if( there is a message for vtx ) {
   *     vprog.recv_message(ctx, vtx, msg);
   *     // Gather Phase: 
   *     vprog::gather_type sum;
   *     ParallelFor(adjacent edges in direction vprog.gather_edges(ctx, vtx) )
   *       sum += vprog.gather(ctx, vtx, edge);
   *     // Apply Phase
   *     vprog.apply(ctx, vtx, sum);
   *     // Scatter Phase
   *     ParallelFor(adjacent edges in direction vprog.scatter_edges(ctx, vtx) )
   *       vprog.scatter(ctx, vtx, edge);
   *     // Vertex program is destroyed
   *     vprog = vertex_program();
   *   }
   * \endcode
   *
   * All user define vertex programs must extend the ivertex_program
   * interface and implement the ivertex_program::apply function.
   * Most vertex programs will also implement the
   * ivertex_program::gather and ivertex_program::scatter functions as
   * well as the ivertex_program::init and
   * ivertex_program::recv_message functions.
   *
   * The state of a vertex program *does not* persist between
   * invocations of receive message.  Moreover on each call to
   * recv_message the vertex program's previous state is
   * cleared. Therefore any persistent state must be saved into the
   * vertex data.  
   *
   * The vertex program depends on several key types which are
   * template arguments to ivertex_program interface. 
   * 
   *   1) Graph: the type of graph used to store the data for this
   *      vertex program.  This is typically distributed_graph.
   *   2) gather_type: the type used in the gather phase
   *   3) message_type: The type used for messaging
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
  template<typename Graph,
           typename GatherType,
           typename MessageType = graphlab::empty> 
  class ivertex_program {    
  public:

    // User defined type members ==============================================
    /**
     * \brief The user defined vertex data type.  
     *
     * The vertex data is the data associated with each vertex in the
     * graph.  Unlike the vertex-program the vertex data of adjacent
     * vetices is visible to vertex programs during the gather and
     * scatter phases. In addition at termination the vertex-data is
     * accessible through the graph while the vertex-program state is
     * not.
     *
     * The vertex data type must be serializable (see \ref
     * serializable)
     */
    typedef typename Graph::vertex_data_type vertex_data_type;

    /**
     * The type of the edge data which must be defined by the
     * vertex-program.
     */
    typedef typename Graph::edge_data_type edge_data_type;

    /**
     * The gather type which must be provided by the vertex program.
     */
    typedef GatherType gather_type;
    
    /**
     * \cond GRAPHLAB_INTERNAL
     * \brief GraphLab Requires that the gather type be default
     * constructible.
     *
     * \code
     * class gather_type {
     * public:
     *   gather_type() { }
     * };  
     * \endcode
     */
    BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<GatherType>));
    /// \endcond

    /**
     * The message type which must be provided by the vertex_program
     */
    typedef MessageType message_type;


    /**
     * \cond GRAPHLAB_INTERNAL
     * GraphLab Requires that the message type be default
     * constructible.
     *
     * \code
     * class message_type {
     * public:
     *   message_type() { }
     * };  
     * \endcode
     * 
     */
    BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<MessageType>));
    /// \endcond 



    // Graph specific type members ============================================
    /**
     * The graph type associative with this vertex program.
     */
    typedef Graph graph_type;

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
    typedef typename graph_type::edge_type edge_type;

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
    typedef icontext<graph_type, gather_type, message_type> icontext_type;
   
    // Functions ==============================================================
    virtual ~ivertex_program() { }

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
                               const vertex_type& vertex, 
                               edge_type& edge) const {
      logstream(LOG_FATAL) << "Gather not implemented!" << std::endl;
      return gather_type();
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
    virtual void scatter(icontext_type& context, const vertex_type& vertex, 
                         edge_type& edge) const { 
      logstream(LOG_FATAL) << "Scatter not implemented!" << std::endl;
    };

  };  // end of ivertex_program
 
}; //end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
