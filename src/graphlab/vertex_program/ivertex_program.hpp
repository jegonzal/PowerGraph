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
   *   1) edge_data_type: the type of the edge data
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

    ////// WORKING HERE:
    ///// ---------------------------------------------------------------
   
    typedef icontext<graph_type, update_functor_type> icontext_type;
    typedef iglobal_context iglobal_context_type;

    typedef graphlab::edge_set           edge_set;
    typedef graphlab::consistency_model  consistency_model;
    

    virtual ~ivertex_program() { }

    /**
     * Gets the context range required by this update functor.  If not
     * implemented by the derived class then the default context range
     * is returned.
     */
    inline virtual consistency_model consistency() const {
      return DEFAULT_CONSISTENCY;
    }

    /**
     * When multiple update functors are scheduled to be run on the
     * same function they are added. The default behavior is to simply
     * ignore the later update functors.
     */
    inline virtual void operator+=(const update_functor_type& other) const { }

    /**
     * Get the priority of the update functor
     */
    inline virtual double priority() const { return double(0); }        

    /**
     * The main part of an update functor
     */
    inline virtual void operator()(icontext_type& context) { 
      logstream(LOG_FATAL) << "Operator() not implemented!" << std::endl;
    } 

    /**
     * Returns true if the factorized (gather, apply, scatter) version
     * of the update functor is to be used.
     */
    inline virtual bool is_factorizable() const { return false; }
    
    /**
     * Returns the set of edges to gather 
     */
    inline virtual edge_set gather_edges() const { return IN_EDGES; }

    /**
     * Returns true of the adjacent edge and vertex are modified
     * during the gather.
     */
    inline virtual consistency_model gather_consistency() const { 
      return DEFAULT_CONSISTENCY;
    }

    /**
     * Returns the set of edges to scatter
     */
    inline virtual edge_set scatter_edges() const { return OUT_EDGES; }

    /**
     * Returns true of the adjacent edge and vertex are modified
     * during the gather.
     */
    inline virtual consistency_model scatter_consistency() const { 
      return DEFAULT_CONSISTENCY;
    }

    
    /**
     * Init gather is called before gathering
     */
    inline virtual void init_gather(icontext_type& context) { };

    /**
     * Gather is called on all gather_edges() and may be called in
     * parallel.  The merge() operation is used to join update
     * functors.
     */
    inline virtual void gather(icontext_type& context, const edge_type& edge) { 
      logstream(LOG_FATAL) << "Gather not implemented!" << std::endl;
    };

    /**
     * Merges update functors during the gather process.
     */
    inline virtual void merge(const update_functor_type& other) {
      logstream(LOG_FATAL) << "Merge not implemented!" << std::endl;
    }

    /**
     * Apply is called within the vertex consistency model on the
     * center vertex after all gathers have completed.
     */
    inline virtual void apply(icontext_type& context) { 
      logstream(LOG_FATAL) << "Apply not implemented!" << std::endl;
    };
    
    
    /**
     * Scatter is invoked on all scatter_edges() after calling
     * init_scatter() and may be called in parallel.
     */
    inline virtual void scatter(icontext_type& context, const edge_type& edge) { 
      logstream(LOG_FATAL) << "Scatter not implemented!" << std::endl;
    }
  };  // end of ivertex_program
 
}; //end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
