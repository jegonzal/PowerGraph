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
 * This file contains an example of graphlab used for MAP inference 
 * in a discrete graphical model (pairwise MRF). The algorithm
 * implemented is the MPLP LP-Relaxation scheme of Globerson & Jaakkola. 
 *
 *  \author Dhruv Batra
 */


#include <vector>
#include <string>
#include <fstream>

#include <distributed_graphlab.hpp>
#include <graphlab/engine/distributed_synchronous_engine.hpp>


#include <graphlab/util/stl_util.hpp>
#include <graphlab/macros_def.hpp>


/**
 * Each GraphLab vertex is a (pairwise) factor from the MRF
 */
struct vertex_data {
  // original pairwise potential
  binary_factor thetaij; 
  // original unary potentials of adjacent nodes. This double storage
  // (each node pot stored twice). TODO: store pointers?
  unary_factor thetai, thetaj; 
  // dual variables being optimized (or messages)
  unary_factor del_fi, del_fj; 
  int i,j; // MRF node id. 
  double lpdg; // priority
  // constructor
  vertex_data(): xi(0), xj(0), lpdg(0) {} 
  void save(graphlab::oarchive& arc) const {
    arc << thetaij << thetai << thetaj << del_fi << del_fj
        << i << j << lpdg;
  }
  void load(graphlab::iarchive& arc) {
    arc >> thetaij >> thetai >> thetaj >> del_fi >> del_fj
        >> i >> j >> lpdg;
  }
}; // End of vertex data

std::ostream& operator<<(std::ostream& out, const vertex_data& vdata) {
  return out << "LPDG_scorev=" << vdata.lpdg;
}


/**
 * The data associated with a pair of factors in a pairwise MRF
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  // primal labelling; We assume pairwise factors, so intersection has
  // a single node
  int xi; 
}; // End of edge data

/**
 * The graph type
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


// GraphLab Update Function ===================================================>

/**
 * The type passed around during the gather phase
 */
struct gather_type {
  // The key quantity is sum_del_fi(xi), ie the sum of messages from
  // all factors f (containing variable i) to i. I find it a little
  // weird that unary_factor does not ask cardinality
  unary_factor sum_del_fi, sum_del_fj;
  gather_type& operator+=(const gather_type& other) {
    // \todo use sum or product?  Currently messages are stored in log
    // form so if most of the math is +/- then this could get
    // inefficient.
    sum_del_fi += other.sum_del_fi;
    sum_del_fj += other.sum_del_fj;
  }
  void save(graphlab::oarchive& arc) const {
    arc << sum_del_fi << sum_del_fj;
  }
  void load(graphlab::iarchive& arc) {
    arc >> sum_del_fi >> sum_del_fj;
  }
}; // end of gather type



/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
class mplp_vertex_program : public ivertex_program<vertex_data,
                                                   edge_data,
                                                   gather_type,
                                                   graphlab::messages::sum> {    
public:

  // void save(graphlab::oarchive& arc) const { /** save members */ }
  // void load(graphlab::iarchive& arc) { /** load members */ }


  /**
   * The init function is called once for each vertex before the
   * start of the GraphLab program.  If the vertex program does not
   * implement this function then the default implementation (NOP)
   * is used.
   */
  void init(icontext_type& context, vertex_type& vertex) { /** NOP */ }



  // I assume that MRF nodes are numbered 1-N. Now if i<k and i<j and
  // (i,j) (j,k) are two factors then I assume that GraphLab vertex
  // (i,j) will have an outgoing edge to (j,k) (edge are small to
  // large number).  This assumption let's me figure out which node i
  // or j to use in vertex data.  
  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 


  /**
   * Recv message is called by the engine to receive a message to this
   * vertex program.  The vertex program can use this to initialize
   * any state before entering the gather phase.  If the vertex
   * program does not implement this function then the default
   * implementation (NOP) is used.
   */
  void recv_message(icontext_type& context, const vertex_type& vertex, 
                    const message_type& msg) { /** NOP */ }

   
  // Run the gather operation over all in edges
   gather_type gather(icontext_type& context, const vertex_type& vertex, 
                      edge_type& edge) const {
    graphlab::unary_factor tmp;
    vertex_data& v1,v2;
    int i,j;
    // I assume I can tell if edge is an incoming edge or outgoing edge
    if (edge.IS_INCOMING()) // made up function
      {
        v1 = context.vertex_data(edge.target()); 
        v2 = context.vertex_data(edge.source());
      }
    else
      {
        v2 = context.vertex_data(edge.target()); 
        v1 = context.vertex_data(edge.source());
      }
        
    // Accumulate message
    if (v1.i == v2.i)
      sum_del_fi += v1.sum_del_fi;
    else if (v1.i == v2.j)
      sum_del_fi += v1.sum_del_fj;
    else if (v1.j == v2.i)
      sum_del_fj += v1.sum_del_fi;
    else if (v1.j == v2.j)
      sum_del_fj += v1.sum_del_fj;
    else 
      error(-1);
  } // end of gather
    
    // Merge two MPLP update accumulators after running gather
  void merge(const mplp_vertex_program& other) 
  { sum_del_fi += other.sum_del_fi; sum_del_fj += other.sum_del_fj; }
    

  // Update the center vertex
  void apply(icontext_type& context) 
  {
    vertex_data& vdata = context.vertex_data();
        
    graphlab::unary_factor new_del_fi, new_del_fj;
        
    // Perform MPLP update
    // Now not sure about how to handle cardinality of del_fi
    // So here is the pseudo-code (in matlab notation)
    for (int xi=0; xi!=del_fi.arity(); ++xi)
      new_del_fi(xi) = -vdata.thetai(xi)/2 - sum_del_fi(xi)/2 + max(vdata.thetaij(xi,:) + sum_del_fj);
    for (int xj=0; xj!=del_fj.arity(); ++xj)
      new_del_fj(xj) = -vdata.thetai(xj)/2 - sum_del_fi(xj)/2 + max(vdata.thetaij(:,xj) + sum_del_fi);
        
    vdata.del_fi = new_del_fi;
    vdata.del_fj = new_del_fj;
        
  } // end of apply

    // Reschedule neighbors 
  void scatter(icontext_type& context, const edge_type& edge) 
  { // Nothing yet. Will hold the LPDG scheduling scheme. 
  } // end of scatter
}

