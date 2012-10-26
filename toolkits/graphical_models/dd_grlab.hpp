/*  
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
 *
 * \brief This application performs MAP inference on Markov Nets 
 * provided in standard UAI file format via Dual-Decomposition. 
 *
 *
 *  \author Dhruv Batra
 */


#ifndef __DD_GRLAB_H__
#define __DD_GRLAB_H__

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>


#include <Eigen/Dense>
#include "eigen_serialization.hpp"



#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>


using namespace std;


/**
 * \brief Eigen library vectors are used to store the potentials (log-space)
 */
typedef Eigen::VectorXd factor_type;
typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;


/**
 * \brief The convergence threshold for each message.  Smaller values
 * imply tighter convergence but slower execution.
 *
 */
double TOLERANCE = 0.01;


/////////////////////////////////////////////////////////////////////////
// Edge and Vertex data and Graph Type
/**
 * \brief There is a vertex for each factor in the graph AND each singleton
 */
struct vertex_data {
    int nvars; // no. of vars in this factor.
    vector<int> cards; // cardinality of each variable.
    vector<int> xmap; // NOTE(afm): this was not used. So I defined this to be the
                      // neighbor vertex ids.    
    vec potential; // vector of potentials for each configuration of the factor.
    
    int best_configuration; // Index of the best configuration at a subgradient step.

    vertex_data(): nvars(0) {}
    
    void load(graphlab::iarchive& arc) {
      arc >> nvars >> cards >> xmap >> potential;
    }
    void save(graphlab::oarchive& arc) const {
      arc << nvars << cards << xmap << potential;
    }
}; // end of vertex_data


/**
 * \brief There is an edge connecting each factor to each singleton
 * in its scope.
 */
struct edge_data 
{
    int varid; // Do we need this? (afm)
    int card; // Do we need this? (afm)
    vec message; // dual variables, i.e. Lagrangian multipliers.
    //vector<double> message;

    edge_data(): varid(0), card(0) {}
    
    void load(graphlab::iarchive& arc) {
      arc >> varid >> card >> message;
    }
    void save(graphlab::oarchive& arc) const {
      arc << varid << card << message;
    }
};


/**
 * The graph type
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;




/**
 * \brief The graph type used to store the Markov Random Field with
 * vertex data containing node potentials and beliefs and edge data
 * containing messages and weights.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



/** 
 * \brief The Loopy Belief Propagation Vertex Program which computes
 * the product of the inbound messages during the gather phase,
 * updates the belief during the apply phase, and then computes the
 * new out-bound messages during the scatter phase.
 *
 * Since the gather phase is computing the product of the inbound
 * messages and the messages are stored in log form the resulting sum
 * operation is actually a vector sum and so the gather type is simply
 * the factor type and the operator+= operation for the factor type is
 * sufficient.
 *
 */
struct bp_vertex_program : 
  public graphlab::ivertex_program< graph_type, factor_type,
                                    graphlab::messages::sum_priority >,
  public graphlab::IS_POD_TYPE {
  
  /////////////////////////////////////////////////////////////////////////
  // Find the configuration index of a factor given the array of states.
  /////////////////////////////////////////////////////////////////////////
  int get_configuration_index(const graph_type::vertex_type& vertex,
                              const std::vector<int>& states) const {
    const vertex_data& vdata = vertex.data();
    int index = states[0];
    for (size_t i = 1; i < states.size(); ++i) {
      index *= vdata.cards[i];
      index += states[i];
    }
    return index;
  }

  /////////////////////////////////////////////////////////////////////////
  // Find the array of states corresponding to a factor configuration index.
  /////////////////////////////////////////////////////////////////////////
  void get_configuration_states(const graph_type::vertex_type& vertex,
                                int index, std::vector<int>* states) const {
    const vertex_data& vdata = vertex.data();
    int tmp = 1;
    for (size_t i = 1; i < states->size(); ++i) {
      tmp *= vdata.cards[i];
    }
    (*states)[0] = index / tmp;
    for (size_t i = 1; i < states->size(); ++i) {
      index = index % tmp;
      tmp /= vdata.cards[i];
      (*states)[i] = index / tmp;
    }
  }
  

  /**
   * \brief Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /**
   * \brief The gather function takes a vertex and an edge as inputs and outputs 
   a vector of numeric values. Vectors of numeric values will later be summed 
   over all edges incident in this vertex. So, if the vertex is a unary factor, 
   we can just return the vector of Lagrange multipliers stored in "edge.messages". 
   Otherwise (if vertex is a general factor), things are a little more tricky. 
   Suppose the factor is linked to K variables, with cardinalities C_1, ..., C_K. 
   Suppose this edge is with respect to the k-th variable. Then, we return a 
   vector of size C_1 + ... + C_K which is zero everywhere except in the 
   k-th slot, where the Lagrange multipliers in "edge.messages" will be copied 
   to. This way, when the "gather sum" takes place, and since all these slots 
   are disjoint, we will just get the Lagrange multipliers of all the variables.    
   */
  factor_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const {
    cout << "gather begin" << endl;
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    const vertex_data& vdata = vertex.data();
    edge_data& edata = edge.data();

    if (vdata.nvars == 1) {
      // Unary factor.
      cout << "This unary factor has " << vertex.num_in_edges() << 
        " in edges and " << vertex.num_out_edges() << " out edges" << endl;
      return edata.message;
    } else {
      // General factor.
      factor_type message;
      message.resize(vdata.potential.size());
      int offset = 0;
      int index_neighbor = -1;
      for (int k = 0; k < vdata.nvars; ++k) {
        int vertex_id = vdata.xmap[k];
        if (vertex_id == other_vertex.id()) {
          index_neighbor = k;
          break;
        }
        offset += vdata.cards[k];
      }
      CHECK_GE(index_neighbor, 0);
      for (int state = 0; state < vdata.cards[index_neighbor]; ++state) {
        message[offset + state] = -edata.message[state];
      }
      return message;
    }
    cout << "gather end" << endl;
  }; // end of gather function

  /**
   * \brief The apply function takes a vertex and a vector of numeric values 
   (a total) as input. For unary vertices, this will be the sum of Lagrange 
   multipliers, and we just need to sum that to the vertex potential and compute 
   the argmax. For general factors, the vector of numeric values, as stated above, 
   will contain all the Lagrange multipliers of the neighboring variables. 
   So we need to loop through all possible factor configurations, get the 
   sequence of states of each configuration, fetch the Lagrange multipliers for 
   those states, and add them to the factor potential. Then we compute the argmax. 
   */
  void apply(icontext_type& context, vertex_type& vertex, 
             const factor_type& total) {
    vertex_data& vdata = vertex.data();
    cout << "begin apply" << endl;
    if (vdata.nvars == 1) {
      // Unary factor.
      ASSERT_EQ(vdata.potential.size(), total.size());
      vec belief = vdata.potential + total;
      // Save the best configuration for this vertex.
      belief.maxCoeff(&vdata.best_configuration);
      cout << "vdata.best_configuration = " << vdata.best_configuration << endl;
    } else {
      // General factor.
      vec belief = vdata.potential;
      int num_configurations = vdata.potential.size();
      for (int index_configuration = 0;
           index_configuration < num_configurations;
           ++index_configuration) {
        vector<int> states(vdata.nvars, -1);
        // This could be made more efficient by defining an iterator over factor
        // configurations.
        get_configuration_states(vertex, index_configuration, &states);
        int offset = 0;
        for (int k = 0; k < vdata.nvars; ++k) {
          belief[index_configuration] += total[offset + states[k]];
          offset += vdata.cards[k];
        }
      }
      // Save the best configuration for this factor.
      belief.maxCoeff(&vdata.best_configuration);
    }
    cout << "end apply" << endl;
  }; // end of apply

  /**
   * \brief Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  /**
   * \brief The scatter function takes a vertex and an edge as input. 
   We just need to update the messages (Lagrange multipliers) by looking at the 
   saved argmaxes.
   */
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {  
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    const vertex_type *unary_vertex;
    const vertex_type *factor_vertex;
    cout << "begin scatter" << endl;
    if (vertex.data().nvars == 1) {
      // Unary factor.
      unary_vertex = &vertex;
      factor_vertex = &other_vertex;
    } else {
      // General factor.
      unary_vertex = &other_vertex;
      factor_vertex = &vertex;
    }
    const vertex_data& vdata = unary_vertex->data();
    const vertex_data& other_vdata = factor_vertex->data();
    edge_data& edata = edge.data();
    double stepsize = 1.0; // TODO: Make this decay over iteration number.

    CHECK_GE(vdata.best_configuration, 0);
    CHECK_LT(vdata.best_configuration, vdata.cards[0]);    
    edata.message[vdata.best_configuration] += stepsize;
    vector<int> states(other_vdata.nvars, -1);
    get_configuration_states(*factor_vertex, other_vdata.best_configuration, &states);
    int offset = 0;
    int index_neighbor = -1;
    for (int k = 0; k < other_vdata.nvars; ++k) {
      int vertex_id = other_vdata.xmap[k];
      if (vertex_id == unary_vertex->id()) {
        index_neighbor = k;
        break;
      }
      offset += other_vdata.cards[k];
    }
    CHECK_GE(index_neighbor, 0);
    CHECK_GE(states[index_neighbor], 0);
    CHECK_LT(states[index_neighbor], other_vdata.cards[index_neighbor]);
    CHECK_EQ(other_vdata.cards[index_neighbor], vdata.cards[0]);
    edata.message[states[index_neighbor]] -= stepsize;
    cout << "end scatter" << endl;
  }; // end of scatter
  
  /**
   * \brief Given an edge and a vertex return the other vertex along
   * that edge. 
   */
  inline vertex_type get_other_vertex(edge_type& edge, 
                                      const vertex_type& vertex) const {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
  }; // end of other_vertex
  
}; // end of class bp_vertex_program



#if 0
/////////////////////////////////////////////////////////////////////////
// Function to perform MAP on each subproblem
/////////////////////////////////////////////////////////////////////////
void subproblem_map(graph_type::vertex_type vertex)
{
    cout << "Called subproblem_map" << endl;
    vertex_data& vdata = vertex.data();
    
    vec rep_pot = vdata.potential; 
    // todo: pull in messages from neighbours (ie reparameterize)
    // for loop on neighbours
    //cout << "Factor with " << vdata.nvars << " variables." << endl;
    if (vdata.nvars == 1) // unary factor
    {
      cout << "This unary factor has " << vertex.num_in_edges() << 
        " in edges and " << vertex.num_out_edges() << " out edges" << endl;
      
      //edge_data& edata = 
      //  rep_pot += edata.message;
    }
    else // general factor
    {
        // subtract a unary potential to a multi-dim factor
        std::vector<int> states(vdata.nvars);
        for (int index = 0; index < rep_pot.size(); ++index)
        {
            get_configuration_states(vertex, index, &states);
            for (int i = 0; i < vdata.nvars; ++i)
            {
                // edge_data& edata = vertex.get_edge(i).data();
                int j = states[i];
                //rep_pot[index] -= edata.message[j]; 
            }
        }
    }
    
    // maximize over potential
    int maxid; 
    rep_pot.maxCoeff(&maxid);
    
    // todo: convert linear index into multi-dim index
    // & save in xmap
    get_configuration_states(vertex, maxid, &vdata.xmap);    
}


/////////////////////////////////////////////////////////////////////////
// Function to perform subgradient descent
void update_duals(graph_type::edge_type edge)
{
    // todo: get two vertices & use xmaps to update message
    // & save in xmap
    vertex_data v1 = edge.source().data(); 
    vertex_data v2 = edge.target().data(); 
    
    
}
#endif

#endif
