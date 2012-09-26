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
 struct vertex_data 
{
    int nvars; // no. of vars in this factor
    vector<int> cards; // cardinality of each var
    vector<int> xmap; 
    
    vec potential;

    vertex_data(): nvars(0)
    {}
    
    void load(graphlab::iarchive& arc) 
    { arc >> cards >> xmap >> potential; }
    void save(graphlab::oarchive& arc) const 
    { arc << cards << xmap << potential; }
}; // end of vertex_data


/**
 * \brief There is an edge connecting each factor to each singleton
 * in its scope.
 */
struct edge_data 
{
    int varid;
    int card;
    vec message; // dual variables, lagrangian multiplier

    edge_data(): varid(0), card(0)
    {}
    
    void load(graphlab::iarchive& arc) 
    { arc >> varid >> card >> message; }
    void save(graphlab::oarchive& arc) const 
    { arc << varid << card << message; }
};


/**
 * The graph type
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/////////////////////////////////////////////////////////////////////////
// Function to perform MAP on each subproblem
void subproblem_map(graph_type::vertex_type vertex)
{
    vertex_data& vdata = vertex.data();
    
    vec rep_pot = vdata.potential; 
    // todo: pull in messages from neighbours (ie reparameterize)
    // for loop on neighbours
//    if (card == 1) // unary factor
//         rep_pot += edata.message;
//    else // general factor
//    {} // subtract a unary potential to a multi-dim factor
    
    // maximize over potential
    int maxid; 
    rep_pot.maxCoeff(&maxid);
    
    // todo: convert linear index into multi-dim index
    // & save in xmap
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
