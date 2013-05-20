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


#ifndef __DD_GRLAB_HPP__
#define __DD_GRLAB_HPP__

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>


#include <Eigen/Dense>
#include "eigen_serialization.hpp"

#include "dd_opts.hpp"
#include "utils.hpp"

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>


using namespace std;


/**
 * \brief Eigen library vectors are used to store the potentials (log-space)
 */
typedef Eigen::VectorXd factor_type;
typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;





double TOLERANCE = 0.01;      // The convergence threshold for each message. Smaller values imply tighter convergence but slower execution.
double old_dual = 100;        // stores the value of dual objective for the previous iteration
double primal_best  = 0;      //  stores the value of bestt primal objective found so far.
bool converged = false ;      // true if dual objective value has converged to required tolerance level, otherwise false
int dual_inc_count = 0;       // keeps track of the number of times the value of dual objective increased
int last_agg_count = 0;        // number of iterations since last aggregate_periodic call
mat history(30,3);
int history_iterator =0;

/////////////////////////////////////////////////////////////////////////
// Edge and Vertex data and Graph Type
/**
 * \brief There is a vertex for each factor in the graph AND each singleton
 */
struct vertex_data 
{
    int nvars;              // Number of variables in this factor.
    int degree;             // Degree of this factor (same as nvars for higher-order factors).
    
    vector<int> cards;      // Cardinality of each variable.
    vector<int> neighbors;  // Vertex ids of the neighbors.    
    vec potentials;         // Potentials for each configuration of the factor.
    
    int best_configuration; // Index of the best configuration at a subgradient step.
                            // TODO: Maybe replace best_configuration by beliefs for the high order variables?
                            // In which case, beliefs would be vector<vec> beliefs.

    double dual_contrib;    // Contribution of this factor to the dual. We can compute this during the applies. 
    double primal_contrib;  // Contribution of this factor to the primal. We can compute this during the applies. 
                            // NOTENOTE: ONLY true for sync engine. For async, we need to write an aggregator function. 
    
    vec beliefs;            // Posterior values for the configurations after averaging (projected DD, unary variables only).
    
    int apply_count;        // No. of times apply has been called on this vertex
    
    vertex_data(): 
    nvars(0), degree(0), 
    dual_contrib(0), primal_contrib(0),
    apply_count(0), best_configuration(0)
    {}
    
    void load(graphlab::iarchive& arc) 
    {
        arc >> nvars >> degree 
            >> cards >> neighbors >> potentials 
            >> dual_contrib >> primal_contrib
            >> best_configuration >> beliefs 
            >> apply_count;
    }
    void save(graphlab::oarchive& arc) const 
    {
        arc << nvars << degree 
            << cards << neighbors << potentials 
            << dual_contrib << primal_contrib
            << best_configuration << beliefs 
            << apply_count;
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
    
    vec potentials; // TODO: Unary potentials distributed evenly through the edges (i.e. unary potentials divided by degree).
    
    vec multiplier_messages; // Dual variables, i.e. Lagrangian multipliers.
    vec local_messages;      // Local MAP variables (for projected DD).
    
    edge_data(): varid(0), card(0) {}
    
    void load(graphlab::iarchive& arc) {
        arc >> varid >> card >> potentials >> multiplier_messages >> local_messages;
    }
    void save(graphlab::oarchive& arc) const {
        arc << varid << card << potentials << multiplier_messages << local_messages;
    }
};


/**
 * \brief gather_type is a structure that will be be used as return type of gather function. It includes 
 *  multiplier messages (used both for unary and factor vertices) and neighbor_best_conf (used only for 
 *  factor vertices).
 */

struct gather_type
{ factor_type messages;
  vector <int> neighbor_conf;
     
    gather_type(){};
    gather_type(factor_type f): messages(f){}; 
    gather_type(factor_type f, vector <int> nc): messages(f), neighbor_conf(nc){};
    void load(graphlab::iarchive& arc) {
        arc >>messages>>neighbor_conf;
    }
    void save(graphlab::oarchive& arc) const {
        arc <<messages<<neighbor_conf;
    }

  gather_type& operator+=(const gather_type& other)
 { messages += other.messages;
   neighbor_conf += other.neighbor_conf;
   return *this;
 }

}; // end of gather_type struct

/**
 * The graph type
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/** 
 * \brief The Dual Decomposition Vertex Program.
 */
struct dd_vertex_program : 
public graphlab::ivertex_program< graph_type, gather_type,
graphlab::messages::sum_priority >,
public graphlab::IS_POD_TYPE 
{
    
    /////////////////////////////////////////////////////////////////////////
    // Find the configuration index of a factor given the array of states.
    /////////////////////////////////////////////////////////////////////////
    int get_configuration_index(const graph_type::vertex_type& vertex,
                                const std::vector<int>& states) const 
    {
        const vertex_data& vdata = vertex.data();
        int index = states[0];
        for (size_t i = 1; i < states.size(); ++i) 
        {
            index *= vdata.cards[i];
            index += states[i];
        }
        return index;
    }
    
    /////////////////////////////////////////////////////////////////////////
    // Find the array of states corresponding to a factor configuration index.
    /////////////////////////////////////////////////////////////////////////
    void get_configuration_states(const graph_type::vertex_type& vertex,
                                  int index, std::vector<int>* states) const 
    {
        const vertex_data& vdata = vertex.data();
        int tmp = 1;
        for (size_t i = 1; i < states->size(); ++i) 
            tmp *= vdata.cards[i];
        
        (*states)[0] = index / tmp;
        for (size_t i = 1; i < states->size(); ++i) 
        {
            index = index % tmp;
            tmp /= vdata.cards[i];
            (*states)[i] = index / tmp;
        }
    }
    
    /**
     * \brief Given an edge and a vertex return the other vertex along
     * that edge. 
     */
    inline vertex_type get_other_vertex(edge_type& edge, 
                                        const vertex_type& vertex) const 
    {
        return vertex.id() == edge.source().id()? edge.target() : edge.source();
    }; // end of other_vertex
    
    
    virtual edge_dir_type gather_edges(icontext_type& context,
                                       const vertex_type& vertex) const = 0;
    virtual gather_type gather(icontext_type& context, const vertex_type& vertex, 
                               edge_type& edge) const = 0;
    virtual void apply(icontext_type& context, vertex_type& vertex, 
                       const gather_type& total) = 0;
    virtual edge_dir_type scatter_edges(icontext_type& context,
                                        const vertex_type& vertex) const = 0; 
    virtual void scatter(icontext_type& context, const vertex_type& vertex, 
                         edge_type& edge) const = 0;
}; // end of class bp_vertex_program



////////////////////////////////////////////////////////////////////////////////
// This class implements the "symmetric" version of dual decomposition described
// in:
// D. Sontag, A. Globerson, T. Jaakkola.
// Introduction to Dual Decomposition for Inference.
// Optimization for Machine Learning, editors S. Sra, S. Nowozin, and S. J.
// Wright: MIT Press, 2011
////////////////////////////////////////////////////////////////////////////////

struct dd_vertex_program_symmetric : public dd_vertex_program {
    /**
     * \brief Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type gather_edges(icontext_type& context,
                               const vertex_type& vertex) const 
    { 
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
    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const 
    {
        if (opts.verbose > 1)
            cout << "gather begin" << endl;
        
        const vertex_type other_vertex = get_other_vertex(edge, vertex);
        const vertex_data& vdata = vertex.data();
        edge_data& edata = edge.data();
        
        if (vdata.nvars == 1) 
        {
            // Unary factor.
            if (opts.verbose > 1)
                cout << "This unary factor has " << vertex.num_in_edges() << 
                " in edges and " << vertex.num_out_edges() << " out edges" << endl;

            if ((opts.verbose>1) && (vertex.id() == 0))
            {
                cout << "Gather on (" << vertex.id() << "," << other_vertex.id() << ") called from " << vertex.id() << "\n";
                cout << "vdata.neighbours = " << vdata.neighbors << "\n";
                cout << "Message: " << edata.multiplier_messages << "\n---\n";
            }
            
            gather_type gather_data(edata.multiplier_messages);
            return gather_data;
        } 
        else 
        {
            // General factor.
            factor_type messages;
            //messages.resize(vdata.potentials.size());
            messages.setZero(vdata.potentials.size());
            int offset = 0;
            int index_neighbor = -1;
            for (int k = 0; k < vdata.nvars; ++k) 
            {
                int vertex_id = vdata.neighbors[k];
                if (vertex_id == other_vertex.id()) 
                {
                    index_neighbor = k;
                    break;
                }
                offset += vdata.cards[k];
            }
            CHECK_GE(index_neighbor, 0);
            vector <int> neighbor_conf(vdata.nvars, 0);
            neighbor_conf[index_neighbor] = other_vertex.data().best_configuration;
            // cout<<vertex.id()<<" "<<other_vertex.id()<<" "<<neighbor_conf[index_neighbor]<<endl;
            for (int state = 0; state < vdata.cards[index_neighbor]; ++state) 
            {
                messages[offset + state] = -edata.multiplier_messages[state];
            }

            
            if ((opts.verbose>1) && (vertex.id() == 9))
            {
                cout << "Gather on (" << vertex.id() << "," << other_vertex.id() << ") called from " << vertex.id() << "\n";
                cout << "vdata.neighbours = " << vdata.neighbors << "\n";
                cout << "estimated offset = " << offset << "\n";
                cout << "Message: " << messages << "\n---\n";
            }
            gather_type gather_data(messages, neighbor_conf);
            return gather_data;
        }
        if (opts.verbose > 2)
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
    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) 
    {        
        vertex_data& vdata = vertex.data();
        
        if (opts.verbose > 1)
            cout << "begin apply" << endl;
        
        ++vdata.apply_count;
        ++last_agg_count;
        if (vdata.nvars == 1) 
        {
            // Unary factor.
            ASSERT_EQ(vdata.potentials.size(), total.messages.size());
            vec belief = vdata.potentials + total.messages;
            // Save the best configuration for this vertex.
            vdata.primal_contrib = vdata.potentials[vdata.best_configuration];
            //cout<<vertex.id()<<" "<<vdata.best_configuration<<endl;
            vdata.dual_contrib = belief.maxCoeff(&vdata.best_configuration);
            //cout<<"unary"<<" "<<vdata.dual_contrib<<"  "<<vdata.primal_contrib<<endl;
            if (opts.verbose > 1)
            {
                cout << "Vertex: " << vertex.id() << "\n";
                cout << "Potential: " << vdata.potentials << "\n";
                cout << "incomming message: " << total.messages << "\n";
                cout << "belief: " << belief << "\n";
                cout << "dual contrib: " << vdata.dual_contrib << "\n";                
                cout << "vdata.best_configuration = " << vdata.best_configuration << "\n---\n";
            }
        } 
        else 
        {
            // General factor.
            vec belief = vdata.potentials;
            int num_configurations = vdata.potentials.size();
            for (int index_configuration = 0;
                 index_configuration < num_configurations;
                 ++index_configuration) 
            {
                vector<int> states(vdata.nvars, -1);
                // This could be made more efficient by defining an iterator over factor
                // configurations.
                get_configuration_states(vertex, index_configuration, &states);
                int offset = 0;
                for (int k = 0; k < vdata.nvars; ++k) 
                {
                    belief[index_configuration] += total.messages[offset + states[k]];
                    offset += vdata.cards[k];
                }
            }
            // Save the best configuration for this factor.
            vdata.dual_contrib = belief.maxCoeff(&vdata.best_configuration);
            int conf_index = get_configuration_index(vertex, total.neighbor_conf);
            vdata.primal_contrib = vdata.potentials[conf_index];
            //cout<<vertex.id()<<" "<<belief<<endl;
            if (opts.verbose > 1)
            {
                cout << "Vertex: " << vertex.id() << "\n";
                cout << "Potential: " << vdata.potentials << "\n";
                cout << "incomming message: " << total.messages << "\n";
                cout << "belief: " << belief << "\n";
                cout << "dual contrib: " << vdata.dual_contrib << "\n";
                cout << "vdata.best_configuration = " << vdata.best_configuration << "\n---\n";
            }
        }

        if (opts.verbose > 1)
            cout << "end apply" << endl;
    }; // end of apply
    
    /**
     * \brief Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const 
    { 
        //return graphlab::ALL_EDGES; 
        // NOTENOTE: This assumes a sync engine. 
        return graphlab::OUT_EDGES; 
    }; // end of scatter edges
    
    /**
     * \brief The scatter function takes a vertex and an edge as input. 
     We just need to update the messages (Lagrange multipliers) by looking at the 
     saved argmaxes.
     */
    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const 
    {  
        const vertex_type other_vertex = get_other_vertex(edge, vertex);
        const vertex_type *unary_vertex;
        const vertex_type *factor_vertex;
                
        if (vertex.data().nvars == 1) 
        {
            // Unary factor.
            unary_vertex = &vertex;
            factor_vertex = &other_vertex;
        } 
        else 
        {
            // General factor.
            unary_vertex = &other_vertex;
            factor_vertex = &vertex;
        }
        const vertex_data& vdata = unary_vertex->data();
        const vertex_data& other_vdata = factor_vertex->data();
        edge_data& edata = edge.data();
        
        if (opts.verbose > 1)
            cout << "begin scatter" << endl;
          // double stepsize = 1;
        vector <double> stepsize(3,0);
        stepsize[0] = 1;
        stepsize[1] = 1.0/max(vdata.apply_count,other_vdata.apply_count);
        stepsize[2] = (1001*(old_dual-primal_best)) / ((last_agg_count)*(dual_inc_count+1000));
        CHECK_GE(vdata.best_configuration, 0);                                                            
        CHECK_LT(vdata.best_configuration, vdata.cards[0]);    
         //cout<< stepsize<<endl;
        // Negative subgradient
        edata.multiplier_messages[vdata.best_configuration] -= stepsize[opts.stepsize_type]; 
        
        vector<int> states(other_vdata.nvars, -1);
        get_configuration_states(*factor_vertex, other_vdata.best_configuration, &states);
        int offset = 0;
        int index_neighbor = -1;
        for (int k = 0; k < other_vdata.nvars; ++k) {
            int vertex_id = other_vdata.neighbors[k];
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

        // Negative subgradient
        edata.multiplier_messages[states[index_neighbor]] += stepsize[opts.stepsize_type];
        
        //if (opts.verbose > 1)
        if ((opts.verbose>1) && (vertex.id() == 15))
        {
            cout << "Scatter on (" << unary_vertex->id() << "," << factor_vertex->id() << ") called from " << vertex.id() << "\n";
            cout << "unary best config = " << vdata.best_configuration << "\n"
                 << "factor best config = " << states[index_neighbor] << "\n"; 
            cout << "Message: " << edata.multiplier_messages << "\n---\n";
        }
        
        if (opts.verbose > 1)
            cout << "end scatter" << endl;
        
        // Signalling the other vertex and yourself to start. 
        if (vertex.data().apply_count < opts.maxiter && converged == false)
        {
            context.signal(vertex);
            context.signal(other_vertex);
        }
    }; // end of scatter

}; // end of class bp_vertex_program_symmetric


/////////////////////////////////////////////////////////////////////////////////
//Aggregator functions to compute primal & dual objective

struct objective
{ double primal, dual;
 
objective(): primal(0), dual(0){};

void load(graphlab::iarchive& arc) {
        arc >>dual>>primal;
    }
    void save(graphlab::oarchive& arc) const {
        arc <<dual<<primal;
    }

objective& operator+=(const objective& other)
{ primal += other.primal;
   dual += other.dual;
   return *this;
 }
};

objective sum(dd_vertex_program_symmetric::icontext_type& context, const dd_vertex_program_symmetric::vertex_type& vertex)
{ objective retval;
   retval.primal = vertex.data().primal_contrib;
  retval.dual = vertex.data().dual_contrib;
  return retval;
}

void print_obj(dd_vertex_program_symmetric::icontext_type& context, objective total) 
{
    cout << "Dual Objective: " << total.dual<< " "<<"Primal Objective: "<<total.primal<<"\n";
    if (std::fabs(total.dual-old_dual) < TOLERANCE)
      { converged = true;
        cout<<" Number of iteration at convergence:"<<context.iteration() +2 <<endl;}
    else old_dual = total.dual;
     if (total.primal> primal_best)
      { primal_best = total.primal;}
      if (total.dual > old_dual)
      { dual_inc_count ++;}
    cout<< "Best Primal so far:"  << primal_best<<endl;
    last_agg_count = 0;

if(opts.save_history>0)
{history(history_iterator,0) = context.iteration() +2;
 history(history_iterator,1) = total.dual;
 history(history_iterator,2) = total.primal;
 history_iterator++;
}
}
/* end of aggregator functions */


struct conf_writer {
    std::string save_vertex(graph_type::vertex_type v) {
        std::stringstream strm;
         if(history_iterator != -1)
       { strm <<history<<endl;
       history_iterator = -1; }
        return strm.str();
    }
    std::string save_edge(graph_type::edge_type e) { return ""; }
};



////////////////////////////////////////////////////////////////////////////////
// This class implements the "projected" version of dual decomposition described
// in:
// Komodakis, N., Paragios, N., and Tziritas, G. (2007).
// MRF optimization via dual decomposition: Message-passing revisited.
// In Proc. of International Conference on Computer Vision.
// 
// The formulation used is the one in Algorithm 1 of:
//
// André F. T. Martins, Mário A. T. Figueiredo, Pedro M. Q. Aguiar,
// Noah A. Smith, and Eric P. Xing.
// "An Augmented Lagrangian Approach to Constrained MAP Inference."
// International Conference on Machine Learning (ICML), 2011.
////////////////////////////////////////////////////////////////////////////////

struct dd_vertex_program_projected : public dd_vertex_program {
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
     over all edges incident in this vertex. 
     If the vertex is a unary factor, compute the sum of all the local MAP variables,
     which in the "apply" function will serve to compute the global MAP. 
     Otherwise (if vertex is a general factor), things are a little more tricky. 
     Suppose the factor is linked to K variables, with cardinalities C_1, ..., C_K. 
     Suppose this edge is with respect to the k-th variable. Then, we return a 
     vector of size C_1 + ... + C_K which is zero everywhere except in the 
     k-th slot, where the Lagrange multipliers in "edge.messages" will be copied 
     to. This way, when the "gather sum" takes place, and since all these slots 
     are disjoint, we will just get the Lagrange multipliers of all the variables.    
     */
    gather_type gather(icontext_type& context, const vertex_type& vertex, 
                       edge_type& edge) const {
        cout << "gather begin" << endl;
        const vertex_type other_vertex = get_other_vertex(edge, vertex);
        const vertex_data& vdata = vertex.data();
        edge_data& edata = edge.data();
        
        if (vdata.nvars == 1) 
        {
            // Unary factor.
            cout << "This unary factor has " << vertex.num_in_edges() << 
            " in edges and " << vertex.num_out_edges() << " out edges" << endl;
            gather_type gatherdata(edata.local_messages);
            return gatherdata; 
        } 
        else 
        {
            // General factor.
            factor_type messages;
            //messages.resize(vdata.potentials.size());
            messages.setZero(vdata.potentials.size());
            int offset = 0;
            int index_neighbor = -1;
            for (int k = 0; k < vdata.nvars; ++k) {
                int vertex_id = vdata.neighbors[k];
                if (vertex_id == other_vertex.id()) {
                    index_neighbor = k;
                    break;
                }
                offset += vdata.cards[k];
            }
            CHECK_GE(index_neighbor, 0);
            //const vec &unary_potential = other_vertex.data().potential;
            //int degree = other_vertex.data().degree;
            
            for (int state = 0; state < vdata.cards[index_neighbor]; ++state) {
                messages[offset + state] = edata.multiplier_messages[state];
                // TODO: somehow set the "edge potential" to be the potential of the
                // unary variable divided by the number of factors in which that 
                // variable appears.
                messages[offset + state] += edata.potentials[state]; 
                //message[offset + state] += unary_potential[state] / static_cast<double>(degree);
            }
            gather_type gather_data(messages);
            return gather_data;
        }
        cout << "gather end" << endl;
    }; // end of gather function
    
    /**
     * \brief The apply function takes a vertex and a vector of numeric values 
     (a total) as input. 
     For a unary vertex, "total" will be the sum of local MAP vectors, and we 
     just need to divide by the vertex degree and save the result as global MAP.
     For higher-order factors, "total" will contain all the Lagrange multipliers 
     of the neighboring variables. So we need to loop through all possible factor 
     configurations, get the sequence of states of each configuration, fetch the 
     Lagrange multipliers for those states, and add them to the factor potential. 
     Then we compute the argmax and save result to local MAP for each variable 
     connected to the factor. 
     */
    void apply(icontext_type& context, vertex_type& vertex, 
               const gather_type& total) {
        vertex_data& vdata = vertex.data();
        cout << "begin apply" << endl;
        if (vdata.nvars == 1) {
            // Unary factor. Divide by vertex degree.
            vdata.beliefs = total.messages / static_cast<double>(vdata.degree);
            return;
        } else {
            // General factor.
            vec beliefs = vdata.potentials;
            int num_configurations = vdata.potentials.size();
            for (int index_configuration = 0;
                 index_configuration < num_configurations;
                 ++index_configuration) {
                vector<int> states(vdata.nvars, -1);
                // This could be made more efficient by defining an iterator over factor
                // configurations.
                get_configuration_states(vertex, index_configuration, &states);
                int offset = 0;
                for (int k = 0; k < vdata.nvars; ++k) {
                    beliefs[index_configuration] += total.messages[offset + states[k]];
                    offset += vdata.cards[k];
                }
            }
            // Save the best configuration for this factor.
            beliefs.maxCoeff(&vdata.best_configuration);
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
     (1) If the vertex is a unary factor, we update the messages (Lagrange multipliers)
     by subtracting the global MAP by the local MAP.
     (2) If the vertex is a higher order factor, this function will take the best
     configuration (obtained at the apply function) and save the local MAP 
     at the corresponding edge.
     */
    void scatter(icontext_type& context, const vertex_type& vertex, 
                 edge_type& edge) const {  
        const vertex_data& vdata = vertex.data();
        edge_data& edata = edge.data();
        cout << "begin scatter" << endl;
        if (vdata.nvars == 1) {
            // Unary factor. Update the messages (Lagrange multipliers).      
            double stepsize = 1.0; // TODO: Make this decay over iteration number.
            edata.multiplier_messages += (vdata.beliefs - edata.local_messages) * stepsize;
        } else {
            // General factor. Update the local MAPs.
            const vertex_type &unary_vertex = get_other_vertex(edge, vertex);
            vector<int> states(vdata.nvars, -1);
            get_configuration_states(vertex, vdata.best_configuration, &states);
            int offset = 0;
            int index_neighbor = -1;
            for (int k = 0; k < vdata.nvars; ++k) {
                int vertex_id = vdata.neighbors[k];
                if (vertex_id == unary_vertex.id()) {
                    index_neighbor = k;
                    break;
                }
                offset += vdata.cards[k];
            }
            CHECK_GE(index_neighbor, 0);
            CHECK_GE(states[index_neighbor], 0);
            CHECK_LT(states[index_neighbor], vdata.cards[index_neighbor]);
            CHECK_EQ(vdata.cards[index_neighbor], unary_vertex.data().cards[0]);
            edata.local_messages.setZero();
            edata.local_messages[states[index_neighbor]] += 1.0;
        }
        cout << "end scatter" << endl;
    }; // end of scatter
}; // end of class dd_vertex_program_projected


#endif
