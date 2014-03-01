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
 *  \authors Dhruv Batra, André Martins, Aroma Mahendru
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



struct dd_global_vars {

double old_dual ;        // stores the value of dual objective for the previous iteration
double primal_best ;     //  stores the value of best primal objective found so far.
bool converged ;         // true if dual objective value has converged to required tolerance level, otherwise false
int dual_inc_count ;     // keeps track of the number of times the value of dual objective increased
vector < vector<double> > history ; // stores dual and primal objective values
int sq_norm_g ;   //  stores the value of the square of the norm of the subgradient vector
int iter_at_aggregate ;  //  iteration number at the time of aggregate
graphlab::timer timer ; //  time object. Helps in finding the time elapsed.
dd_global_vars(): old_dual(200), primal_best(-1e10),
                  converged(false), dual_inc_count(1),
                  history(4,vector<double>()), 
                  sq_norm_g(100), 
                  iter_at_aggregate(0){};
} global_vars;

/* end of struct dd_global_vars */







/////////////////////////////////////////////////////////////////////////
// Edge and Vertex data and Graph Type
/**
 * \brief There is a vertex for each factor in the graph AND each singleton
 */
struct vertex_data 
{
    int nvars;                 // Number of variables in this factor.
    int factor_type;           //type of  factor : dense(0), budget(1)
    int degree;                // Degree of this factor (same as nvars for higher-order factors).
    
    vector<int> cards;         // Cardinality of each variable.
    vector<int> neighbors;     // Vertex ids of the neighbors.    
    vec potentials;            // Potentials for each configuration of the factor.
    int budget;                // Only for Budget factors
    vector <int> bound_states;  // Only for Budget Factors
    vector<int> unary_degree;   // Only for unary vertices
    
    int best_configuration;    // Index of the best configuration at a subgradient step.
                               // TODO: Maybe replace best_configuration by beliefs for the high order variables?
                               // In which case, beliefs would be vector<vec> beliefs.

    double dual_contrib;       // Contribution of this factor to the dual. We can compute this during the applies. 
    double primal_contrib;     // Contribution of this factor to the primal. We can compute this during the applies. 
                               // NOTENOTE: ONLY true for sync engine. For async, we need to write an aggregator function. 
    double dual_res_contrib;   // Contribution of this factor to the dual residual. (Used only for ADMM)                    
    double primal_res_contrib; // Contribution of this factor to the primal residual. (Used only for ADMM)   
    double primal_rel_contrib; // Contribution of this factor to the relaxed primal. (Used only for ADMM)   
    
    vec beliefs;               // Posterior values for the configurations after averaging (projected DD, unary variables only).
    vec factor_beliefs;        // Posterior value for factor variables
    int apply_count;           // No. of times apply has been called on this vertex
    int sum_sq_norm_g;         // Sum of square of norm of subgradient for each vertex (used only for factor vertices)
    
    bool schedule_vertex;      // Decides if vertex is to be scheduled for further iterations or not
    

    vertex_data(): 
    nvars(0), factor_type(0), degree(0),
    budget(0),
    best_configuration(0),
    dual_contrib(0), primal_contrib(0),
    dual_res_contrib(0), primal_res_contrib(0), primal_rel_contrib(0),
    sum_sq_norm_g(0),
    apply_count(0),
    schedule_vertex(true)
    {}
    
    void load(graphlab::iarchive& arc) 
    {
        arc >> nvars >> degree 
            >> cards >> neighbors >> potentials 
            >> dual_contrib >> primal_contrib
            >> best_configuration >> beliefs 
            >> apply_count>>factor_beliefs
            >>sum_sq_norm_g>>primal_res_contrib
            >>primal_rel_contrib>>dual_res_contrib
            >>schedule_vertex>>factor_type
            >>budget>>bound_states
            >>unary_degree;
    }
    void save(graphlab::oarchive& arc) const 
    {
        arc << nvars << degree 
            << cards << neighbors << potentials 
            << dual_contrib << primal_contrib
            << best_configuration << beliefs 
            << apply_count <<factor_beliefs
            << sum_sq_norm_g<<primal_res_contrib
            <<primal_rel_contrib<<dual_res_contrib
            <<schedule_vertex<<factor_type
            <<budget<<bound_states
            <<unary_degree;
    }
}; // end of vertex_data


/**
 * \brief There is an edge connecting each factor to each singleton
 * in its scope.
 */
struct edge_data 
{ 
    vec potentials; 
    
    vec multiplier_messages; // Dual variables, i.e. Lagrangian multipliers.
    vec local_messages;      // Local MAP variables (for projected DD).
    
    void load(graphlab::iarchive& arc) {
        arc >> potentials >> multiplier_messages >> local_messages;
    }
    void save(graphlab::oarchive& arc) const {
        arc << potentials << multiplier_messages << local_messages;
    }
};  //end of edge_data


/**
 * \brief gather_type is a structure that will be used as the return type of gather function. It includes 
 *  messages (used both for unary and factor vertices), neighbor_best_conf, neighbor_distribution (used only for 
 *  factor vertices) and sq_norm_g (for storing square of norm of subgradient) for each edge.
 */

struct gather_type
{ factor_type messages;
  factor_type multipliers;
  vector <int> neighbor_conf;
  vec neighbor_distribution;
  int sq_norm_g;
     
    gather_type():sq_norm_g(0){};
    
    
    gather_type(factor_type f, vector <int> nc = vector <int> (), int sg = 0, 
                             vec nd = vec() ): messages(f), neighbor_conf(nc),
                                   neighbor_distribution(nd), sq_norm_g(sg){};
    void load(graphlab::iarchive& arc) {
        arc >>messages>>neighbor_conf
            >>sq_norm_g>>neighbor_distribution
            >>multipliers;
    }
    void save(graphlab::oarchive& arc) const {
        arc <<messages<<neighbor_conf
            <<sq_norm_g<<neighbor_distribution
            <<multipliers;
    }

  gather_type& operator+=(const gather_type& other)
 { messages += other.messages;
   neighbor_conf += other.neighbor_conf;
   neighbor_distribution += other.neighbor_distribution;
   sq_norm_g += other.sq_norm_g;
   multipliers += multipliers;
   
   return *this;
 }

}; // end of gather_type struct


/**
 * \brief objective is a structure that is used as the summable data type for computing 
 * dual, primal objectives, residuals etc with aggregator map and reduce functions.
 */

struct objective
{ double primal, dual, primal_rel, sum_sq_norm_g, total_confs, dual_res, primal_res;
 
objective(): primal(0), dual(0), primal_rel(0), sum_sq_norm_g(0), total_confs(1),dual_res(1), primal_res(1){};

void load(graphlab::iarchive& arc) {
        arc >>dual>>primal>>sum_sq_norm_g
            >>primal_rel>>total_confs
            >>dual_res>>primal_res;
    }
    void save(graphlab::oarchive& arc) const {
        arc <<dual<<primal<<sum_sq_norm_g
            <<primal_rel<<total_confs
            <<dual_res<<primal_res;
    }

objective& operator+=(const objective& other)
{ primal += other.primal;
   dual += other.dual;
   sum_sq_norm_g += other.sum_sq_norm_g;
   primal_rel += other.primal_rel;
   total_confs += other.total_confs;
   dual_res += other.dual_res;
   primal_res += other.primal_res;
   return *this;
 }
}; // end of structure objective

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

    ///////////////////////////////////////////////////////////
    // Updates stepsize according to different stepsize rules. 
    ///////////////////////////////////////////////////////////

    double update_stepsize(icontext_type& context, int type, double old_dual, double primal_best,int norm_g_sq,
                                                                             int dual_inc_count,int iter_since_aggregate) const
   {  switch (type) {
      case 0: return opts.step_size;
              break;
      case 1: return(opts.step_size/(context.iteration()+2));
              break;
      case 2: return(2* opts.step_size *(old_dual-primal_best)/((norm_g_sq+1) * (iter_since_aggregate + dual_inc_count + 1)));
              break;
      case 3: return(opts.step_size/dual_inc_count + 1);
              break;
      case 4: if(context.iteration()+2 < 300)
              return( opts.step_size/(context.iteration()+2));
              else return(opts.step_size/(300));
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
    {   if(!opts.debug){
        return graphlab::ALL_EDGES; }
        else 
        return graphlab::NO_EDGES;
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
     It also gathers the best_configuration of neighbors for the factors and norm
     of subgradient value for each edge.    
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

            if (opts.verbose>1)
            {
                cout << "Gather on (" << vertex.id() << "," << other_vertex.id() 
                                     << ") called from " << vertex.id() << "\n";
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
            

            for (int state = 0; state < vdata.cards[index_neighbor]; ++state) 
            {
                messages[offset + state] = -edata.multiplier_messages[state];
            }

     
            if (opts.verbose>1) 
            {
                cout << "Gather on (" << vertex.id() << "," << other_vertex.id() << ") called from " << vertex.id() << "\n";
                cout << "vdata.neighbours = " << vdata.neighbors << "\n";
                cout << "estimated offset = " << offset << "\n";
                cout << "Message: " << messages << "\n---\n";
            }
            vector<int> states(vdata.nvars, -1);
            get_configuration_states(vertex, vdata.best_configuration, &states);
            int sq_norm_g = (states[index_neighbor] == other_vertex.data().best_configuration)?0:2;
            gather_type gather_data(messages, neighbor_conf);
            gather_data.sq_norm_g = sq_norm_g;
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
     It also computes dual and primal contribution for finding dual and primal
     objective values.
     */
    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) 
    {   if (!opts.debug){     
        vertex_data& vdata = vertex.data();
                
        if (vdata.nvars == 1) 
        {
            // Unary factor.
            ASSERT_EQ(vdata.potentials.size(), total.messages.size());
            
            vec belief = vdata.potentials + total.messages;
             // Find primal contrib
            vdata.primal_contrib = vdata.potentials[vdata.best_configuration];
             // Save the best configuration for this vertex and find dual contrib
            vdata.dual_contrib = belief.maxCoeff(&vdata.best_configuration);
            
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
            // Save the best configuration for this factor and find dual contrib
            vdata.dual_contrib = belief.maxCoeff(&vdata.best_configuration);
            //Find primal contrib
            int conf_index = get_configuration_index(vertex, total.neighbor_conf);
            vdata.primal_contrib = vdata.potentials[conf_index];
            //Find contribution fir sum of square of gradient
            vdata.sum_sq_norm_g = total.sq_norm_g;
            
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
       }
      else usleep(1);

      if (opts.verbose > 1)
            cout << "end apply" << endl;
    }; // end of apply
    
    /**
     * \brief Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const 
    { 
         return graphlab::ALL_EDGES; 
        // NOTENOTE: This assumes a sync engine. 
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
     if (!opts.debug){       
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
                        
        int iter_since_aggregate = (context.iteration()+2) - global_vars.iter_at_aggregate ;
        double stepsize = update_stepsize(context, 1, global_vars.old_dual, global_vars.primal_best, 
                                      global_vars.sq_norm_g, global_vars.dual_inc_count, iter_since_aggregate);
         
        CHECK_GE(vdata.best_configuration, 0);                                                            
        CHECK_LT(vdata.best_configuration, vdata.cards[0]);    
        
        // Negative subgradient
       
        edata.multiplier_messages[vdata.best_configuration] -= stepsize; 
        
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
        edata.multiplier_messages[states[index_neighbor]] += stepsize;
        
        //if (opts.verbose > 1)
        if (opts.verbose>1) 
        {
            cout << "Scatter on (" << unary_vertex->id() << "," << factor_vertex->id() << ") called from " << vertex.id() << "\n";
            cout << "unary best config = " << vdata.best_configuration << "\n"
                 << "factor best config = " << states[index_neighbor] << "\n"; 
            cout << "Message: " << edata.multiplier_messages << "\n---\n";
        }
       }
        if (opts.verbose > 1)
            cout << "end scatter" << endl;
        
        // Signalling the other vertex and yourself to start. 
        if ((context.iteration()+1) < opts.maxiter && global_vars.converged == false)
        {
            context.signal(vertex);
            context.signal(other_vertex);
        }

    }; // end of scatter

}; // end of class dd_vertex_program_symmetric



////////////////////////////////////////////////////////////////////////////////
// This class implements the "projected" version of dual decomposition described
// in:
// Komodakis, N., Paragios, N., and Tziritas, G. (2007).
// "MRF optimization via dual decomposition: Message-passing revisited"
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

     graphlab::timer vertex_timer;
    /**
     * \brief Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type gather_edges(icontext_type& context,
                               const vertex_type& vertex) const { 
        if(!opts.debug){
        return graphlab::ALL_EDGES; }
        else 
        return graphlab::NO_EDGES;
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
        
        const vertex_type other_vertex = get_other_vertex(edge, vertex);
        const vertex_data& vdata = vertex.data();
        edge_data& edata = edge.data();
        if (vdata.nvars == 1 ) 
        {
            // Unary factor.
            if (opts.verbose > 1) {
            cout << "This unary factor has " << vertex.num_in_edges() << 
            " in edges and " << vertex.num_out_edges() << " out edges" << endl; }
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
            vector <int> neighbor_conf(vdata.nvars, 0);
            neighbor_conf[index_neighbor] = other_vertex.data().best_configuration;
            
            for (int state = 0; state < vdata.cards[index_neighbor]; ++state) {
                messages[offset + state] = edata.multiplier_messages[state] + edata.potentials[state];               
            }
            gather_type gather_data(messages,neighbor_conf);
            return gather_data;
        }
       // cout << "gather end" << endl;
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
     Note that since global MAP is computed from local MAP , locla MAP is needed 
     to be computed before global MAP. Hence in even iterations local MAP is 
     computed ( in scatter step) and in the subsequent iteration (which is hence 
     odd) global MAP (in apply step) and multiplier messages (in scatter)
     are updated.
     */
    void apply(icontext_type& context, vertex_type& vertex, 
               const gather_type& total) {
        vertex_data& vdata = vertex.data();
        //cout << "begin apply" << endl;
     if (!opts.debug){
        if (vdata.nvars == 1 ) {
           if (context.iteration()%2 != 0) 
        {   
                                               
            vdata.beliefs = total.messages / static_cast<double>(vdata.degree);
            vdata.beliefs.maxCoeff(&vdata.best_configuration);
            //Find primal contrib
            vdata.primal_contrib = vdata.potentials[vdata.best_configuration];
                   
            }
        } 
        else 
         {if(context.iteration()%2 == 0)
           {
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
            // Save the best configuration for this factor and find dual contrib             
            vdata.dual_contrib = beliefs.maxCoeff(&vdata.best_configuration);
            //Find primal contrib
            int conf_index = get_configuration_index(vertex, total.neighbor_conf);
            vdata.primal_contrib = vdata.potentials[conf_index]; 
                    
            }
        }
      }
      else usleep(1);
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
        const vertex_type other_vertex = get_other_vertex(edge, vertex);
        if(!opts.debug){
        if (vdata.nvars == 1 ) {
            if (context.iteration()%2 != 0) {
            // Unary factor. Update the messages (Lagrange multipliers).      
             double stepsize = update_stepsize(context, 1, global_vars.old_dual, global_vars.primal_best, 
                                              global_vars.sq_norm_g, global_vars.dual_inc_count, 0);
            
             edata.multiplier_messages += (vdata.beliefs - edata.local_messages) * stepsize;
          }  
        } 
       else 
       {   if (context.iteration()%2 == 0) {
            //General factor. Update the local MAPs.
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
            edata.local_messages[states[index_neighbor]] += 1.0; }
           }
        }
        if ((context.iteration()+1) < opts.maxiter && global_vars.converged == false)
        {
            context.signal(vertex);
            context.signal(other_vertex);
        }
    }; // end of scatter
}; // end of class dd_vertex_program_projected





///////////////////////////////////////////////////////////////////////////////////
// This class implements the general Alternating Directions Method of Multipliers.
//  
//  The formulation used is the one in Algorithm 2 of:
//
// André F. T. Martins, Mário A. T. Figueiredo, Pedro M. Q. Aguiar,
// Noah A. Smith, and Eric P. Xing.
// "Alternating Directions Dual Decomposition"
// Arxiv preprint arXiv:1212.6550, 2012.

///////////////////////////////////////////////////////////////////////////////////
 
struct admm_vertex_program:public dd_vertex_program {
  
  typedef int Configuration;
  
  // Function to solve each quadratic programming sub problem 
  virtual void compute_beliefs(vertex_type& vertex,const gather_type& total,
                 vec& variable_posteriors, vec& additional_posteriors) = 0;
  virtual void SolveMAP(vertex_type& vertex,const gather_type& total,
                 vec& variable_posteriors, vec& additional_posteriors, double& value) = 0;
                                 
                                 
   /**
     * \brief Since the MRF is undirected we will use all edges for gather and
     * scatter
     */
    edge_dir_type gather_edges(icontext_type& context,
                               const vertex_type& vertex) const { 
        if(!opts.debug){
        return graphlab::ALL_EDGES; }
        else 
        return graphlab::NO_EDGES;
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
        //cout << "gather begin" << endl;
        const vertex_type other_vertex = get_other_vertex(edge, vertex);
        const vertex_data& vdata = vertex.data();
        edge_data& edata = edge.data();
        if (vdata.nvars == 1 ) {  
            // Unary factor.
            if (opts.verbose > 1){
            cout << "This unary factor has " << vertex.num_in_edges() << 
            " in edges and " << vertex.num_out_edges() << " out edges" << endl; 
            }
            gather_type gatherdata(edata.local_messages);
            return gatherdata; 
        } 
        else if(vdata.factor_type ==0){  
            // Dense factor.
            factor_type messages, neighbor_distribution, multipliers;        
            messages.setZero(vdata.potentials.size());
            neighbor_distribution.setZero(vdata.potentials.size());
            multipliers.setZero(vdata.potentials.size());
            int offset = 0;
            int index_neighbor = -1;
            for (int k = 0; k < vdata.nvars; ++k) {
                int vertex_id = vdata.neighbors[k];
                if (vertex_id == other_vertex.id()) {
                    index_neighbor = k;
                    break;}
                offset += vdata.cards[k];
            }
            CHECK_GE(index_neighbor, 0);
            
            vector <int> neighbor_conf(vdata.nvars, 0);
            neighbor_conf[index_neighbor] = other_vertex.data().best_configuration;
           

            for (int state = 0; state < vdata.cards[index_neighbor]; ++state) {
                messages[offset + state] = edata.multiplier_messages[state] + edata.potentials[state]; 
                multipliers[offset + state] = edata.multiplier_messages[state]; 
                neighbor_distribution[offset + state] = other_vertex.data().beliefs[state];
            }
            
            gather_type gather_data(messages,neighbor_conf);
            gather_data.neighbor_distribution = neighbor_distribution;
            gather_data.multipliers = multipliers;
            return gather_data;
        }
            
        else if(vdata.factor_type ==1) {
            //Budget factor
            factor_type messages, neighbor_distribution;        
            messages.setZero(vdata.nvars);
            neighbor_distribution.setZero(vdata.nvars);
            int index_neighbor = -1;
            for (int k = 0; k < vdata.nvars; ++k) {
                int vertex_id = vdata.neighbors[k];
                if (vertex_id == other_vertex.id()) {
                    index_neighbor = k;
                    break;}
            }
           messages[index_neighbor] =  edata.multiplier_messages[0] + edata.potentials[vdata.bound_states[index_neighbor]];

           neighbor_distribution[index_neighbor] = other_vertex.data().beliefs[vdata.bound_states[index_neighbor]];
           
           gather_type gather_data(messages);
           gather_data.neighbor_distribution = neighbor_distribution;
           return gather_data;
          
        }
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
     Then we compute solution of Quadratic subproblem and save result to local 
     MAP for each variable connected to the factor. 
     Note that since global MAP is computed from local MAP , locla MAP is needed 
     to be computed before global MAP. Hence in even iterations local MAP is 
     computed ( in scatter step) and in the subsequent iteration (which is hence 
     odd) global MAP (in apply step) and multiplier messages (in scatter)
     are updated.
     */  
  
  void apply(icontext_type& context, vertex_type& vertex, 
               const gather_type& total) {
        vertex_data& vdata = vertex.data();
     if (!opts.debug){
        if (vdata.nvars == 1 ) {   
             if (context.iteration()%2 != 0) {
            // Unary factor. 
            //Find dual residual contrib
                vec dual_res_contrib;
                dual_res_contrib.setZero(vdata.cards[0]);
                for(int i=0; i<vdata.cards[0]; i++){
                   dual_res_contrib[i] =  (total.messages[i] / static_cast<double>(vdata.unary_degree[i]))
                                                                                - vdata.beliefs[i] ;
                   dual_res_contrib[i] = pow(dual_res_contrib[i],2); 
            // update global MAP     
                   vdata.beliefs[i] = total.messages[i] / static_cast<double>(vdata.unary_degree[i]);
                }
                vdata.dual_res_contrib = dual_res_contrib.sum(); 
            // Find best configuration
                vdata.beliefs.maxCoeff(&vdata.best_configuration);
            //Find relaxed primal contribution
                vdata.primal_rel_contrib = vdata.potentials.dot(vdata.beliefs);
            // Find primal contribution
                vdata.primal_contrib = vdata.potentials[vdata.best_configuration];
            }
        } 
        else{  
            if(context.iteration()%2 == 0){
            // Dense and Budget factors
               vec additional_posteriors, variable_posteriors;
               additional_posteriors.setZero(vdata.potentials.size());
               variable_posteriors.setZero(vdata.potentials.size());
               if(vdata.factor_type == 1){
                 additional_posteriors.setZero(vdata.nvars);
                 variable_posteriors.setZero(vdata.nvars); 
                }
            //Find dual contrib
                SolveMAP(vertex, total, variable_posteriors, additional_posteriors, vdata.dual_contrib);
            // Find relaxed primal contribution
                if(vdata.factor_type == 0) {
                  vdata.primal_rel_contrib = vdata.potentials.dot(vdata.factor_beliefs);
                }
            //Find primal residual contribution
                vec primal_res_contrib = vdata.beliefs - total.neighbor_distribution;
                for(int i=0;i< vdata.beliefs.size(); i++){
                   primal_res_contrib[i] = pow(primal_res_contrib[i],2.0);
                }
                vdata.primal_res_contrib = primal_res_contrib.sum();
            // Compute QP subproblem solution
                compute_beliefs(vertex, total, vdata.beliefs, vdata.factor_beliefs); 
            //Find primal contrib
                int conf_index = get_configuration_index(vertex, total.neighbor_conf);
                vdata.primal_contrib = vdata.potentials[conf_index]; 
            } 
          }
       }
      else usleep(1);
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
        const vertex_type other_vertex = get_other_vertex(edge, vertex);
        if (!opts.debug){
         // Unary factor. Update the messages (Lagrange multipliers).   
            if (vdata.nvars == 1 ){ 
                if (context.iteration()%2 != 0) {
           
                    double stepsize = update_stepsize(context, 0 , global_vars.old_dual, global_vars.primal_best, 
                                                       global_vars.sq_norm_g, global_vars.dual_inc_count, 0);
                    if(other_vertex.data().factor_type == 1){
                       int index_neighbor = -1;
                       for(int i=0; i < other_vertex.data().nvars; i++) {
                           if(other_vertex.data().neighbors[i] == vertex.id()){
                              index_neighbor =i;
                              break; 
                           }
                       }
                      edata.multiplier_messages[0] +=(vdata.beliefs[other_vertex.data().bound_states[index_neighbor]]
                         - edata.local_messages[other_vertex.data().bound_states[index_neighbor]]) * stepsize;
                     }
                     else if(other_vertex.data().factor_type == 0){
                            edata.multiplier_messages += (vdata.beliefs - edata.local_messages) * stepsize; 
                     } 
                }    
             } 
             else if(vdata.factor_type == 0){   
                     if (context.iteration()%2 == 0) {
            //General factor. Update the local MAPs.
                        const vertex_type &unary_vertex = get_other_vertex(edge, vertex);
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
                        CHECK_EQ(vdata.cards[index_neighbor], unary_vertex.data().cards[0]);
            
                        for(int state = 0; state < vdata.cards[index_neighbor]; state++){
                           edata.local_messages[state] = vdata.beliefs[offset+state];
                        }
                   } 
              }
              else if(vdata.factor_type == 1){ 
               // Budget factor. Update local MAPs
                      if (context.iteration()%2 == 0) {
                         const vertex_type &unary_vertex = get_other_vertex(edge, vertex);
                         edata.local_messages.setZero();
                         int index_neighbor = -1;
                         for (int k = 0; k < vdata.nvars; ++k) {
                             int vertex_id = vdata.neighbors[k];
                             if (vertex_id == unary_vertex.id()) {
                                index_neighbor = k;
                                break;
                             }
                         }
                         edata.local_messages[vdata.bound_states[index_neighbor]] = 
                                                    vdata.beliefs[index_neighbor];
                      }
               } 
        }
        if ((context.iteration()+1) < opts.maxiter && global_vars.converged == false) {
            context.signal(vertex);
            context.signal(other_vertex);
        }
    }; // end of scatter
    
 }; /* end of admm_vertex_program */
 
 


/////////////////////////////////////////////////////////////////////////////////
//Aggregator functions to compute primal & dual objectives and residuals

objective sum(dd_vertex_program::icontext_type& context, const dd_vertex_program::vertex_type& vertex){
  objective retval;
  retval.primal = vertex.data().primal_contrib;
  retval.dual = vertex.data().dual_contrib;
  retval.sum_sq_norm_g = vertex.data().sum_sq_norm_g;
  retval.primal_rel = vertex.data().primal_rel_contrib;
  retval.dual_res  = vertex.data().dual_res_contrib;
  retval.primal_res = vertex.data().primal_res_contrib;
  retval.total_confs =(vertex.data().nvars ==1)?(std::accumulate(vertex.data().unary_degree.begin(),
                                                  vertex.data().unary_degree.end(),0)):0;
  global_vars.iter_at_aggregate = (context.iteration() +2); 
  
  return retval;
}

void print_obj(dd_vertex_program::icontext_type& context, objective total) {   
      if (context.iteration() % 2 == 0 || opts.algorithm == 0) { 
         if (total.dual > global_vars.old_dual){
             global_vars.dual_inc_count ++;
         }
     
         if (total.primal> global_vars.primal_best) {
            global_vars.primal_best = total.primal;
         }

         if(opts.verbose >0) {
            cout<<"iteration: "<<context.iteration()<<" Dual Objective: " << total.dual<<
                   " "<<"Primal Objective: "<<total.primal<<endl;
            cout<< "Best Primal so far: "  << global_vars.primal_best<<" "<<endl;   
      
            if(opts.algorithm == 2) { 
              cout<<"Relaxed Primal Objective:"<<total.primal_rel<<endl;
              cout<<"Dual Residual:"<< sqrt(total.dual_res/total.total_confs)<<" "
              <<"Primal Residual:"<<sqrt(total.primal_res/total.total_confs)<<endl;
            }
        }
     
        if ((std::fabs(total.dual-global_vars.old_dual) < opts.dualimprovthres) && opts.algorithm != 2) { 
           global_vars.converged = true;
           cout<< "Dual Objective: " << total.dual<< " "<<"Primal Objective: "<<total.primal<<endl;
           cout<<" Number of iteration at convergence:"<<context.iteration() +2 <<endl;
        }
        
        if((sqrt(total.dual_res/total.total_confs) < opts.dualimprovthres
             && sqrt(total.primal_res/total.total_confs) < opts.dualimprovthres)
      &&(std::fabs(total.dual-global_vars.old_dual) < opts.dualimprovthres)&& opts.algorithm == 2){ 
           global_vars.converged = true;
           cout<< "Dual Objective: " << total.dual<< " "<<"Primal Objective: "<<total.primal<<endl;
           cout<<" Number of iteration at convergence:"<<context.iteration() +2 <<endl;       
        }
       
        global_vars.old_dual = total.dual;
        global_vars.sq_norm_g = total.sum_sq_norm_g;
    
        if (opts.history_file != "\0"){
            global_vars.history[0].push_back(context.iteration()+2);
            global_vars.history[1].push_back(global_vars.timer.current_time());
            global_vars.history[2].push_back(total.dual);
            global_vars.history[3].push_back(total.primal);
        }
    } 
 
}

/* end of aggregator functions */

 

#endif
