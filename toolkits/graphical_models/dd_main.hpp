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


#ifndef __DD_MAIN_H__
#define __DD_MAIN_H__

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <sys/time.h>
//#include <ctime
#include <getopt.h>

//#include "utils.h"
#include "utils.hpp"
#include "dd_grlab.hpp"
#include "ad3_qp.hpp"


/////////////////////////////////////////////////////////////////////////
// Load the UAI file. Each factor as a different vertex
void loadUAIfile(graphlab::distributed_control& dc, graph_type& graph, string graph_file, int& nodes) 
{  
    // Not sure why this is needed
    dc.barrier();
    // Open file
    ifstream in(graph_file.c_str());
     
    //CHECK(in.good(),"Could not open file: "+graph_file);
    CHECK(in.good());
    
    // Read type of network
    string name; 
    
    in >> name; 
    //CHECK(name.compare("MARKOV")==0, "Only Markov networks are supported. Are you sure this is a typeUAI energy file?");
    CHECK(name.compare("MARKOV")==0);
    
    // Read size of graph
    int nnodes, nfactors;
    in >> nnodes;
    nodes = nnodes;
    //CHECK(nnodes>0, "No. of nodes can't be negative. Are you sure this is a typeUAI energy file?");
    CHECK(nnodes>0);    
    // Read node cardinalities
    vector<int> cardinalities(nnodes,0);
    int cardinality_i, sum_of_cardinalities = 0;
    for (int i = 0; i != nnodes; ++i) 
    {
        in >> cardinality_i;
        cardinalities[i] = cardinality_i;
        sum_of_cardinalities += cardinality_i;
       
        //CHECK(in.good(), "Could not finish reading cardinalities. Are you sure this is a typeUAI energy file?");
        CHECK(in.good());
    }

    int vid = 0;
    if(opts.algorithm != 0){
       for(int i = 0; i < nnodes; i++){                      //temporary .. put condition
           vertex_data vdata;
           vdata.factor_type = VAR; 
           vdata.nvars = 1;
           vdata.cards.resize(1, cardinalities[i]);
           vdata.potentials.setZero(cardinalities[i]);
           vdata.beliefs.setConstant(cardinalities[i], 0.5);
           graph.add_vertex(vid, vdata);
           vid++;
       }
    }
    // Read no. of factors
    in >> nfactors;
    
    //factor_size.resize(nfactors); factor_id.resize(nfactors);
    vector<int> factor_size(nfactors,0); //vector<int> factor_id(nfactors,0); 
    vector< vector<int> > factor_memb; factor_memb.resize(nfactors);
    int temp1, temp2;
    
    // Loop and read factor members
    for (int i=0; i!=nfactors; ++i) 
    {
        in >> temp1;
        factor_size[i] = temp1; 
        
        factor_memb[i].resize(temp1);
        for (int j=0; j!=temp1; ++j) 
        {
            in >> temp2;
            factor_memb[i][j] = temp2;
        }
        
        //CHECK(in.good(), "Could not finish reading cardinalities. Are you sure this is a typeUAI energy file?");
        CHECK(in.good());
    }
    
    if (opts.verbose > 1)
        cout 
        << "Finished Reading UAI-Preamble:"
        << " #Nodes = " << nnodes 
        << ", #Factors = "<< nfactors 
        << ", Average Cardinality = " << double(sum_of_cardinalities)/nfactors
        << "\n";
        
        
    // Now read factor potentials
    for (int i=0; i!=nfactors; ++i) 
    {
        int cardprod; double potential_value; //, energy;
        in >> cardprod;
        
        vertex_data vdata;        
        vdata.nvars = factor_size[i];
        if (vdata.nvars > 1) {
          vdata.degree = vdata.nvars; // Factor degree.
          vdata.factor_type = DENSE;
        }
        else {
          vdata.degree = 1; // Factor degree.
          vdata.factor_type = XOR;
        }
        vdata.cards.resize(factor_size[i]);
        vdata.neighbors.resize(factor_size[i]);
        
        vector<edge_data> edata(factor_size[i]);
        vector<int> varid(factor_size[i]);
        vector<int> card(factor_size[i]);
        
        int cardprod2 = 1;
        for (int j=0; j!=factor_size[i]; ++j) 
        {
            vdata.cards[j] = cardinalities[factor_memb[i][j]];
            vdata.neighbors[j] = factor_memb[i][j]; // afm (check if this was intended!)
            cardprod2 *= vdata.cards[j];
                      
            // Also create edge structs here
            //if (factor_size[i]>1)
           // {
                varid[j] = factor_memb[i][j];
                card[j] = cardinalities[varid[j]];
                edata[j].multiplier_messages.setZero(card[j]);
                edata[j].local_messages.setZero(card[j]);
                edata[j].potentials.setZero(card[j]);
          //  }
        }
        
        //CHECK_EQ(cardprod, cardprod2, "Incorrectly sized factor");
        CHECK_EQ(cardprod, cardprod2);
        
        // Read factor potentials
        vdata.potentials.resize(cardprod);
        vdata.beliefs.resize(cardprod);
        int x_offset = 0;
        for(int x=0; x< vdata.nvars; x++){
            for(int y=0; y<vdata.cards[x]; y++){
               vdata.beliefs[x_offset+y] = 1.0/vdata.cards[x];
            }
            x_offset += vdata.cards[x];
        }         
            
        vdata.factor_beliefs.setConstant(cardprod, 1.0/cardprod);
        for (int k = 0; k != cardprod; ++k) 
        {
            in >> potential_value;
            //energy = Potential2Energy(potential_value);
            
            vdata.potentials[k] = log10(potential_value) ;
        }
        
        //CHECK(in.good(), "Could not finish reading factor tables. Are you sure this is a typeUAI energy file?");

        CHECK(in.good());
         

         vdata.potentials.maxCoeff(&vdata.best_configuration);
        // allocate factors evenly to different machines.
        if (i%dc.numprocs() != dc.procid()) 
            continue;
        
        // If all is well, add vertex and edge
        graph.add_vertex(vid ,vdata);

        if (factor_size[i] > 1 || opts.algorithm > 0) // if not a unary, add edges to unaries
        for (int j=0; j!=factor_size[i]; ++j) 
            graph.add_edge(vid,varid[j],edata[j]);
        
        //after adding everything increment vertex id
        vid++; 
        
        if (opts.verbose > 1)
        {
            cout << "Machine #" << dc.procid() << ", Vertex Id = " << i
            << " with " << vdata.nvars << " variables."; 
            if (factor_size[i] > 1)
            {
                cout << ", Edges = ";
                for (int j=0; j!=factor_size[i]; ++j)             
                    cout << ", (" << i << "," << varid[j] << ")";
            }
            cout << "\n";
            cout << "potential: " << vdata.potentials << "\n";
        }
        
    } // End of reading factors     
   
    dc.barrier();
} // end of loading UAI file

/////////////////////////////////////////////////////////////////////////
// Load the distributed UAI file
bool line_parser(graph_type& graph, const std::string& filename, const std::string& textline) {
    std::stringstream strm(textline);
    graphlab::vertex_id_type vid;
    vertex_data vdata;
    vdata.dual_contrib = 0.0;
    string type;
    strm >> type;
 
     if(type == "v") { 
      vdata.factor_type = VAR;
      vdata.nvars = 1;
      vdata.cards.resize(1);
      strm>>vid;
      strm >> vdata.cards[0];
      vdata.potentials.resize(vdata.cards[0]);
      //vdata.beliefs.setOnes(vdata.cards[0]);
      //vdata.beliefs /= vdata.cards[0];
      vdata.beliefs.setConstant(vdata.cards[0], 0.5);
      vdata.unary_degree.resize(vdata.cards[0], 0);
      //for(int i=0; i< vdata.cards[0]; i++){
      // strm>>vdata.potentials[i];
      //   vdata.potentials[i] = log10(vdata.potentials[i]);
         
      //   }
     //    vdata.potentials.maxCoeff(&vdata.best_configuration);
      graph.add_vertex(vid,vdata);
    }
   else if(type == "d" || type == "u") {
     vdata.factor_type = (type=="d")?DENSE:XOR;
     if(vdata.factor_type == DENSE)
     strm>>vdata.nvars;
     else 
     vdata.nvars = 1;
     strm>>vid;
     vdata.neighbors.resize(vdata.nvars);
     vdata.cards.resize(vdata.nvars);
     int cardprod = 1;
     int cardsum =0;
     for(int i=0; i<vdata.nvars; i++){
        strm>>vdata.neighbors[i]; 
        }
     for(int i=0; i<vdata.nvars; i++){
        strm>>vdata.cards[i]; 
        cardprod *=vdata.cards[i];
        cardsum +=vdata.cards[i];}
     vdata.potentials.setZero(cardprod);
     vdata.beliefs.setOnes(cardprod);
     vdata.beliefs /=cardsum;
     //vdata.beliefs.setConstant(cardprod, 0.5);
     vdata.factor_beliefs.setOnes(cardprod);
     vdata.factor_beliefs /= cardprod;
     for(int i=0; i<cardprod; i++){
        strm>>vdata.potentials[i]; 
        vdata.potentials[i] = log10(vdata.potentials[i]);
        }
        vdata.potentials.maxCoeff(&vdata.best_configuration);
     graph.add_vertex(vid, vdata);
     edge_data edata;
     for(int i=0; i<vdata.nvars; i++)  {
        edata.multiplier_messages.setZero(vdata.cards[i]);
        edata.local_messages.setZero(vdata.cards[i]);
        edata.potentials.setZero(vdata.cards[i]);
        graph.add_edge(vid, vdata.neighbors[i], edata);
     }  
   }
   else if(type == "b") {
   vdata.factor_type = BUDGET;
     strm>>vdata.nvars;
     strm>>vid;
     vdata.neighbors.resize(vdata.nvars);
     vdata.bound_states.resize(vdata.nvars);
     vdata.cards.resize(vdata.nvars);
     vdata.beliefs.setZero(vdata.nvars);
     for(int i=0; i<vdata.nvars; i++){
        strm>>vdata.neighbors[i]; }
        for(int i=0; i<vdata.nvars; i++){
        strm>>vdata.cards[i]; }
     for(int i=0; i<vdata.nvars; i++){
        strm>>vdata.bound_states[i]; }
     strm>>vdata.budget;
     graph.add_vertex(vid, vdata);
     edge_data edata;
     for(int i=0; i<vdata.nvars; i++)  {
        edata.multiplier_messages.setZero(1);
        edata.local_messages.setZero(vdata.cards[i]);
        edata.potentials.setZero(vdata.cards[i]);
        graph.add_edge(vid, vdata.neighbors[i], edata);
     } 
    }
    return true;
 }
/* end of graph loading functions */


////////////////////////////////////////////////////////////////////////////
// Graph transform functions for computing degree and dividing potentials
void compute_degree(graph_type::vertex_type& vertex)
{ 
   vertex.data().degree = vertex.num_out_edges() + vertex.num_in_edges();
   
}

void dist_unary_potentials(graph_type::edge_type& edge)
{ vertex_data& vdata = (edge.source().data().factor_type == VAR)?edge.source().data():edge.target().data();
  edge.data().potentials = vdata.potentials/vdata.degree;
 
}
struct gather_potentials{
   vector<int> degree;
   vec potentials;

    void load(graphlab::iarchive& arc) {
        arc >>degree>>potentials;
    }
    void save(graphlab::oarchive& arc) const {
        arc <<degree<<potentials;
    }

gather_potentials& operator+=(const gather_potentials& other)
{ degree += other.degree;
   potentials += other.potentials;
   return *this;
 } 

};

/* Brief In case of graphs with budget factors degree cannot be determined 
* simply by transform function. A separate vertex program iteration is 
* required. compute_degree_budget computes degree of each vertex and 
* divides unary potentials accordingly. */
struct compute_degree_budget : 
public graphlab::ivertex_program< graph_type, gather_potentials,
graphlab::messages::sum_priority >,
public graphlab::IS_POD_TYPE {

edge_dir_type gather_edges(icontext_type& context,
                               const vertex_type& vertex) const 
    {  return graphlab::ALL_EDGES; 
     };

gather_potentials gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const 
     { 
       const vertex_data& vdata = vertex.data();
       const vertex_type& other_vertex = (edge.source().id() == vertex.id())?edge.target():edge.source();
       vector <int> degree;
       vec potentials;
       
       if(vdata.factor_type == VAR)  {
          potentials.resize(vdata.cards[0]);
          potentials.setZero();
          switch(other_vertex.data().factor_type){
          case XOR : potentials = other_vertex.data().potentials;
       
          case DENSE : degree.resize(vdata.potentials.size(),1);
       
                       break;
          
          case BUDGET : degree.resize(vdata.potentials.size(), 0);
                        int index_neighbor = -1;
                        for(int i=0; i< other_vertex.data().nvars; i++){
                            if(other_vertex.data().neighbors[i] == vertex.id()){
                               index_neighbor = i;
                               break;
                              }
                        }
       
                       degree[other_vertex.data().bound_states[index_neighbor]] = 1;
          }
       }
       else {
         degree.resize(1);
         potentials.resize(1);
       
       }   
       gather_potentials gather_data;
       gather_data.degree  = degree;
       gather_data.potentials = potentials;
       return gather_data; 
     };
       
void apply(icontext_type& context, vertex_type& vertex, const gather_potentials& total)
     { vertex_data& vdata =  vertex.data();
       if(vdata.factor_type == VAR) {
          vdata.unary_degree = total.degree;
          vdata.potentials = total.potentials;
         } 
    
      };

edge_dir_type scatter_edges(icontext_type& context,
                               const vertex_type& vertex) const 
    { 
     return graphlab::ALL_EDGES; };

void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const 
     { const vertex_data& vdata = vertex.data();
       const vertex_type& other_vertex = (edge.source().id() == vertex.id())?edge.target():edge.source();
       if(vdata.factor_type == VAR) {  
          if(other_vertex.data().factor_type != BUDGET) {
             for(int i =0; i< vdata.potentials.size(); i++){
               edge.data().potentials[i] = vdata.potentials[i]/vdata.unary_degree[i];
               }
           }
          else if(other_vertex.data().factor_type == BUDGET) {
               int index_neighbor = -1;
               edge.data().potentials.setZero();
               for(int i=0; i< other_vertex.data().nvars; i++){
                  if(other_vertex.data().neighbors[i] == vertex.id()){
                     index_neighbor = i;
                      break;}
                }
               int state_index = other_vertex.data().bound_states[index_neighbor];
               edge.data().potentials[state_index]  = vdata.potentials[state_index]/vdata.unary_degree[state_index];
              }
            }  
           //cout<<"complete scatter"<<endl;          
       };
};


////////////////////////////////////////////////////////////////////////////
// Graph writer class for saving MAP values. Only unary vertices are saved.
class graph_writer {
public:
std::string save_vertex(graph_type::vertex_type v) {
std::stringstream strm;
if(v.data().factor_type == VAR)
strm << v.id() << "\t" << v.data().best_configuration<< "\n";
return strm.str();
 }
std::string save_edge(graph_type::edge_type e) { return ""; }
 }; /* end of graph_writer */



////////////////////////////////////////////////////////////////////////////
// Functions for running dd , admm

void run_dd_symmetric(graphlab::distributed_control& dc, graph_type& graph, 
               std::string exec_type, graphlab::command_line_options clopts){
    // Define the engine.   
    typedef graphlab::omni_engine<dd_vertex_program_symmetric> engine_type;
    // Instantiate the engine object  
    engine_type engine(dc, graph, opts.exec_type, clopts);
    engine.signal_all();
    graphlab::timer timer;
    // Attach an aggregator to compute primal/dual objective, with periodic interval specified in cmdline argument.
    engine.add_vertex_aggregator<objective>("pd_obj",sum, print_obj); 
    if(!opts.debug){
     engine.aggregate_periodic("pd_obj",opts.agg_time); }
    //The main command. Run graphlab
    engine.start();  
    engine.aggregate_now("pd_obj");
    const double runtime = timer.current_time();    
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    }
    /* end of run_dd_symmetric */
    
    void run_dd_projected(graphlab::distributed_control& dc, graph_type& graph, 
                   std::string exec_type, graphlab::command_line_options clopts){
   
    // Instantiate the engine object
    graph.transform_vertices(compute_degree);
    graph.transform_edges(dist_unary_potentials);
     // Define the engine.    
    typedef graphlab::omni_engine<dd_vertex_program_projected> engine_type;
    engine_type engine(dc, graph, opts.exec_type, clopts);
    engine.signal_all();
    graphlab::timer timer;    
    // Attach an aggregator to compute primal/dual objective, with periodic interval specified in cmdline argument.
    engine.add_vertex_aggregator<objective>("pd_obj",sum, print_obj); 
    if (!opts.debug){
     engine.aggregate_periodic("pd_obj", opts.agg_time);}
    // The main command. Run graphlab
    engine.start();  
    engine.aggregate_now("pd_obj");
    const double runtime = timer.current_time();    
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    }
     /* end of run_dd_projected */
    
    void run_ad3(graphlab::distributed_control& dc, graph_type& graph, 
              std::string exec_type, graphlab::command_line_options clopts){
    // Define the engine.
    typedef  graphlab::omni_engine<compute_degree_budget> transform_engine;
    transform_engine distribute_potentials(dc, graph, opts.exec_type, clopts);
    distribute_potentials.signal_all();
    distribute_potentials.start();
     typedef graphlab::omni_engine<ad3_vertex_program> engine_type;
    // Instantiate the engine object
    engine_type engine(dc, graph, opts.exec_type, clopts);
    engine.signal_all();
    graphlab::timer timer;
    // Attach an aggregator to compute primal/dual objective, with periodic interval specified in cmdline argument.
    engine.add_vertex_aggregator<objective>("pd_obj",sum, print_obj);
    if(!opts.debug){ 
     engine.aggregate_periodic("pd_obj",opts.agg_time); }
    // The main command. Run graphlab
    engine.start();  
    engine.aggregate_now("pd_obj");
    const double runtime = timer.current_time();    
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    
    }
    /* end of run_admm */

 void run_bethe_admm(graphlab::distributed_control& dc, graph_type& graph, 
              std::string exec_type, graphlab::command_line_options clopts){
    // Define the engine.
    typedef  graphlab::omni_engine<compute_degree_budget> transform_engine;
    transform_engine distribute_potentials(dc, graph, opts.exec_type, clopts);
    distribute_potentials.signal_all();
    distribute_potentials.start();
    typedef graphlab::omni_engine<bethe_admm_vertex_program> engine_type;
    // Instantiate the engine object
    engine_type engine(dc, graph, opts.exec_type, clopts);
    engine.signal_all();
    graphlab::timer timer;
    // Attach an aggregator to compute primal/dual objective, with periodic interval specified in cmdline argument.
    engine.add_vertex_aggregator<objective>("pd_obj",sum, print_obj);
    if(!opts.debug){ 
     engine.aggregate_periodic("pd_obj",opts.agg_time); }
    // The main command. Run graphlab
    engine.start();  
    engine.aggregate_now("pd_obj");
    const double runtime = timer.current_time();    
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    }
    /* end of run_bethe_admm */
#endif
