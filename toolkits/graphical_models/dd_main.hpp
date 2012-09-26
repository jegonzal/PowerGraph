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
//#include "utils.hpp"
#include "dd_grlab.hpp"

/////////////////////////////////////////////////////////////////////////
// Option Struct
struct Options 
{
    double dualimprovthres;
    double pdgapthres;
    int maxiter;
    int verbose;
    
    Options(): dualimprovthres(1e-5), pdgapthres(1e-1), maxiter(100), verbose(0)
    {}
};



/////////////////////////////////////////////////////////////////////////
// Load the UAI file. Each factor as a different vertex
void loadUAIfile(graphlab::distributed_control& dc, graph_type& graph,       
                 string graph_file, Options opts) 
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
        
        //cout << cardinalities[i] << " ";
        //CHECK(in.good(), "Could not finish reading cardinalities. Are you sure this is a typeUAI energy file?");
        CHECK(in.good());
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
    
    if (opts.verbose > 0)
        cout 
        << "Finished Reading UAI-Preamble:"
        << " #Nodes = " << nnodes 
        << ", #Factors = "<< nfactors 
        << ", Average Cardinality = " << sum_of_cardinalities/nfactors
        << "\n";
        
        
    // Now read factor potentials
    for (int i=0; i!=nfactors; ++i) 
    {
        // allocate factors evenly to different machines.
        if (i%dc.numprocs() != dc.procid()) 
            continue;
        
        int cardprod; double potential_value; //, energy;
        in >> cardprod;
        
        vertex_data vdata;
        
        vdata.nvars = factor_size[i];
        vdata.cards.resize(factor_size[i]);
        vdata.xmap.resize(factor_size[i]);
        
        vector<edge_data> edata(factor_size[i]);
        
        int cardprod2 = 1;
        for (int j=0; j!=factor_size[i]; ++j) 
        {
            vdata.cards[j] = cardinalities[factor_memb[i][j]];
            cardprod2 *= vdata.cards[j];
            
            // Also create edge structs here
            if (factor_size[i]>1)
            {
                edata[j].varid = factor_memb[i][j];
                edata[j].card = cardinalities[edata[j].varid];
                edata[j].message.resize(edata[j].card,0);
            }
        }
        
        //CHECK_EQ(cardprod, cardprod2, "Incorrectly sized factor");
        CHECK_EQ(cardprod, cardprod2);
        
        // Read factor potential
        vdata.potential.resize(cardprod);
        for (int k = 0; k != cardprod; ++k) 
        {
            in >> potential_value;
            //energy = Potential2Energy(potential_value);
            
            vdata.potential[k] = potential_value;
        }
        
        //CHECK(in.good(), "Could not finish reading factor tables. Are you sure this is a typeUAI energy file?");
        CHECK(in.good());
        
        // If all is well, add vertex and edges
        graph.add_vertex(i,vdata);
        if (factor_size[i] > 1) // if not a unary, add edges to unaries
            for (int j=0; j!=factor_size[i]; ++j) 
                graph.add_edge(i,edata[j].varid,edata[j]);
        
        if (opts.verbose > 0)
        {
            cout << "Machine #" << dc.procid() << ", Vertex Id = " << i; 
            if (factor_size[i] > 1)
            {
                cout << ", Edges = ";
                for (int j=0; j!=factor_size[i]; ++j)             
                    cout << ", (" << i << "," << edata[j].varid << ")";
            }
            cout << "\n";
        }
        
    } // End of reading factors   
    
    
    
    dc.barrier();
} // end of loading UAI file

#endif