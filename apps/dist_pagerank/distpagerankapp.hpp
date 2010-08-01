/*
 *  pagerankapp.hpp
 *  GraphLab_saved
 *
 *  Created by Aapo Kyrola on 1/21/10.
 *  Copyright 2010 Carnegie Mellon University. All rights reserved.
 *
 */

#ifndef EIGENVECTORAPP_HPP
#define EIGENVECTORAPP_HPP

#include<graphlab.hpp>
#include <graphlab/distributed/graph/distributed_graph.hpp>

using namespace graphlab;

typedef types<blob_distributed_graph> gl_types;


class pagerankapp : public iapp {

public:
  pagerankapp(distributed_control * _dc, std::string inputfile, std::string binoutfile, bool optimize);
  ~pagerankapp();
    
  void start();
  
  
  
private:
  blob_distributed_graph distgraph;
  distributed_control*  dc; 
  
  std::string inputfile;
  std::string binoutfile;
  bool graphoptimize;
  void loadGraph(std::string filename, blob_graph &mgraph);
};


#endif
