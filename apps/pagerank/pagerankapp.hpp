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

using namespace graphlab;

struct edge_data {
  float weight;
  float lastusedval;
  float srcvalue;
  float residual;
  edge_data() {}
  edge_data(float _w) {
    weight = _w;
    lastusedval = 0;
    srcvalue = 0;
    residual = 0;
  }  
  void save(graphlab::oarchive &oarc) const{
    serialize(oarc, this, sizeof(edge_data));  
  }
  void load(graphlab::iarchive &iarc) {
    deserialize(iarc, this, sizeof(edge_data));
  }
};

struct vertex_data {
  float value;
  float selfweight;
  vertex_data() {}
  vertex_data(float v) {
    value = v;
  }
  void save(graphlab::oarchive &oarc) const{
    serialize(oarc, this, sizeof(vertex_data));  
  }
  void load(graphlab::iarchive &iarc) {
    deserialize(iarc, this, sizeof(vertex_data));
  }
};


typedef graph<vertex_data, edge_data> pagerank_graph;
typedef types<pagerank_graph> gl_types;



using namespace graphlab;

class pagerankapp : public iapp<pagerank_graph> {

public:
  pagerankapp(std::string inputfile, std::string binoutfile, bool optimize);
  ~pagerankapp();
    
  void start();
  
  /* Used with analyzer_listener */
  virtual global_dumper dump_function();  
  virtual std::vector<std::string> dump_headers();  
  virtual int dump_frequency();
    
    
  pagerank_graph g;

  
  
private:
  gl_types::iengine * graphlab;
  std::string inputfile;
    std::string binoutfile;
  bool graphoptimize;
  void loadGraph(std::string filename);
};


#endif
