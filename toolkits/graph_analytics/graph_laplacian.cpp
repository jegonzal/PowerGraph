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

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <time.h>

#include <graphlab.hpp>
#include <graphlab/graph/distributed_graph.hpp>

std::vector<float> vertex_degrees;

struct vdata{
  float degree;
  vdata() : degree(0.0){}
  void save(graphlab::oarchive& oarc) const {
    oarc << degree;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> degree;
  }
};

struct edata{
  float weight;
  edata() : weight(1.0){}
  explicit edata(const float in_w):
    weight(in_w){}
  void save(graphlab::oarchive& oarc) const {
    oarc << weight;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> weight;
  }
};

typedef graphlab::distributed_graph<vdata, edata> graph_type;


// [vertex_id1] [vertex_id2] [weight]
// NOTE: vertex id should start from 1.
bool line_parser(graph_type& graph, const std::string& filename,
    const std::string& textline) {
  std::stringstream strm(textline);
  size_t source = 0;
  size_t target = 0;
  float weight = 0.0;
  strm >> source;
  strm.ignore(1);
  strm >> target;
  strm.ignore(1);
  strm >> weight;
  if(source != target)
	  graph.add_edge(source, target, edata(weight));

  return true;
}

//calculate vertex degree
class add_rows: public graphlab::ivertex_program<graph_type, float>, public graphlab::IS_POD_TYPE {
public:
	add_rows() {
  }

  void init(icontext_type& context, const vertex_type& vertex,
      const message_type& msg) {
  }

  //gather on all the edges
  edge_dir_type gather_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }

  //for each edge gather the weighted of the edge
  float gather(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
    return edge.data().weight;
  }

  //take inverse square root
  void apply(icontext_type& context, vertex_type& vertex,
      const gather_type& total) {
    vertex.data().degree = total;
  }

  edge_dir_type scatter_edges(icontext_type& context,
      const vertex_type& vertex) const {
      return graphlab::NO_EDGES;
  }
  void scatter(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
  }
};

//take inverse square root
void inverse_square_root(graph_type::vertex_type& v) {
  v.data().degree = 1.0 / sqrt(v.data().degree);
}

//multiply D^-1/2
void multiply_D(graph_type::edge_type& e) {
  const float& d1 = e.source().data().degree;
  const float& d2 = e.target().data().degree;
  e.data().weight = e.data().weight * d1 * d2;
}

//needed for normalization for ratio cut
struct max_degree{
  float degree;
  max_degree(): degree(0.0){}
  explicit max_degree(float in_degree): degree(in_degree){}

  max_degree& operator+=(const max_degree& other){
    degree = std::max(degree, other.degree);
    return *this;
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << degree;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> degree;
  }
};
max_degree create_max_degree(const graph_type::vertex_type& vertex) {
  return max_degree(vertex.data().degree);
}

//for normalization for ratio cut
float normalize_factor = 1.0;
void normalize_weight(graph_type::edge_type& e) {
  e.data().weight /= normalize_factor;
}
void normalize_degree(graph_type::vertex_type& v) {
  v.data().degree /= normalize_factor;
}

struct max_vid{
  size_t vid;
  max_vid(): vid(0){}
  explicit max_vid(size_t in_vid): vid(in_vid){}

  max_vid& operator+=(const max_vid& other){
    vid = std::max(vid, other.vid);
    return *this;
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << vid;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> vid;
  }
};
max_vid absolute_vertex_data(const graph_type::vertex_type& vertex) {
  return max_vid(vertex.id());
}

class graph_writer_normalized_cut {
public:
	graph_writer_normalized_cut(){}
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    size_t vid = v.id();
    if(vid == 0)
    	return "";
    strm << vid << " " << vid << " 2.0\n";
    return strm.str();
  }

  std::string save_edge(graph_type::edge_type e) {
    std::stringstream strm;
    size_t source = e.source().id();
    size_t target = e.target().id();
    float weight = e.data().weight;
    strm << source << " " << target << " " << weight << "\n";
    strm << target << " " << source << " " << weight << "\n";
    return strm.str();
  }
};

class graph_writer_ratio_cut {
public:
	graph_writer_ratio_cut(){}
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    size_t vid = v.id();
    if(vid == 0)
    	return "";
    strm << vid << " " << vid << " " << 5.0 - v.data().degree << "\n";
    return strm.str();
  }

  std::string save_edge(graph_type::edge_type e) {
    std::stringstream strm;
    size_t source = e.source().id();
    size_t target = e.target().id();
    float weight = e.data().weight;
    strm << source << " " << target << " " << weight << "\n";
    strm << target << " " << source << " " << weight << "\n";
    return strm.str();
  }
};


int main(int argc, char** argv) {
  std::cout << "Construct graph Laplacian for graph partitioning.\n\n";

  //parse command line
  std::string graph_dir;
  std::string format = "adj";
  bool normalized_cut = true;
  bool ratio_cut = false;
  graphlab::command_line_options clopts
	("Constructing graph Laplacian for graph partitioning");
  clopts.attach_option("graph", graph_dir,
                       "The graph file. This is not optional. Vertex ids must start from 1 "
                       "and must not skip any numbers.");
  clopts.attach_option("format", format,
                       "The graph file format. If \"weight\" is set, the program will read "
                       "the data file where each line holds [id1] [id2] [weight].");
//  clopts.attach_option("normalized-cut", normalized_cut,
//                       "construct graph laplacian for normalized cut");
//  clopts.attach_option("ratio-cut", ratio_cut,
//                       "construct graph laplacian for ratio cut");
  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (graph_dir == "") {
	std::cout << "--graph is not optional\n";
	return EXIT_FAILURE;
  }
//  if(normalized_cut == true && ratio_cut == true){
//    std::cout << "Both normalized-cut and ratio-cut are true. Ratio cut is selected.\n";
//    normalized_cut = false;
//  }else if(normalized_cut == false && ratio_cut == false){
//    std::cout << "Both normalized-cut and ratio-cut are false. Ratio cut is selected.\n";
//    ratio_cut = true;
//  }

  //load graph
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  graph_type graph(dc);
  dc.cout() << "Loading graph in format: "<< format << std::endl;
  if(format == "weight")
    graph.load(graph_dir, line_parser);
  else
    graph.load_format(graph_dir, format);
  graph.finalize();

  time_t start, end;
  graphlab::omni_engine<add_rows> engine(dc, graph, "sync", clopts);
  engine.signal_all();
  time(&start);
  engine.start();
  if(normalized_cut == true){
	  graph.transform_vertices(inverse_square_root);
	  graph.transform_edges(multiply_D);
  }else if(ratio_cut == true){//normalize weight for ratio cut
//    normalize_factor = graph.map_reduce_vertices<max_degree>(create_max_degree).degree;
//    graph.transform_edges(normalize_weight);
//    graph.transform_vertices(normalize_degree);
  }
  time(&end);

  dc.cout() << "graph calculation time is " << (end - start) << " sec\n";
  dc.cout() << "writing...\n";

  const std::string outputname = graph_dir + ".glap";
  if(normalized_cut == true)
    graph.save(
	  outputname,
	  graph_writer_normalized_cut(), false, //set to true if each output file is to be gzipped
	  true, //whether vertices are saved
	  true);//whether edges are saved
  else if(ratio_cut == true)
    graph.save(
	  outputname,
	  graph_writer_ratio_cut(), false, //set to true if each output file is to be gzipped
	  true, //whether vertices are saved
	  true);//whether edges are saved

  size_t data_num = graph.map_reduce_vertices<max_vid>(absolute_vertex_data).vid;
  graphlab::mpi_tools::finalize();

  //write the number of data
  const std::string datanum_filename = graph_dir + ".datanum";
  std::ofstream ofs(datanum_filename.c_str());
  if(!ofs) {
    std::cout << "can't create file for number of data" << std::endl;
    return EXIT_FAILURE;
  }
  ofs << data_num;

  return EXIT_SUCCESS;
}

