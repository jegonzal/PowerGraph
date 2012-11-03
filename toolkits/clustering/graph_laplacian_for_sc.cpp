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
#include <algorithm>
#include <vector>
#include <map>
#include <time.h>

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <graphlab.hpp>
#include <graphlab/graph/distributed_graph.hpp>

//shared parameters
float gaussian_kernel_scale_parameter = 0.1;
float threshold_to_discard_small_similarities = 0.0;
size_t number_of_nearest_neighbors= 20;


//data point
struct vertex_data {
  std::vector<float> x;
  float D_ii;
  vertex_data():x(), D_ii(0.0) {};
  explicit vertex_data(const std::vector<float>& x_in) :
      x(x_in), D_ii(0.0) {}

  void save(graphlab::oarchive& oarc) const {
    oarc << x.size();
    for(size_t i=0;i<x.size();++i)
      oarc << x[i];
    oarc << D_ii;
  }

  void load(graphlab::iarchive& iarc) {
    size_t size = 0;
    iarc >> size;
    for(size_t i=0;i<size;++i){
      float temp = 0.0;
      iarc >> temp;
      x.push_back(temp);
    }
    iarc >> D_ii;
  }
};

//similarity
struct edge_data{
  float A_ij;
  bool nearest;
  edge_data() : A_ij(0.0), nearest(false){}
  void save(graphlab::oarchive& oarc) const {
    oarc << A_ij << nearest;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> A_ij >> nearest;
  }
};

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


//[vertex_id] [element1] [element2] [element3] ...
bool line_parser(graph_type& graph, const std::string& filename,
    const std::string& line) {
  if (line.empty()) return true;
  size_t id = 0;
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  vertex_data vtx;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(id) = qi::_1] >> -qi::char_(",") >>
      (qi::double_[phoenix::push_back(phoenix::ref(vtx.x), qi::_1)] % -qi::char_(",") )
      )
     ,
     //  End grammar
     ascii::space);
  if (!success) return false;
  graph.add_vertex(id, vtx);

  for(size_t i=1;i<id;++i){
    graph.add_edge(i, id);
  }

  return true;
}

// helper function to compute similarity between points
float similarity(const std::vector<float>& v1, const std::vector<float>& v2) {
  float ret = 0.0;
  for (size_t i = 0; i < v1.size(); ++i) {
    float tmp = v1[i] - v2[i];
    ret += tmp * tmp;
  }
  return exp(-ret / gaussian_kernel_scale_parameter);
}

//calculate similarities between data points
void calc_similarities(graph_type::edge_type& edata) {
  edata.data().A_ij = similarity(edata.source().data().x, edata.target().data().x);
}


//discard small similarities (Optional)
void discard_small_similarity(graph_type::edge_type& edata) {
  if(edata.data().A_ij < threshold_to_discard_small_similarities)
    edata.data().A_ij = 0.0;
}


//gather T-nearest neighbor (Optional)
struct top_t_similarity{
  std::vector<size_t> ids;
  std::vector<float> sims;
  top_t_similarity(): ids(number_of_nearest_neighbors, std::numeric_limits<size_t>::max()),
      sims(number_of_nearest_neighbors, -1.0){}
  top_t_similarity(size_t id, float sim): ids(number_of_nearest_neighbors, std::numeric_limits<size_t>::max()),
      sims(number_of_nearest_neighbors, -1.0){
    ids[0] = id;
    sims[0] = sim;
  }

  top_t_similarity& operator+=(const top_t_similarity& other){
    std::vector<size_t> new_ids;
    std::vector<float> new_sims;
    size_t pos1=0;
    size_t pos2=0;
    while(pos1+pos2 < number_of_nearest_neighbors){
      if(sims[pos1] >= other.sims[pos2]){
        new_ids.push_back(ids[pos1]);
        new_sims.push_back(sims[pos1]);
        pos1++;
      }else{
        new_ids.push_back(other.ids[pos2]);
        new_sims.push_back(other.sims[pos2]);
        pos2++;
      }
    }
    ids = new_ids;
    sims = new_sims;
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << ids.size();
    for(size_t i=0;i<ids.size();++i)
      oarc << ids[i];
    for(size_t i=0;i<sims.size();++i)
      oarc << sims[i];
  }
  void load(graphlab::iarchive& iarc) {
    ids.clear();
    sims.clear();
    size_t size = 0;
    iarc >> size;
    for(size_t i=0;i<size;++i){
      size_t id = 0;
      iarc >> id;
      ids.push_back(id);
    }
    for(size_t i=0;i<size;++i){
      float sim = 0;
      iarc >> sim;
      sims.push_back(sim);
    }
  }
};

//get T-nearest neighbor and discard others (Optional)
class t_nearest: public graphlab::ivertex_program<graph_type,
  top_t_similarity>, public graphlab::IS_POD_TYPE {
private:
  float threshold;

public:
  t_nearest():threshold(0.0){}

  edge_dir_type gather_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }
  top_t_similarity gather(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
    if(edge.target().id() == vertex.id()){//in edge
      return top_t_similarity(edge.source().id(), edge.data().A_ij);
    }else{//out edge
      return top_t_similarity(edge.target().id(), edge.data().A_ij);
    }
  }

  //assign a cluster, considering the clusters of neighbors
  void apply(icontext_type& context, vertex_type& vertex,
      const gather_type& total) {
    threshold = total.sims[number_of_nearest_neighbors-1];
//    std::cout << vertex.id() << "\t" << total.ids[0] << "-" << total.sims[0] << ", "
//        << total.ids[1] << "-" << total.sims[1] << std::endl;
  }

  edge_dir_type scatter_edges(icontext_type& context,
      const vertex_type& vertex) const {
      return graphlab::ALL_EDGES;
  }
  void scatter(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
    if(edge.data().A_ij >= threshold)
      edge.data().nearest = true;
  }
};

//discard small similarities (Optional)
void make_other_similarities_zero(graph_type::edge_type& edata) {
  if(edata.data().nearest == false)
    edata.data().A_ij = 0.0;
}


//compute sums over rows and then take inverse square root
class calc_degrees: public graphlab::ivertex_program<graph_type,
    float>, public graphlab::IS_POD_TYPE {
public:
  //gather A_ij
  edge_dir_type gather_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }
  float gather(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
    return edge.data().A_ij;
  }

  //assign a cluster, considering the clusters of neighbors
  void apply(icontext_type& context, vertex_type& vertex,
      const gather_type& total) {
    vertex.data().D_ii = 1.0 / sqrt(total);
  }

  edge_dir_type scatter_edges(icontext_type& context,
      const vertex_type& vertex) const {
      return graphlab::NO_EDGES;
  }
  void scatter(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
  }
};

//multiply D^-1/2
void mult_D(graph_type::edge_type& edata) {
  edata.data().A_ij = edata.data().A_ij * edata.source().data().D_ii * edata.target().data().D_ii;
}

struct max_min_similarity{
  float max_sim;
  float min_sim;

  max_min_similarity(): max_sim(0.0), min_sim(0.0){}
  explicit max_min_similarity(float similarity): max_sim(similarity),
      min_sim(similarity){}

  max_min_similarity& operator+=(const max_min_similarity& other){
    if(max_sim < 1.0 && other.max_sim < 1.0){
      max_sim = std::max(max_sim, other.max_sim);
    }else if(other.max_sim < 1.0){
      max_sim = other.max_sim;
    }
    if(min_sim > 0.0 && other.min_sim > 0.0){
      min_sim = std::min(min_sim, other.min_sim);
    }else if(other.min_sim > 0.0){
      min_sim = other.min_sim;
    }
    return *this;
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << max_sim << min_sim;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> max_sim >> min_sim;
  }
};
max_min_similarity absolute_edge_data(const graph_type::edge_type& edge) {
  return max_min_similarity(edge.data().A_ij);
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

class graph_writer {
public:
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    size_t vid = v.id();
    if(vid == 0)
      return "";
    strm << vid << " " << vid << " 1.0\n";
    return strm.str();
  }

  std::string save_edge(graph_type::edge_type e) {
    const float& A_ij = e.data().A_ij;
    std::stringstream strm;
    if(A_ij > 0.0){
      strm << e.source().id() << " " << e.target().id() << " " <<
          A_ij << "\n";
      strm << e.target().id() << " " << e.source().id() << " " <<
          A_ij << "\n";
    }
    return strm.str();
  }
};

int main(int argc, char** argv) {
  std::cout << "construct graph Laplacian for spectral clustering.\n\n";

  //parse command line
  std::string datafile;
  graphlab::command_line_options clopts
    ("Constructing graph Laplacian for spectral clustering");
  clopts.attach_option("data", datafile,
                       "Input file. Each line hold a sample id followed by a white-space or "
                       "comma separated numeric vector. Id should start from 1");
  clopts.attach_option("sigma",  gaussian_kernel_scale_parameter,
                       "Scale parameter for Gaussian kernel.");
  clopts.attach_option("similarity-thres", threshold_to_discard_small_similarities,
                       "Threshold to discard small similarities. ");
  clopts.attach_option("t-nearest", number_of_nearest_neighbors,
                      "Number of nearest neighbors (=t). Will use only the t-nearest similarities "
                      "for each datapoint. If set at 0, will use all similarities.");
  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (datafile == "") {
    std::cout << "--data is not optional\n";
    return EXIT_FAILURE;
  }
  gaussian_kernel_scale_parameter *= 2.0*gaussian_kernel_scale_parameter;

  //construct graph
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  graph_type graph(dc, clopts);
  graph.load(
      datafile,
      line_parser);
  graph.finalize();

  time_t start, end;
  time(&start);
  size_t data_num = graph.map_reduce_vertices<max_vid>(absolute_vertex_data).vid;

  //calculate similarities
  graph.transform_edges(calc_similarities);

  //show the max similarity less than 1 and the min similarity grater than 0
  max_min_similarity stat = graph.map_reduce_edges<max_min_similarity>(absolute_edge_data);
  dc.cout() << "max squared distance(min similarity): "
      << -log(stat.min_sim)*gaussian_kernel_scale_parameter
      << "(" << stat.min_sim << ")\n"
      << "min squared distance(max similarity):"
      << -log(stat.max_sim)*gaussian_kernel_scale_parameter
      << "(" << stat.max_sim << ")\n";


  //if t is set, use only t-nearest similarities
  if(number_of_nearest_neighbors > 0){
    if(number_of_nearest_neighbors > data_num-1)
      number_of_nearest_neighbors = data_num-1;
    dc.cout() << "use only the " << number_of_nearest_neighbors
        << "-nearest similarities for each datapoint\n";
    graphlab::omni_engine<t_nearest> engine_nearest(dc, graph, "sync", clopts);
    engine_nearest.signal_all();
    engine_nearest.start();
    graph.transform_edges(make_other_similarities_zero);
  }
  //if threshold is set, discard similarities less then the threshold
  if(threshold_to_discard_small_similarities > 0.0){
    dc.cout() << "discard small similarities less than "
        << threshold_to_discard_small_similarities << "\n";
    graph.transform_edges(discard_small_similarity);
  }

  //sum elements over rows (calculate the degree matrix D)
  graphlab::omni_engine<calc_degrees> engine(dc, graph, "sync", clopts);
  engine.signal_all();
  engine.start();
  //multiply D
  graph.transform_edges(mult_D);
  time(&end);

  dc.cout() << "graph calculation time is " << (end - start) << " sec\n";
  dc.cout() << "writing...\n";

  //write results
  const std::string outputname = datafile + ".glap";
  graph.save(
      outputname + "_diag",
      graph_writer(), false, //set to true if each output file is to be gzipped
      true, //whether vertices are saved
      false,1); //whether edges are saved
  graph.save(
      outputname + "_other",
      graph_writer(), false, //set to true if each output file is to be gzipped
      false, //whether vertices are saved
      true,1); //whether edges are saved

  //write the number of data
  const std::string datanum_filename = datafile + ".datanum";
  std::ofstream ofs(datanum_filename.c_str());
  if(!ofs) {
    std::cout << "can't create file for number of data" << std::endl;
    return EXIT_FAILURE;
  }
  ofs << data_num;

  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
}

