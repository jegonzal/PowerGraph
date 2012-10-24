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

//helper function
float myrand() {
  return static_cast<float>(rand() / (RAND_MAX + 1.0));
}

//helper function to return hash for Flajolet & Martin bitmask
size_t hash_value() {
  size_t ret = 0;
  while (myrand() < 0.5) {
    ret++;
  }
  return ret;
}

const size_t DUPULICATION_OF_BITMASKS = 20;

struct vdata {
  //use two bitmasks for consistency
  std::vector<std::vector<bool> > bitmask1;
  std::vector<std::vector<bool> > bitmask2;
  //indicate which is the bitmask for reading (or writing)
  bool odd_iteration;
  vdata() :
      bitmask1(), bitmask2(), odd_iteration(true) {
  }
  //for exact counting (but needs large memory)
  void create_bitmask(size_t id) {
    std::vector<bool> mask1(id + 2, 0);
    mask1[id] = 1;
    bitmask1.push_back(mask1);
    std::vector<bool> mask2(id + 2, 0);
    mask2[id] = 1;
    bitmask2.push_back(mask2);
  }
  //for approximate Flajolet & Martin counting
  void create_hashed_bitmask(size_t id) {
    for (size_t i = 0; i < DUPULICATION_OF_BITMASKS; ++i) {
      size_t hash_val = hash_value();
      std::vector<bool> mask1(hash_val + 2, 0);
      mask1[hash_val] = 1;
      bitmask1.push_back(mask1);
      std::vector<bool> mask2(hash_val + 2, 0);
      mask2[hash_val] = 1;
      bitmask2.push_back(mask2);
    }
  }
  void save(graphlab::oarchive& oarc) const {
    size_t num = bitmask1.size();
    oarc << num;
    for (size_t a = 0; a < num; ++a) {
      size_t size = bitmask1[a].size();
      oarc << size;
      for (size_t i = 0; i < size; ++i)
        oarc << bitmask1[a][i];
      for (size_t i = 0; i < size; ++i)
        oarc << bitmask2[a][i];
    }
  }
  void load(graphlab::iarchive& iarc) {
    size_t num = 0;
    iarc >> num;
    for (size_t a = 0; a < num; ++a) {
      size_t size = 0;
      iarc >> size;
      std::vector<bool> mask1;
      for (size_t i = 0; i < size; ++i) {
        bool element = true;
        iarc >> element;
        mask1.push_back(element);
      }
      bitmask1.push_back(mask1);
      std::vector<bool> mask2;
      for (size_t i = 0; i < size; ++i) {
        bool element = true;
        iarc >> element;
        mask2.push_back(element);
      }
      bitmask2.push_back(mask2);
    }
  }
};

typedef graphlab::distributed_graph<vdata, graphlab::empty> graph_type;

//[vertex_id1] [vertex_id2]
// NOTE: vertex id must start from 0.
//       A vertex that has no edge must not exist.
bool line_parser(graph_type& graph, const std::string& filename,
    const std::string& textline) {
  std::stringstream strm(textline);
  size_t source = 0;
  size_t target = 0;
  strm >> source;
  strm.ignore(1);
  strm >> target;
  if (source != target)
    graph.add_edge(source, target);

  return true;
}

//initialize bitmask
void initialize_vertex(graph_type::vertex_type& v) {
  v.data().create_bitmask(v.id());
}
//initialize bitmask
void initialize_vertex_for_approximation(graph_type::vertex_type& v) {
  v.data().create_hashed_bitmask(v.id());
}

//helper function to compute OR
void bitwise_or(std::vector<std::vector<bool> >& v1,
    const std::vector<std::vector<bool> >& v2) {
  for (size_t a = 0; a < v1.size(); ++a) {
    while (v1[a].size() < v2[a].size()) {
      v1[a].push_back(false);
    }
    for (size_t i = 0; i < v2[a].size(); ++i) {
      v1[a][i] = v1[a][i] || v2[a][i];
    }
  }
}

struct bitmask_gatherer {
  std::vector<std::vector<bool> > bitmask;

  bitmask_gatherer() :
      bitmask() {
  }
  explicit bitmask_gatherer(const std::vector<std::vector<bool> > & in_b) {
    bitmask = in_b;
  }

  //plus is bitwise-or
  bitmask_gatherer& operator+=(const bitmask_gatherer& other) {
    bitwise_or(bitmask, other.bitmask);
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    size_t num = bitmask.size();
    oarc << num;
    for (size_t a = 0; a < num; ++a) {
      size_t size = bitmask[a].size();
      oarc << size;
      for (size_t i = 0; i < size; ++i)
        oarc << bitmask[a][i];
    }
  }
  void load(graphlab::iarchive& iarc) {
    size_t num = 0;
    iarc >> num;
    for (size_t a = 0; a < num; ++a) {
      size_t size = 0;
      iarc >> size;
      std::vector<bool> mask1;
      for (size_t i = 0; i < size; ++i) {
        bool element = true;
        iarc >> element;
        mask1.push_back(element);
      }
      bitmask.push_back(mask1);
    }
  }
};

//The next bitmask b(h + 1; i) of i at the hop h + 1 is given as:
//b(h + 1; i) = b(h; i) BITWISE-OR {b(h; k) | source = i & target = k}.
class one_hop: public graphlab::ivertex_program<graph_type, bitmask_gatherer>,
    public graphlab::IS_POD_TYPE {
private:
public:
  one_hop() {
  }

  void init(icontext_type& context, const vertex_type& vertex,
      const message_type& msg) {
  }

  //we are going to gather on out edges
  edge_dir_type gather_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::OUT_EDGES;
  }

  //for each edge gather the bitmask of the edge
  bitmask_gatherer gather(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
    if (vertex.data().odd_iteration) {
      return bitmask_gatherer(edge.target().data().bitmask2);
    } else {
      return bitmask_gatherer(edge.target().data().bitmask1);
    }
  }

  //get bitwise-ORed bitmask and switch bitmasks
  void apply(icontext_type& context, vertex_type& vertex,
      const gather_type& total) {
    if (vertex.data().odd_iteration) {
      if (total.bitmask.size() > 0)
        bitwise_or(vertex.data().bitmask1, total.bitmask);
      vertex.data().odd_iteration = false;
    } else {
      if (total.bitmask.size() > 0)
        bitwise_or(vertex.data().bitmask2, total.bitmask);
      vertex.data().odd_iteration = true;
    }
  }

  edge_dir_type scatter_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
  void scatter(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
  }
};

// count number of pairs reached in the current hop
// by looking at bitmasks
struct pair_counter {
  size_t count;

  pair_counter() :
      count(0) {
  }
  explicit pair_counter(size_t in_count) :
      count(in_count) {
  }

  pair_counter& operator+=(const pair_counter& other) {
    count += other.count;
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << count;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> count;
  }
};

//count the number of notes reached in the current hop
pair_counter absolute_vertex_data(const graph_type::vertex_type& vertex) {
  if (vertex.data().odd_iteration == false) { //odd_iteration has just finished
    size_t count = 0;
    for (size_t i = 0; i < vertex.data().bitmask1[0].size(); ++i)
      if (vertex.data().bitmask1[0][i])
        count++;
    return pair_counter(count - 1);
  } else {
    size_t count = 0;
    for (size_t i = 0; i < vertex.data().bitmask2[0].size(); ++i)
      if (vertex.data().bitmask2[0][i])
        count++;
    return pair_counter(count - 1);
  }
}

//count the number of notes reached in the current hop with Flajolet & Martin counting method
size_t approximate_pair_number(std::vector<std::vector<bool> > bitmask) {
  float sum = 0.0;
  for (size_t a = 0; a < bitmask.size(); ++a) {
    for (size_t i = 0; i < bitmask[a].size(); ++i) {
      if (bitmask[a][i] == 0) {
        sum += (float) i;
        break;
      }
    }
  }
  return (size_t) (pow(2.0, sum / (float) (bitmask.size())) / 0.77351);
}
//count the number of notes reached in the current hop
pair_counter absolute_vertex_data_for_approximation(
    const graph_type::vertex_type& vertex) {
  if (vertex.data().odd_iteration == false) { //odd_iteration has just finished
    size_t count = approximate_pair_number(vertex.data().bitmask1);
    return pair_counter(count);
  } else {
    size_t count = approximate_pair_number(vertex.data().bitmask2);
    return pair_counter(count);
  }
}

class graph_writer {
public:
  graph_writer() {
  }
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id();
    strm << " " << approximate_pair_number(v.data().bitmask1);
    strm << " " << approximate_pair_number(v.data().bitmask2);
    strm << "\n";
    return strm.str();
  }

  std::string save_edge(graph_type::edge_type e) {
    std::stringstream strm;
    size_t source = e.source().id();
    size_t target = e.target().id();
    strm << source << " " << target << "\n";
    return strm.str();
  }
};

//return number of data
int main(int argc, char** argv) {
  std::cout << "Approximate graph diameter\n\n";
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  std::string datafile;
  float termination_criteria = 0.001;
  bool approximation = true;
  //parse command line
  graphlab::command_line_options clopts(
                "Approximate graph diameter. The input data file is "
                "provided by the --graph argument which is non-optional. "
                "Directions of edges are considered.");
  std::string graph_dir;
  std::string format = "adj";
  std::string exec_type = "synchronous";
  clopts.attach_option("graph", graph_dir,
                       "The graph file. This is not optional");
  clopts.add_positional("graph");
  clopts.attach_option("engine", exec_type,
                       "The engine type synchronous or asynchronous");
  clopts.attach_option("tol", termination_criteria,
                       "The permissible change at convergence.");
  clopts.attach_option("format", format,
                       "The graph file format");
  clopts.attach_option("approximation", approximation,
                       "If true, use Flajolet & Martin bitmask with "
                       "a smaller memory");

  if (!clopts.parse(argc, argv)){
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (graph_dir == "") {
    std::cout << "--graph is not optional\n";
    return EXIT_FAILURE;
  }

  //load graph
  graph_type graph(dc, clopts);
  dc.cout() << "Loading graph in format: "<< format << std::endl;
  graph.load_format(graph_dir, format);
//  graph.load(datafile, line_parser);
  graph.finalize();

  time_t start, end;
  //initialize vertices
  time(&start);
  if (approximation == false)
    graph.transform_vertices(initialize_vertex);
  else
    graph.transform_vertices(initialize_vertex_for_approximation);

  graphlab::graphlab_options ops;
  graphlab::omni_engine<one_hop> engine(dc, graph, exec_type, ops);

  //main iteration
  size_t previous_count = 0;
  size_t diameter = 0;
  for (size_t iter = 0; iter < 100; ++iter) {
    engine.signal_all();
    engine.start();
    pair_counter stat;
    if (approximation == false)
      stat = graph.map_reduce_vertices<pair_counter>(absolute_vertex_data);
    else
      stat = graph.map_reduce_vertices<pair_counter>(
          absolute_vertex_data_for_approximation);
    size_t current_count = stat.count;
    dc.cout() << iter + 1 << "-th hop: " << current_count
        << " edge pairs are reached\n";
    if (iter > 0
        && (float) current_count
            < (float) previous_count * (1.0 + termination_criteria)) {
      diameter = iter;
      dc.cout() << "converge\n";
      break;
    }
    previous_count = current_count;
  }
  time(&end);

  dc.cout() << "graph calculation time is " << (end - start) << " sec\n";
  dc.cout() << "approximate diameter is " << diameter << "\n";

//  const std::string outputname = datafile + "_out";
//  graph.save(
//      outputname,
//      graph_writer(), false, //set to true if each output file is to be gzipped
//      true, //whether vertices are saved
//      false)//whether edges are saved

  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
}

