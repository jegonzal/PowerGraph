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
#include <algorithm>
#include <vector>
#include <map>
#include <boost/unordered_map.hpp>
#include <time.h>

#include <graphlab.hpp>
#include <graphlab/graph/distributed_graph.hpp>

struct vdata {
  std::vector<size_t> vids;
  void save(graphlab::oarchive& oarc) const {
    oarc << vids;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> vids;
  }
};

typedef graphlab::distributed_graph<vdata, graphlab::empty> graph_type;


void vertex_combine(vdata& a, const vdata& b) {
  for (size_t i = 0;i < b.vids.size(); ++i) a.vids.push_back(b.vids[i]);
}


bool ccoutput_parser(graph_type& graph, const std::string& filename, const std::string& textline) {
  size_t split = textline.find_first_of(",");
  if (split == std::string::npos) return true;
  else {
    std::string t = textline;
    t[split] = 0;
    vdata data;
    data.vids.push_back(atol(t.c_str()));
    graph.add_vertex(atol(t.c_str() + split + 1), data);
    return true;
  }
}

struct size_counter {
  // a map from size to count
  boost::unordered_map<size_t, size_t> counts;

  size_counter() { }

  explicit size_counter(size_t size) {
    counts[size] = 1;
  }

  size_counter& operator+=(const size_counter& other) {
    boost::unordered_map<size_t, size_t>::const_iterator iter = other.counts.begin();
    while(iter != other.counts.end()) {
      counts[iter->first] += iter->second;
      ++iter;
    }
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << counts;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> counts;
  }
};

size_counter absolute_vertex_data(const graph_type::vertex_type& vertex) {
  return size_counter(vertex.data().vids.size());
}

class graph_writer {
public:
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << ":";
    for (size_t i = 0;i < v.data().vids.size(); ++i) {
      strm << v.data().vids[i] << " ";
    }
    strm << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) {
    return "";
  }
};




int main(int argc, char** argv) {
  std::cout << "Connected Component\n\n";

  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  //parse options
  graphlab::command_line_options clopts("Connected Component Stats.");
  std::string graph_dir;
  std::string saveprefix;
  std::string format = "adj";
  clopts.attach_option("graph", graph_dir,
                       "The graph file. This is not optional");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "The graph file format");
  clopts.attach_option("saveprefix", saveprefix,
                       "save location");
  if (!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (graph_dir == "") {
    std::cout << "--graph is not optional\n";
    return EXIT_FAILURE;
  }

  graph_type graph(dc, clopts);

  //load graph
  dc.cout() << "Loading graph in format: "<< format << std::endl;
  graphlab::timer ti;
  graph.set_duplicate_vertex_strategy(vertex_combine);
  graph.load(graph_dir, ccoutput_parser);
  graph.finalize();
  dc.cout() << "Complete Finalization in " << ti.current_time() << std::endl;

  ti.start();
  //take statistics
  size_counter stat = graph.map_reduce_vertices<size_counter>(
      absolute_vertex_data);

  dc.cout() << "graph calculation time is " << ti.current_time() << " sec\n";
  dc.cout() << "RESULT:\nsize\tcount\n";
  for (boost::unordered_map<size_t, size_t>::const_iterator iter = stat.counts.begin();
      iter != stat.counts.end(); iter++) {
    dc.cout() << iter->first << "\t" << iter->second << "\n";
  }
  
  //write results
  if (saveprefix.size() > 0) {
    graph.save(saveprefix, graph_writer(),
        false, //set to true if each output file is to be gzipped
        true, //whether vertices are saved
        false); //whether edges are saved
  }

  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
}

