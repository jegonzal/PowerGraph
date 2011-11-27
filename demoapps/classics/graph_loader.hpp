/**  
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


#include <iomanip>



#include <graphlab/macros_def.hpp>

template<typename Graph>
void save_graph_as_edge_list(const std::string& fname,
                             const Graph& graph) {
  typedef typename Graph::vertex_id_type vertex_id_type;
  typedef typename Graph::edge_id_type   edge_id_type;
  std::ofstream fout; 
  fout.open(fname.c_str()); 
  fout << std::setprecision(10);
  for(vertex_id_type vid = 0; 
      vid < graph.num_vertices(); ++vid) {
    fout << vid << '\t' << vid << '\t' 
         << graph.vertex_data(vid).self_weight << "\n";
    foreach(edge_id_type eid, graph.out_edge_ids(vid)) {
      fout << vid << '\t' << graph.target(eid) << '\t' 
           << graph.edge_data(eid).weight << "\n";
    }
  }
  fout.close();
} // end of save graph as edge list




template<typename Graph>
bool load_graph(const std::string& filename,
                const std::string& format,
                Graph& graph) {
  if(filename.empty()) { 
    std::cout << "No graph file was provided so we will make a toy graph." 
              << std::endl;
    make_toy_graph(graph);
    return true;
  } else {
    // load the graph from the file
    bool success = false;
    std::cout << "Loading: " << filename << std::endl;
    if(format == "metis") {
      std::cout << "\t using format: " << format << std::endl;
      success = load_graph_from_metis_file(filename, graph);
    } else if(format == "jure") {
      std::cout << "\t using format: " << format << std::endl;
      success = load_graph_from_jure_file(filename, graph);
    } else if(format == "tsv") {
      std::cout << "\t using format: " << format << std::endl;
      success = load_graph_from_tsv_file(filename, graph);
    } else {
      logstream(LOG_WARNING) 
        << "Unsupported format \"" << format << "\"!" << std::endl;
    }
    graph.finalize();
    return success;
  }
} // end of generic load graph


/**
 *
 *           0 ------- 1
 *         /   \       | \
 *        2     3      4  5
 *         \   /       | /
 *           6 ------- 7
 *
 */

template<typename Graph>
void make_toy_graph(Graph& graph) {
  graph.resize(8);
  graph.add_edge(0,1); graph.add_edge(1,0);
  graph.add_edge(0,2); graph.add_edge(2,0);
  graph.add_edge(0,3); graph.add_edge(3,0);
  graph.add_edge(1,4); graph.add_edge(4,1);
  graph.add_edge(1,5); graph.add_edge(5,1);
  graph.add_edge(2,6); graph.add_edge(6,2);
  graph.add_edge(3,6); graph.add_edge(6,3);
  graph.add_edge(4,7); graph.add_edge(7,4);
  graph.add_edge(5,7); graph.add_edge(7,5);
  graph.add_edge(6,7); graph.add_edge(7,6);
} // end of make_toy_graph





inline void skip_newline(std::ifstream& fin) {
  char next_char = ' ';
  fin.get(next_char);
  ASSERT_EQ(next_char, '\n');  
}

template<typename Graph>
bool load_graph_from_metis_file(const std::string& filename,
                                Graph& graph) {
  typedef typename Graph::vertex_id_type vertex_id_type;
  typedef typename Graph::edge_id_type   edge_id_type;
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  size_t nverts = 0, nedges = 0;
  fin >> nverts >> nedges;
  std::cout << "Loading graph with " 
            << nverts << " vertices and "
            << nedges << " edges." << std::endl;
  skip_newline(fin);
  std::string line_buffer;
  graph.resize(nverts);
  vertex_id_type source = 0 ;
  size_t edges_processed = 0;
  for(; source < nverts; ++source) {
    while(fin.peek() != '\n') {
      ASSERT_TRUE(fin.good());
      vertex_id_type target = 0;
      fin >> target; 
      edges_processed++;
      ASSERT_GT(target, 0);
      // decrement the value since starting value is 1 not zero
      target--; 
      ASSERT_LT(target, graph.num_vertices());     
      if(source != target) {
        graph.add_edge(source, target);
      } 
    }
    skip_newline(fin);
  }
  fin.close();
  if(edges_processed == nedges) {
    std::cout << "\tDirected Graph." << std::endl;
  } else {
    std::cout << "\tUndirected Graph." << std::endl;
  } 
  return true;
} // end of load graph from metis file.

template<typename Graph>
bool load_graph_from_jure_file(const std::string& filename,
                               Graph& graph) {
  typedef typename Graph::vertex_id_type vertex_id_type;
  typedef typename Graph::edge_id_type   edge_id_type;
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  // Loop through file reading each line
  while(fin.good() && !fin.eof()) {
    if(fin.peek() == '#') {
      std::string str;
      std::getline(fin, str);
      std::cout << str << std::endl;
      continue;
    }
    size_t source = 0;
    size_t target = 0;
    fin >> source;
    if(!fin.good()) break;
    //  fin.ignore(1); // skip comma
    fin >> target;
    assert(fin.good());
    // Ensure that the number of vertices is correct
    if(source >= graph.num_vertices() ||
       target >= graph.num_vertices())
      graph.resize(std::max(source, target) + 1);
    if(source != target)  graph.add_edge(source, target);
  }
  std::cout 
    << "Finished loading graph with: " << std::endl
    << "\t Vertices: " << graph.num_vertices() << std::endl
    << "\t Edges: " << graph.num_edges() << std::endl;
  return true;
} // end of load graph from jure file




/**
 * Load a graph file specified in the format:
 *
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
template<typename Graph>
bool load_graph_from_tsv_file(const std::string& filename,
                              Graph& graph) {
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  // Loop through file reading each line
  while(fin.good() && !fin.eof()) {
    size_t source = 0;
    size_t target = 0;
    float weight = -1;
    fin >> source;
    if(!fin.good()) break;
    //  fin.ignore(1); // skip comma
    fin >> target;
    assert(fin.good());
    //  fin.ignore(1); // skip comma
    fin >> weight;
    assert(fin.good());
    // Ensure that the number of vertices is correct
    if(source >= graph.num_vertices() ||
       target >= graph.num_vertices())
      graph.resize(std::max(source, target) + 1);
    if(source != target) {
      // Add the edge
      graph.add_edge(source, target);
    } 
  }
  std::cout 
    << "Finished loading graph with: " << std::endl
    << "\t Vertices: " << graph.num_vertices() << std::endl
    << "\t Edges: " << graph.num_edges() << std::endl;  
  return true;
} // end of load graph
  
  
  
// file format from gtgraph file
template<typename Graph>
bool load_graph_from_gtgraph_file(const std::string& fname, 
                                  Graph& graph) {
  typedef typename Graph::vertex_id_type vertex_id_type;
  typedef typename Graph::edge_id_type   edge_id_type;

  const size_t MAGIC_WORD = 0x10102048;
  graphlab::binary_input_stream bin(fname.c_str());
  if(!bin.good()) return false;

  // write it 4B wise?
  uint32_t key = 0; bin.read(key);
  if (key != MAGIC_WORD) {
    std::cout << "Invalid file format!" << std::endl;
    return false;
  }

  bin.read(key); const bool need_back = key;
  assert(!need_back);
  
  uint32_t nverts = 0;
  bin.read(nverts);
  uint32_t nedges = 0;
  bin.read(nedges);

  std::cout << "Nverts:  " << nverts << std::endl;
  std::cout << "Nedges:  " << nedges << std::endl;  

  std::vector<uint32_t> offsets(nverts + 1, -1);
  bin.read_vector(offsets);

  graph.resize(nverts);
  for(size_t i = 0; i < nverts; ++i) {
    const size_t from = offsets[i];
    const size_t to = offsets[i+1]; 
    for(size_t j = from; j < to; ++j) {
      const vertex_id_type from_node = i;
      uint32_t to_node = -1; bin.read(to_node);
      if ((from_node == 2187551) && (to_node == 2868359)) {
        std::cout << "invalid node id!" << std::endl;
      }
      graph.add_edge(from_node, to_node);
    }
  }
  return true;
} // end of load graph from file

#include <graphlab/macros_undef.hpp>
