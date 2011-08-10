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
#include "pagerank.hpp"



#include <graphlab/macros_def.hpp>




std::ostream& operator<<(std::ostream& out, const edge_data& edata) {
  out << "(w: " << edata.weight << ", ov: " << edata.old_source_value << ")";
  return out;
}


std::ostream& operator<<(std::ostream& out, const vertex_data& vdata) {
  out << "(v: " << vdata.value << ", w: " << vdata.self_weight << ")";
  return out;
}



void save_pagerank(const std::string& fname,
                   const pagerank_graph& graph) {
  std::ofstream fout;
  fout.open(fname.c_str());
  fout << std::setprecision(10);
  for(gl::vertex_id vid = 0; vid < graph.num_vertices(); ++vid) {
    fout << graph.vertex_data(vid).value << "\n";
  }
  fout.close();
} // end of save_pagerank

void save_edges_as_tsv(const std::string& fname,
                       const pagerank_graph& graph) {
  std::ofstream fout; 
  fout.open(fname.c_str()); 
  fout << std::setprecision(10);
  for(gl::vertex_id vid = 0; vid < graph.num_vertices(); ++vid) {
    fout << vid << '\t' << vid << '\t' 
         << graph.vertex_data(vid).self_weight << "\n";
    foreach(gl::edge_id eid, graph.out_edge_ids(vid)) {
      fout << vid << '\t' << graph.target(eid) << '\t' 
           << graph.edge_data(eid).weight << "\n";
    }
  }
  fout.close();
} // end of write graph as tsv



void normalize_graph(pagerank_graph& graph) {
  std::cout 
    << "Finalizing graph." << std::endl
    << "\t This is required for the locking protocol to function correctly"
    << std::endl;
  graph.finalize();
  std::cout << "Finished finalization!" << std::endl;

  std::cout << "Normalizing out edge weights." << std::endl;
  // This could be done in graphlab but the focus of this app is
  // demonstrating pagerank
  for(gl::vertex_id vid = 0; 
      vid < graph.num_vertices(); ++vid) {
    vertex_data& vdata = graph.vertex_data(vid);
    // Initialze with self out edge weight
    double sum = vdata.self_weight;
    const gl::edge_list& out_eids = graph.out_edge_ids(vid);
    // Sum up weight on out edges
    for(size_t i = 0; i < out_eids.size(); ++i) {
      const gl::edge_id out_eid = out_eids[i];
      sum += graph.edge_data(out_eid).weight;      
    }
    if (sum == 0) {
      vdata.self_weight = 1.0;
      sum = 1.0; // Dangling page
    }
    assert(sum > 0);
    // divide everything by sum
    vdata.self_weight /= sum;
    for(size_t i = 0; i < out_eids.size(); ++i) {
      const gl::edge_id out_eid = out_eids[i];
      graph.edge_data(out_eid).weight /= sum;
    } 
  }
  std::cout << "Finished normalizing edges." << std::endl;

  //Normalize the vertices
  //  app_metrics.start_time("normalize");
  std::cout << "Randomizing initial pagerank values." << std::endl;
  double sum = 0;
  for(size_t i = 0; i < graph.num_vertices(); ++i) {
    // graph.vertex_data(i).value = 
    //   graphlab::random::rand01() + 1;
    sum += graph.vertex_data(i).value;
  }
  for(size_t i = 0; i < graph.num_vertices(); ++i) {
    graph.vertex_data(i).value = 
      graph.vertex_data(i).value / sum;
  }
  std::cout << "Finished Randomization." << std::endl;
  // app_metrics.stop_time("normalize");
} // end of normalize_graph



// Creates simple 5x5 graph
void make_toy_graph(pagerank_graph& graph) {
  // Create 5 vertices
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());

	
  // Page 0 links to page 3 only, so weight is 1
  graph.add_edge(0, 3, edge_data(1));
	
  // Page 1 links to 0 and 2
  graph.add_edge(1, 0, edge_data(0.5));
  graph.add_edge(1, 2, edge_data(0.5));
	
  // ... and so on
  graph.add_edge(2, 0, edge_data(1.0/3));
  graph.add_edge(2, 1, edge_data(1.0/3));
  graph.add_edge(2, 3, edge_data(1.0/3));

  graph.add_edge(3, 0, edge_data(0.25));
  graph.add_edge(3, 1, edge_data(0.25));
  graph.add_edge(3, 2, edge_data(0.25));
  graph.add_edge(3, 4, edge_data(0.25));

  graph.add_edge(4, 0, edge_data(0.2));
  graph.add_edge(4, 1, edge_data(0.2));
  graph.add_edge(4, 2, edge_data(0.2));
  graph.add_edge(4, 3, edge_data(0.2));
  // and self edge which must be handled specially from 4 to 4
  graph.vertex_data(4).self_weight = 0.2;

  normalize_graph(graph);
} // end of make_toy_graph





inline void skip_newline(std::ifstream& fin) {
  char next_char = ' ';
  fin.get(next_char);
  ASSERT_EQ(next_char, '\n');  
}


bool load_graph_from_metis_file(const std::string& filename,
                                pagerank_graph& graph) {
  
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
  gl::vertex_id source = 0 ;
  size_t edges_processed = 0;
  for(; source < nverts; ++source) {
    while(fin.peek() != '\n') {
      ASSERT_TRUE(fin.good());
      gl::vertex_id target = 0;
      fin >> target; 
      edges_processed++;
      ASSERT_GT(target, 0);
      // decrement the value since starting value is 1 not zero
      target--; 
      ASSERT_LT(target, graph.num_vertices());     
      if(source != target) {
        graph.add_edge(source, target);
      } else {
        // add the self edge by updating the vertex weight
        graph.vertex_data(source).self_weight = 1;
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



  std::cout << "Normalizing the Graph." << std::endl;
  normalize_graph(graph);  
  std::cout << "Finished normalizing." << std::endl;
  

  // { // Simple debugging code
  //   const gl::vertex_id vid = 1;
  //   std::cout << "Inspecting: " << vid << std::endl;
  //   std::cout << graph.vertex_data(vid) << std::endl;
  //   std::cout << "In edges: " << std::endl;
  //   foreach(gl::edge_id ineid, graph.in_edge_ids(vid)) {
  //     std::cout << '\t' << graph.source(ineid) << ", " 
  //               << graph.edge_data(ineid) << std::endl;
  //   }
  //   std::cout << "Out edges: " << std::endl;
  //   foreach(gl::edge_id outeid, graph.out_edge_ids(vid)) {
  //     std::cout << '\t' << graph.target(outeid) << ", " 
  //               << graph.edge_data(outeid) << std::endl;
  //   }
  //   getchar();
  // }
  

  return true;
} // end of load graph from metis file.



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
bool load_graph_from_tsv_file(const std::string& filename,
                              pagerank_graph& graph) {
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
      const edge_data edata(weight);
      graph.add_edge(source, target, edata);
    } else {
      // add the self edge by updating the vertex weight
      graph.vertex_data(source).self_weight = weight;
    }       
  }
  std::cout 
    << "Finished loading graph with: " << std::endl
    << "\t Vertices: " << graph.num_vertices() << std::endl
    << "\t Edges: " << graph.num_edges() << std::endl;



  normalize_graph(graph);

  return true;
} // end of load graph
