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


#include <iostream>
#include <algorithm>
#include <set>
#include <fstream>


#include <graphlab/util/fast_multinomial.hpp>

#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>



typedef uint32_t vertex_id_type;
typedef std::pair<vertex_id_type, vertex_id_type> edge_type;
typedef std::vector<edge_type> edge_list_type;

void load_edge_list(const std::string& fname,
                    edge_list_type& edge_list) {
  std::ifstream fin(fname.c_str());
  assert(fin.good());
  std::cout << "Loading edge list" << std::endl;
  size_t max_vid = 0;
  while(fin.good()) {
    size_t source(-1), target(-1);
    fin >> source >> target;    
    if(!fin.good()) break;
    if(source != target) 
      edge_list.push_back(edge_type(source, target));
    max_vid = std::max(max_vid, std::max(source, target) );
  }
  fin.close();
  std::cout << "Nverts: " << (max_vid+1) << std::endl;
  std::cout << "Nedges: " << edge_list.size() << std::endl;  
} // end of load edge list

void cumsum(std::vector<size_t>& vec) {
  for(size_t i = 1; i < vec.size(); ++i) vec[i] += vec[i-1];
}

void compute_out_degree(const edge_list_type& edge_list,
                        std::vector<size_t>& degree) {
  std::cout << "Computing max vertex id" << std::endl;
  size_t max_vid = 0;
  for(size_t i = 0; i < edge_list.size(); ++i) 
    max_vid = std::max(max_vid, size_t(std::max(edge_list[i].first, 
                                                edge_list[i].second)));  
  ASSERT_LT(max_vid, size_t(uint32_t(-1)));
  std::cout << "compute degree table" << std::endl;
  degree.clear(); degree.reserve(max_vid+1); degree.resize(max_vid+1, 0);
  for(size_t i = 0; i < edge_list.size(); ++i) ++degree[edge_list[i].first];
} // end of compute out degree

void  build_csr(const edge_list_type& edge_list,
                std::vector<size_t>& degree,
                std::vector<vertex_id_type>& neighbors) {
  std::cout << "Computing the cumulative sum of the degree" << std::endl;
  cumsum(degree);
  std::cout << "nedges: " << degree.back() << std::endl;
  std::cout << "Allocating the neighbor table" << std::endl;
  neighbors.clear(); neighbors.reserve(degree.back()); 
  neighbors.resize(degree.back(),-1);
  std::cout << "Allocating the local counters" << std::endl;
  std::vector<size_t> counters(degree.size(),0);
  for(size_t i = 1; i < counters.size(); ++i) counters[i] = degree[i-1];
  std::cout << "Populating neighbors table." << std::endl;
  for(size_t i = 0; i < edge_list.size(); ++i) {
    const vertex_id_type source(edge_list[i].first), 
      target(edge_list[i].second);
    const size_t index = counters[source]++;
    ASSERT_LT(index, neighbors.size());
    ASSERT_EQ(neighbors[index], size_t(-1));
    neighbors[index] = target;
  }
  std::cout << "validating neighbors table." << std::endl;
  for(size_t i = 0; i < neighbors.size(); ++i) 
    ASSERT_NE(neighbors[i], size_t(-1));  
} // end of compute out degree


void sort_neighbors(const std::vector<size_t>& offsets,
                    std::vector<vertex_id_type>& neighbors) {
  std::cout << "Sorting neighbors" << std::endl;
  std::sort(&neighbors[0], &neighbors[offsets[0]]);
  for(size_t source = 1; source < offsets.size(); ++source) {    
    std::sort(&neighbors[offsets[source-1]], 
              &neighbors[offsets[source]]);
  }
} // end of sort neighbors

void save_adj(const std::string& fname,
              const std::vector<size_t>& offsets,
              const std::vector<vertex_id_type>& neighbors) {
  std::cout << "Saving adjacency information" << std::endl;
  std::ofstream fout(fname.c_str());
  assert(fout.good());
  size_t index = 0;
  for(size_t source = 0; source < offsets.size(); ++source) {    
    fout << source << ' ' << (offsets[source] - index) << ' ';
    for( ; index < offsets[source]; ++index) {
      ASSERT_LT(index, neighbors.size());
      fout << neighbors[index];
      if((index + 1) < offsets[source]) fout <<  ' ';
    }    
    fout << '\n';
  }
  fout.close();
} // end of save adj


void edge_to_adj(const std::string& in_fname,
                const std::string& out_fname) {

  std::vector<size_t> offsets;
  std::vector<vertex_id_type> neighbors;
  edge_list_type edge_list;
  load_edge_list(in_fname, edge_list);
  compute_out_degree(edge_list, offsets);
  build_csr(edge_list, offsets, neighbors);
  sort_neighbors(offsets, neighbors);
  save_adj(out_fname, offsets, neighbors);
  
} // end of edge to adj


void reverse(const std::string& in_fname,
             const std::string& out_fname) {
  std::cout << "Input must be in edge list format" << std::endl;
  std::ifstream fin(in_fname.c_str()); assert(fin.good());
  std::ofstream fout(out_fname.c_str()); assert(fout.good());
  size_t max_vid = 0;
  size_t nedges = 0;
  while(fin.good()) {
    size_t source(-1), target(-1);
    fin >> source >> target;    
    if(!fin.good()) break;
    fout << target << '\t' << source << '\n';
    max_vid = std::max(max_vid, std::max(source, target) );
    nedges++;
  }
  fin.close();
  fout.close();
  std::cout << "Nvertices: " << (max_vid+1) << std::endl;
  std::cout << "Nedges:    " << nedges << std::endl;
} // end of edge to adj





int main(int argc, char** argv) {
  
  graphlab::command_line_options 
    clopts("Generate synthetic graph data.", true);
  std::string in_fname, action, out_fname;
  clopts.attach_option("in", &in_fname, in_fname,
                       "The name of the input graph file."); 
  clopts.add_positional("in");
  clopts.attach_option("action", &action, action,
                       "{edge2adj, reverse}."); 
  clopts.add_positional("action");
  clopts.attach_option("out", &out_fname, out_fname,
                       "The name of the output graph file."); 
  clopts.add_positional("out");

  if(!clopts.parse(argc, argv) || 
     in_fname.empty() || action.empty() || out_fname.empty()) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  if(action == "edge2adj") edge_to_adj(in_fname, out_fname);
  else if(action == "reverse") reverse(in_fname, out_fname);
  else {
    std::cout << "Invalid action: " << action << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // end of main
