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

template<typename Source>
void load_edge_list(Source& fin, edge_list_type& edge_list) {
  assert(fin.good());
  size_t max_vid = 0;
  size_t counter = 0;
  while(fin.good()) {
    size_t source(-1), nneighbors(-1);
    fin >> source >> nneighbors;    
    if(!fin.good()) break;
    for(size_t i = 0; i < nneighbors; ++i) {
      size_t target(-1); 
      fin >> target;
      if(source != target) 
        edge_list.push_back(edge_type(source, target));
      max_vid = std::max(max_vid, std::max(source, target) );
    }

    if(counter++ % 100000 == 0) {
      std::cout << counter << "\t" << max_vid << std::endl;
    }

  }
  std::cout << "Nverts: " << (max_vid+1) << std::endl;
  std::cout << "Nedges: " << edge_list.size() << std::endl;  
} // end of load edge list


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


void compute_in_degree(const edge_list_type& edge_list,
                       std::vector<size_t>& degree) {
  std::cout << "Computing max vertex id" << std::endl;
  size_t max_vid = 0;
  for(size_t i = 0; i < edge_list.size(); ++i) 
    max_vid = std::max(max_vid, size_t(std::max(edge_list[i].first, 
                                                edge_list[i].second)));  
  ASSERT_LT(max_vid, size_t(uint32_t(-1)));
  std::cout << "compute degree table" << std::endl;
  degree.clear(); degree.reserve(max_vid+1); degree.resize(max_vid+1, 0);
  for(size_t i = 0; i < edge_list.size(); ++i) ++degree[edge_list[i].second];
} // end of compute in degree


void write_degrees(const std::string& fname, 
                   std::vector<size_t>& degree) {
  std::cout << "Saving degree information: " << fname << std::endl;
  std::ofstream fout(fname.c_str());
  for(size_t i = 0; i < degree.size(); ++i) {
    fout << i << '\t' << degree[i] << '\n';
  }
  fout.close();
} // end of write degree


int main(int argc, char** argv) {  
  graphlab::command_line_options 
    clopts("Generate degree counts.", true);
  std::string graph_fname;
  clopts.attach_option("in", &graph_fname, graph_fname,
                       "The name of the input graph file."); 
  clopts.add_positional("in");

  if(!clopts.parse(argc, argv) || graph_fname.empty()) {
    std::cout << "Error in parsing command line arguments." 
              << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "loading file" << std::endl;
  edge_list_type edges;
  std::ifstream in_file(graph_fname.c_str(), 
                        std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
  // Using gzip filter
  const bool gzip = boost::ends_with(graph_fname, ".gz");      
  if (gzip) fin.push(boost::iostreams::gzip_decompressor());
  fin.push(in_file);
  if(!fin.good()) {
    std::cout << "Error opening file: " << graph_fname << std::endl;
    return false;
  }
  load_edge_list(fin, edges);
  if (gzip) fin.pop();
  fin.pop();
  in_file.close();


  std::cout << "computing out degree" << std::endl;
  std::vector<size_t> degree;
  compute_out_degree(edges, degree);
  write_degrees("out_degree.tsv", degree);
  degree.clear();
  std::cout << "computing in degree" << std::endl;
  compute_in_degree(edges, degree);
  write_degrees("in_degree.tsv", degree);



  


  return EXIT_SUCCESS;
} // end of main
