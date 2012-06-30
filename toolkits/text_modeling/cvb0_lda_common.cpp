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


#include <vector>
#include <algorithm>
#include <graphlab.hpp>


#include "cvb0_lda_common.hpp"

#include <graphlab/macros_def.hpp>

double ALPHA    = 0.1;
double BETA     = 0.1;
size_t NTOPICS  = 50;
size_t NWORDS   = 0;
size_t TOPK     = 5;
size_t INTERVAL = 10;
factor_type GLOBAL_TOPIC_COUNT;
std::vector<std::string> DICTIONARY;
size_t MAX_COUNT = 100;

bool graph_loader(graph_type& graph, const std::string& fname, 
                  const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  const int BASE = 10;
  char* next_char_ptr = NULL;
  graph_type::vertex_id_type doc_id = 
    strtoul(line.c_str(), &next_char_ptr, BASE);
  if(next_char_ptr == NULL) return false;
  const graph_type::vertex_id_type word_id = 
    strtoul(next_char_ptr, &next_char_ptr, BASE);
  if(next_char_ptr == NULL) return false;
  size_t count = 
    strtoul(next_char_ptr, &next_char_ptr, BASE);
  if(next_char_ptr == NULL) return false;
  
  count = std::min(count, MAX_COUNT);

  // since this is a bipartite graph I need a method to number the
  // left and right vertices differently.  To accomplish I make sure
  // all vertices have non-zero ids and then negate the right vertex.
  doc_id += 2; 
  ASSERT_GT(doc_id, 1); 
  doc_id = -doc_id;
  ASSERT_NE(doc_id, word_id);

  // Create an edge and add it to the graph
  graph.add_edge(doc_id, word_id, edge_data(count));
  return true; // successful load
}; // end of graph loader



/** populate the global dictionary */
bool load_dictionary(const std::string& fname) {
  const bool gzip = boost::ends_with(fname, ".gz");
  // test to see if the graph_dir is an hadoop path
  if(boost::starts_with(fname, "hdfs://")) {
    graphlab::hdfs hdfs;
    graphlab::hdfs::fstream in_file(hdfs, fname);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    fin.set_auto_close(false);
    if(gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good()) {
      logstream(LOG_ERROR) << "Error loading dictionary: "
                           << fname << std::endl;
      return false;
    }
    std::string term;
    while(std::getline(fin,term).good()) DICTIONARY.push_back(term);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } else {
    std::cout << "opening: " << fname << std::endl;
    std::ifstream in_file(fname.c_str(), 
                          std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    if (gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good()) {
      logstream(LOG_ERROR) << "Error loading dictionary: "
                           << fname << std::endl;
      return false;
    }
    std::string term;
    std::cout << "Loooping" << std::endl;
    while(std::getline(fin, term).good()) DICTIONARY.push_back(term);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
  std::cout << "Dictionary Size: " << DICTIONARY.size() << std::endl;
  return true;
} // end of load dictionary



bool load_and_initialize_graph(graphlab::distributed_control& dc,
                               graph_type& graph,
                               const std::string& matrix_dir) {  
  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; timer.start();
  graph.load(matrix_dir, graph_loader); 
  dc.cout() << ": Loading graph. Finished in " 
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " 
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Initializing Vertex Data" << std::endl;
  timer.start();
  graph.transform_vertices(initialize_vertex_data);
  dc.cout() << "Finished initializing Vertex Data in " 
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Verivying dictionary size." << std::endl;
  NWORDS = graph.map_reduce_vertices<size_t>(is_word);
  dc.cout()  << "Number of words: " << NWORDS;
  //ASSERT_LT(NWORDS, DICTIONARY.size());

  return true;
} // end of load and initialize graph





























