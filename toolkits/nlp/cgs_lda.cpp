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

// #include <graphlab/util/stl_util.hpp>

#include "cgs_lda.hpp"
#include "cgs_lda_vertex_program.hpp"
#include "fast_cgs_lda_vertex_program.hpp"

#include <graphlab/macros_def.hpp>

double ALPHA    = 0.1;
double BETA     = 0.1;
size_t NITERS   = -1;
size_t NTOPICS  = 50;
size_t NWORDS   = 0;
size_t TOPK     = 5;
size_t INTERVAL = 10;
factor_type GLOBAL_TOPIC_COUNT;
std::vector<std::string> DICTIONARY;


bool graph_loader(graph_type& graph, const std::string& fname, 
                  const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  std::stringstream strm(line);
  graph_type::vertex_id_type doc_id(-1), word_id(-1);
  size_t count(0);
  strm >> doc_id >> word_id >> count;
  // since this is a bipartite graph I need a method to number the
  // left and right vertices differently.  To accomplish I make sure
  // all vertices have non-zero ids and then negate the right vertex.
  doc_id += 2; 
  ASSERT_GT(doc_id, 1); 
  doc_id = -doc_id;
  ASSERT_NE(doc_id, word_id);
  // Create an edge and add it to the graph
  graph.add_edge(doc_id, word_id, edge_data(count, NULL_TOPIC)); 
  return true; // successful load
}; // end of graph loader


void initialize_vertex_data(graph_type::vertex_type& vertex) {
  vertex.data().factor.resize(NTOPICS);
}


/** populate the global dictionary */
void load_dictionary(const std::string& fname) {
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
      logstream(LOG_FATAL) << "Error loading dictionary: "
                           << fname << std::endl;
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
      logstream(LOG_FATAL) << "Error opening file:" << fname << std::endl;
    }
    std::string term;
    std::cout << "Loooping" << std::endl;
    while(std::getline(fin, term).good()) DICTIONARY.push_back(term);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
  std::cout << "Dictionary Size: " << DICTIONARY.size() << std::endl;
} // end of load dictionary





int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Run the asynchronous collapsed Gibbs Sampler.";
  graphlab::command_line_options clopts(description);
  std::string matrix_dir; 
  std::string dictionary_fname;
  std::string algorithm = "fast";
  clopts.attach_option("dictionary", &dictionary_fname, dictionary_fname,
                       "The file containing the list of unique words");
  clopts.add_positional("dictionary");
  clopts.attach_option("matrix", &matrix_dir, matrix_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("ntopics", &NTOPICS, NTOPICS,
                       "Number of topics to use.");
  clopts.attach_option("niters", &NITERS, NITERS,
                       "Maximum number of iterations.");
  clopts.attach_option("alpha", &ALPHA, ALPHA,
                       "The document hyper-prior");
  clopts.attach_option("beta", &BETA, BETA,                       
                       "The word hyper-prior");
  clopts.attach_option("topk", &TOPK, TOPK,
                       "The number of words to report");
  clopts.attach_option("interval", &INTERVAL, INTERVAL,
                       "statistics reporting interval");
  clopts.attach_option("algorithm", &algorithm, algorithm,
                       "algorithm to use (fast, consistent)");
  if(!clopts.parse(argc, argv) || dictionary_fname.empty() || matrix_dir.empty()) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  ///! Initialize global variables
  GLOBAL_TOPIC_COUNT.resize(NTOPICS);
  load_dictionary(dictionary_fname); 

  
  ///! load the graph
  dc.cout() << ": Loading graph." << std::endl;
  graphlab::timer timer; timer.start();
  graph_type graph(dc, clopts);  
  graph.load(matrix_dir, graph_loader); 
  dc.cout() << ": Loading graph. Finished in " 
            << timer.current_time() << std::endl;
  dc.cout() << ": Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  graph.transform_vertices(initialize_vertex_data);
  dc.cout() << ": Finalizing graph. Finished in " 
            << timer.current_time() << std::endl;
  ///! compute the number of words
  NWORDS = graph.map_reduce_vertices<size_t>(is_word);
  dc.cout()  << "Number of words: " << NWORDS;
  ASSERT_EQ(NWORDS, DICTIONARY.size());

  ///! Run the algorithm
  if(algorithm == "fast") {
    run_fast_cgs_lda(dc, graph, clopts);
  } else {
    run_cgs_lda(dc, graph, clopts);
  }


  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main

























