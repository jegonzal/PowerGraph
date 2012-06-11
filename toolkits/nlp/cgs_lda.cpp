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


#include "cgs_lda_vertex_program.hpp"

#include <graphlab/macros_def.hpp>

std::vector<std::string> dictionary;


// Determine the engine type
typedef graphlab::omni_engine<cgs_lda_vertex_program> engine_type;


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




size_t is_word(const graph_type::vertex_type& vertex) {
  return vertex.num_in_edges() > 0 ? 1 : 0;
}

size_t is_doc(const graph_type::vertex_type& vertex) {
  return vertex.num_out_edges() > 0 ? 1 : 0;
}

void initialize_vertex_data(graph_type::vertex_type& vertex) {
  vertex.data().factor.resize(cgs_lda_vertex_program::NTOPICS);
}



graphlab::empty signal_docs(cgs_lda_vertex_program::icontext_type& context, 
                            graph_type::vertex_type& vertex) {
  if(is_doc(vertex)) context.signal(vertex);
  return graphlab::empty();
} // end of signal_docs

graphlab::empty signal_words(cgs_lda_vertex_program::icontext_type& context, 
                             graph_type::vertex_type& vertex) {
  if(is_word(vertex)) context.signal(vertex);
  return graphlab::empty();
} // end of signal_docs




struct topk {
  typedef std::pair<float, graphlab::vertex_id_type> cw_pair_type;
  std::vector< std::set<cw_pair_type> > top_words;
  size_t nchanges;
  static size_t KVALUE;

  topk(size_t nchanges = 0) : nchanges(nchanges) { }

  topk& operator+=(const topk& other) {
    nchanges += other.nchanges;
    if(other.top_words.empty()) return *this;
    if(top_words.size() < other.top_words.size())
      top_words.resize(other.top_words.size());
    for(size_t i = 0; i < top_words.size(); ++i) {
      // Merge the topk
      top_words[i].insert(other.top_words[i].begin(), 
                          other.top_words[i].end());
      // Remove excess elements        
      while(top_words[i].size() > KVALUE) 
        top_words[i].erase(top_words[i].begin());
    }
    return *this;
  } // end of operator +=

  void save(graphlab::oarchive& arc) const { arc << top_words << nchanges; }
  void load(graphlab::iarchive& arc) { arc >> top_words >> nchanges; }

  static topk map(cgs_lda_vertex_program::icontext_type& context, 
                  const graph_type::vertex_type& vertex) {
    if(is_word(vertex) && !vertex.data().factor.empty()) {
      const graphlab::vertex_id_type wordid = vertex.id();
      const vertex_data& vdata = vertex.data();
      topk ret_value;
      ret_value.top_words.resize(vdata.factor.size());
      float normalizer = 0;
      foreach(double d, vdata.factor) 
        normalizer += (d + cgs_lda_vertex_program::BETA);
      normalizer=1;
      for(size_t i = 0; i < vdata.factor.size(); ++i) {
        const float value = 
          (vdata.factor[i] + cgs_lda_vertex_program::BETA) / normalizer;
        //        const float value = vdata.factor[i];
        const cw_pair_type pair(value, wordid);
        ret_value.top_words[i].insert(pair);
      }
      return ret_value;
    } else { 
      //      context.signal(vertex);
      return topk(vertex.data().nchanges);       
    }
  } // end of map function

  static void finalize(cgs_lda_vertex_program::icontext_type& context,
                       const topk& total) {
    if(context.procid() != 0) return;
    std::cout << "Number of changes: " << total.nchanges << std::endl;
    for(size_t i = 0; i < total.top_words.size(); ++i) {
      std::cout << "Topic " << i << ": ";
      rev_foreach(cw_pair_type pair, total.top_words[i])  {
        ASSERT_LT(pair.second, dictionary.size());
        std::cout << dictionary[pair.second] 
                  << "(" << pair.first << ")" << ", "; 
      }
      std::cout << std::endl;
    }
  } // end of finalize

}; // end of topk struct
size_t topk::KVALUE = 5;










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
    while(std::getline(fin,term).good()) dictionary.push_back(term);
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
    while(std::getline(fin, term).good()) dictionary.push_back(term);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
  std::cout << "Dictionary Size: " << dictionary.size() << std::endl;
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
  size_t interval = 10;
  size_t ntopics = 50;
  clopts.attach_option("dictionary", &dictionary_fname, dictionary_fname,
                       "The file containing the list of unique words");
  clopts.add_positional("dictionary");
  clopts.attach_option("matrix", &matrix_dir, matrix_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");

  clopts.attach_option("ntopics", &ntopics, ntopics,
                       "Number of topics to use.");
  clopts.attach_option("niters", 
                       &cgs_lda_vertex_program::NITERS,
                       cgs_lda_vertex_program::NITERS,
                       "Maximum number of iterations.");
  clopts.attach_option("alpha", 
                       &cgs_lda_vertex_program::ALPHA,
                       cgs_lda_vertex_program::ALPHA,
                       "The document hyper-prior");
  clopts.attach_option("beta",
                       &cgs_lda_vertex_program::BETA,
                       cgs_lda_vertex_program::BETA,                       
                       "The word hyper-prior");
  clopts.attach_option("interval", &interval, interval,
                       "statistics reporting interval");
  if(!clopts.parse(argc, argv) || dictionary_fname.empty() || matrix_dir.empty()) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  ///! set the ntopics global variable
  cgs_lda_vertex_program::set_ntopics(ntopics);

  ///! load the dictionary
  std::cout << dc.procid() << ": Loading the dictionary." << std::endl;
  load_dictionary(dictionary_fname);

  
  ///! load the graph
  std::cout << dc.procid() << ": Loading graph." << std::endl;
  graphlab::timer timer; timer.start();
  graph_type graph(dc, clopts);  
  graph.load(matrix_dir, graph_loader); 
  std::cout << dc.procid() << ": Loading graph. Finished in " 
            << timer.current_time() << std::endl;
  std::cout << dc.procid() << ": Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  graph.transform_vertices(initialize_vertex_data);
  std::cout << dc.procid() << ": Finalizing graph. Finished in " 
            << timer.current_time() << std::endl;
  
  if(dc.procid() == 0){
    std::cout
      << "========== Graph statistics on proc " << dc.procid() 
      << " ==============="
      << "\n Num vertices: " << graph.num_vertices()
      << "\n Num edges: " << graph.num_edges()
      << "\n Num replica: " << graph.num_replicas()
      << "\n Replica to vertex ratio: " 
      << float(graph.num_replicas())/graph.num_vertices()
      << "\n --------------------------------------------" 
      << "\n Num local own vertices: " << graph.num_local_own_vertices()
      << "\n Num local vertices: " << graph.num_local_vertices()
      << "\n Replica to own ratio: " 
      << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
      << "\n Num local edges: " << graph.num_local_edges()
      //<< "\n Begin edge id: " << graph.global_eid(0)
      << "\n Edge balance ratio: " 
      << float(graph.num_local_edges())/graph.num_edges()
      << std::endl;
  }
   
  ///! compute the number of words
  cgs_lda_vertex_program::NWORDS = graph.map_reduce_vertices<size_t>(is_word);
  if(dc.procid() == 0) 
    std::cout << "Number of words: " << cgs_lda_vertex_program::NWORDS;
  ASSERT_EQ(cgs_lda_vertex_program::NWORDS, dictionary.size());

  ///! Setup the engine
  std::cout << dc.procid() << ": Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts, "asynchronous");
  ///! Add an aggregator
  bool success = false;
  success = engine.add_vertex_aggregator<topk>("topk", topk::map, topk::finalize);
  assert(success);
  std::cout << "interval: " << interval << std::endl;
  success = engine.aggregate_periodic("topk", interval);
  assert(success);

  ///! schedule only documents
  engine.map_reduce_vertices<graphlab::empty>(signal_docs);

    

  std::cout << "Running The Collapsed Gibbs Sampler" << std::endl;
  timer.start();
  engine.start();  

  const double runtime = timer.current_time();
  if(dc.procid() == 0) {
    std::cout << "----------------------------------------------------------"
              << std::endl;
    std::cout << "Final Runtime (seconds):   " << runtime 
              << std::endl
              << "Updates executed: " << engine.num_updates() << std::endl
              << "Update Rate (updates/second): " 
              << engine.num_updates() / runtime << std::endl;
  }

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main

























