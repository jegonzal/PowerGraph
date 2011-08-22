#include <iostream>
#include <iomanip>
#include <fstream>


#include <stdint.h>
#include <vector>
#include <map>

#include <graphlab.hpp>



#include "corpus.hpp"
#include "lda.hpp"



#include <graphlab/macros_def.hpp>

class lda_update : public gl::iupdate_functor {
  void operator()(gl::iscope& scope, gl::icallback& callback) {   
    // Do the LDA update for all the words in this document
    
  }
}; 








void load_graph(graph_type& graph, const corpus& data);

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");


  std::string dictionary_fname("dictionary.txt");
  std::string counts_fname("counts.tsv");
  size_t ntopics(50);
  size_t niters(10);
  double alpha(50.0/double(ntopics));
  double beta(0.1);
  std::string alg_str("gibbs");
  std::string llik_fname("llik.txt");
  

  // Setup the parser
  graphlab::command_line_options
    clopts("Apply the LDA model to estimate topic "
           "distributions for each document.");
  clopts.attach_option("dictionary",
                       &dictionary_fname, dictionary_fname,
                       "Dictionary file");
  clopts.attach_option("counts", 
                       &counts_fname, counts_fname, 
                       "Counts file");
  clopts.attach_option("alg",
                       &alg_str, alg_str, 
                       "Algorithm {gibbs, annealing, opt}");
  clopts.attach_option("ntopics", 
                       &ntopics, ntopics, "Number of topics");
  clopts.attach_option("niters",
                       &niters, niters, "Number of iterations");
  clopts.attach_option("alpha",
                       &alpha, alpha, "Alpha prior");
  clopts.attach_option("beta", 
                       &beta, beta, "Beta prior");
  clopts.attach_option("llik_fname",
                       &llik_fname, llik_fname, 
                       "Log-likelihood file.");
  // Parse the command line input
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing input." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Loading the corpus." << std::endl;
  corpus corpus(dictionary_fname, counts_fname);
  std::cout << "Number of words:   " << corpus.nwords << std::endl
            << "Number of docs:    " << corpus.ndocs << std::endl
            << "Number of tokens:  " << corpus.ntokens << std::endl
            << "Algorithm:         " << alg_str << std::endl
            << "Ntopics:           " << ntopics << std::endl
            << "Alpha:             " << alpha   << std::endl
            << "Beta:              " << beta    << std::endl;

  std::cout << "Seeding Generator: " << std::endl;
  graphlab::random::nondet_seed();
  std::cout << "Shuffling corpus: " << std::endl;
  corpus.shuffle_tokens();

  // Setup the core
  gl::core core;
  core.set_options(clopts);
  std::cout << "Building Graph" << std::endl;
  load_graph(core.graph(), corpus);
  std::cout << "Finished loading graph" << std::endl;
  return EXIT_SUCCESS;
} // end of main




void load_graph(graph_type& graph, const corpus& data) {
  // Construct all the vertices
  const gl::vertex_id nverts = data.nwords + data.ndocs;
  graph.resize(nverts);
  for(gl::vertex_id vid = 0; vid < data.ndocs; ++vid) {
    graph.vertex_data(vid).type = DOCUMENT;    
  }
  for(gl::vertex_id vid = data.ndocs; vid < nverts; ++vid) {
    graph.vertex_data(vid).type = WORD;
  }
  const size_t word_offset = data.ndocs;
  typedef corpus::token token_type;
  foreach(const token_type& tok, data.tokens) {
    document_vertex_data::word_topic_pair wt_pair(tok.word);    
    graph.vertex_data(tok.doc).doc.words.push_back(wt_pair);
    graph.add_edge(tok.doc, tok.word + word_offset);
  }
}; // end of load_graph







