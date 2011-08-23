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

size_t ntopics = 50;
size_t nwords = 0;
double alpha(1/ntopics);
double beta(0.1);





void load_graph(graph_type& graph, const corpus& data);

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");


  std::string dictionary_fname("dictionary.txt");
  std::string counts_fname("counts.tsv");
  
  size_t niters(10);
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
  for(gl::vertex_id vid = 0; vid < data.nwords; ++vid) {
    graph.vertex_data(vid).type = WORD;
  }
  for(gl::vertex_id vid = data.nwords; vid < nverts; ++vid) {
    graph.vertex_data(vid).type = DOCUMENT;
  }
  const size_t doc_offset = data.nwords;
  typedef corpus::token token_type;

  // Compute word counts
  std::vector< std::map< word_id_type, count_type> > 
    word_counts(data.ndocs);
  foreach(const token_type& tok, data.tokens)
    word_counts[tok.doc][tok.word]++;
  
  typedef std::pair<word_id_type, count_type> wc_pair_type;
  for(doc_id_type doc = 0; doc < data.ndocs; ++doc) {
    foreach(const wc_pair_type& wc_pair, word_counts[doc]) {
      graph.add_edge(doc, wc_pair.first, wc_pair.second); 
    }
  }
}; // end of load_graph







