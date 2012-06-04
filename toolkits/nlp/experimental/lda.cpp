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
size_t ndocs = 0;
size_t global_lag = 100;
double alpha(1.0/double(ntopics));
double beta(0.1);



bool lda_update::use_factorized = false;



void load_graph(graph_type& graph, const corpus& data, 
                size_t iterations);

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "LDA starting\n");


  std::string dictionary_fname("dictionary.txt");
  std::string counts_fname("counts.tsv");
  
  size_t niters(2);
  std::string llik_fname("llik.txt");
  

  // Setup the parser
  graphlab::command_line_options
    clopts("Apply the LDA model to estimate topic "
           "distributions for each document.");
  clopts.set_scheduler_type("sweep");
  clopts.set_scope_type("full");
  clopts.attach_option("dictionary",
                       &dictionary_fname, dictionary_fname,
                       "Dictionary file");
  clopts.attach_option("counts", 
                       &counts_fname, counts_fname, 
                       "Counts file");
  clopts.attach_option("ntopics", 
                       &ntopics, ntopics, "Number of topics");
  clopts.attach_option("factorized", 
                       &lda_update::use_factorized, lda_update::use_factorized, 
                       "use factorized update functor");
  clopts.attach_option("niters",
                       &niters, niters, "Number of iterations");
  clopts.attach_option("lag",
                       &global_lag, global_lag, "Lag between commit");
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
            << "Ntopics:           " << ntopics << std::endl
            << "Alpha:             " << alpha   << std::endl
            << "Beta:              " << beta    << std::endl;

  std::cout << "Seeding Generator: " << std::endl;
  graphlab::random::nondet_seed();
  std::cout << "Shuffling corpus: " << std::endl;
  corpus.shuffle_tokens();

  std::cout << "Setting global variables" << std::endl;
  ndocs = corpus.ndocs;
  nwords = corpus.nwords;

  // Setup the core
  graphlab::core<graph_type, lda_update> core;
  core.set_options(clopts);
  std::cout << "Building Graph" << std::endl;
  load_graph(core.graph(), corpus, niters );
  std::cout << "Finished loading graph" << std::endl;

  // Schedule only the document vertices
  const size_t doc_offset = corpus.nwords;  
  for(size_t i = 0; i < corpus.ndocs; ++i) {
    const graph_type::vertex_id_type doc_vid = doc_offset + i;
    ASSERT_LT(doc_vid, core.graph().num_vertices());
    core.schedule(doc_vid, lda_update(niters));
  }

  // Initialize the global variables
  core.add_global("n_t", size_t(0), ntopics);

  double runtime = core.start();
  std::cout << "Runtime: " << runtime << std::endl;

  return EXIT_SUCCESS;
} // end of main




void load_graph(graph_type& graph, const corpus& data, size_t niters) {
  // Construct all the vertices
  const graph_type::vertex_id_type nverts = data.nwords + data.ndocs;
  std::cout << "Initializing vertices" << std::endl;
  graph.resize(nverts);
  for(graph_type::vertex_id_type word_vid = 0; word_vid < data.nwords; ++word_vid) {
    graph.vertex_data(word_vid).type = WORD;
    graph.vertex_data(word_vid).n_t.resize(ntopics, 0);
  }
  for(graph_type::vertex_id_type doc_vid = data.nwords; doc_vid < nverts; ++doc_vid) {
    graph.vertex_data(doc_vid).type = DOCUMENT;
    graph.vertex_data(doc_vid).n_t.resize(ntopics,0);
  }
  const size_t doc_offset = data.nwords;
  typedef corpus::token token_type;

  // Compute word counts
  std::cout << "Compute word counts." << std::endl;
  std::vector< std::map< word_id_type, count_type> > 
    word_counts(data.ndocs);
  foreach(const token_type& tok, data.tokens)
    word_counts[tok.doc][tok.word]++;

  std::cout << "Adding edges" << std::endl;
  typedef std::pair<word_id_type, count_type> wc_pair_type;
  for(doc_id_type doc = 0; doc < data.ndocs; ++doc) {
    foreach(const wc_pair_type& wc_pair, word_counts[doc]) {
      const edge_data edata(wc_pair.second);
      graph.add_edge(doc + doc_offset, wc_pair.first, edata);
    }
  }
  std::cout << "Finalizing the graph." << std::endl;
  graph.finalize();
}; // end of load_graph







