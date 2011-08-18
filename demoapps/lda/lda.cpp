#include <iostream>
#include <iomanip>
#include <fstream>


#include <stdint.h>
#include <vector>
#include <map>

#include <boost/foreach.hpp>



#include <boost/program_options.hpp>

#define foreach BOOST_FOREACH
#define rev_foreach BOOST_REVERSE_FOREACH

#include <random.hpp>
#include <corpus.hpp>
#include <collapsed_gibbs.hpp>








int main(int argc, char** argv) {

  std::string dictionary_fname("dictionary.txt");
  std::string counts_fname("counts.tsv");
  size_t ntopics(50);
  size_t niters(10);
  double alpha(50.0/double(ntopics));
  double beta(0.1);
  std::string alg_str("gibbs");
  std::string llik_fname("llik.txt");
  


  // Parse command line options
  namespace po = boost::program_options;
  po::options_description desc("LDA sampler code");
  desc.add_options()
    ("help", "produce help message")
    ("dictionary", po::value<std::string>(&dictionary_fname)->
     default_value(dictionary_fname), "Dictionary file")
    ("counts", po::value<std::string>(&counts_fname)->
     default_value(counts_fname), "Counts file")
    ("alg", po::value<std::string>(&alg_str)->
     default_value(alg_str), "Algorithm {gibbs, annealing, opt}")
    ("ntopics", po::value<size_t>(&ntopics)->
     default_value(ntopics), "Number of topics")
    ("niters", po::value<size_t>(&niters)->
     default_value(niters), "Number of iterations")
    ("alpha", po::value<double>(&alpha)->
     default_value(alpha), "Alpha prior")
    ("beta", po::value<double>(&beta)->
     default_value(beta), "Beta prior")
    ("llik_fname", po::value<std::string>(&llik_fname)->
     default_value(llik_fname), "Log-likelihood file.");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  if (vm.count("help")) {
    std::cout << desc << "\n";
    return EXIT_FAILURE;
  }

  // Process the inference type string
  inference_type alg = string_to_inference_type(alg_str);

  
  std::cout << "Loading the corpus." << std::endl;
  corpus_type corpus(dictionary_fname, counts_fname);
  std::cout << "Number of words:   " << corpus.nwords << std::endl
            << "Number of docs:    " << corpus.ndocs << std::endl
            << "Number of tokens:  " << corpus.ntokens << std::endl
            << "Algorithm:         " 
            << inference_type_to_string(alg) << std::endl
            << "Ntopics:           " << ntopics << std::endl
            << "Alpha:             " << alpha   << std::endl
            << "Beta:              " << beta    << std::endl;

  std::cout << "Seeding Generator: " << std::endl;
  seed_nondet();
  std::cout << "Shuffling corpus: " << std::endl;
  corpus.shuffle_tokens();


  std::cout << "Constructing Gibbs Sampler: " << std::endl;
  collapsed_gibbs gibbs(corpus, ntopics, alpha, beta);
  

  
  std::ofstream llik_fout(llik_fname.c_str());
  llik_fout.precision(16);
  for(size_t i = 0; i < niters; ++i) {
    // Run gibbs sampling on the first run and then the chosen
    // algorithm on subsequent runs.
    if(i == 0) gibbs.alg() = GIBBS;
    else gibbs.alg() = alg;

    std::cout << "Sampling iteration: " << i << std::endl;
    std::cout << "Using algorithm: " 
              << inference_type_to_string(gibbs.alg()) << std::endl;
    if(alg == ANNEALING) {
      std::cout << "Temperature: " << gibbs.temperature() << std::endl;
    }
    gibbs.iterate();
    gibbs.display_top(5);
    size_t nchanges = gibbs.get_nchanges();
    std::cout << "Number of changes: " << gibbs.get_nchanges() << std::endl
              << "Prop. Changes:     " 
              << double(nchanges)/ corpus.ntokens << std::endl;
    double llik = gibbs.get_log_likelihood();
    std::cout << "Log-likelihood:    " // std::setprecision(8) <<
              <<  llik << std::endl;
    llik_fout << llik << '\t' << nchanges << '\t'
              << gibbs.alg() << std::endl;

  }
  llik_fout.close();


  return EXIT_SUCCESS;
}
