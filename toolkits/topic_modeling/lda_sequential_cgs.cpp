/*  
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
#include <iomanip>
#include <fstream>
#include <algorithm>


#include <stdint.h>
#include <vector>
#include <map>

#include <boost/program_options.hpp>
#include <boost/math/special_functions/gamma.hpp>


#include <graphlab.hpp>


#include <graphlab/macros_def.hpp>

typedef uint32_t word_id_type;
typedef uint32_t doc_id_type;
typedef uint16_t topic_id_type;
typedef uint32_t count_type;
#define NULL_TOPIC topic_id_type(-1)

struct token_type {
  word_id_type word;
  doc_id_type doc;
  token_type(const word_id_type& word = 0, const doc_id_type& doc = 0) : 
    word(word), doc(doc) { }
};
std::ostream& operator<<(std::ostream& out, const token_type& tok) {
  return out << "(" << tok.word << ", " << tok.doc << ")";
}


struct corpus_type {
  size_t nwords, ndocs, ntokens;
  std::vector< token_type > tokens;
  std::vector<std::string> dictionary;
  std::vector< word_id_type > ntokens_in_doc;
  
  corpus_type(const std::string& dictionary_fname, 
              const std::string& counts_fname ) : 
    nwords(0), ndocs(0), ntokens(0) {
    dictionary.reserve(20000);
    ntokens_in_doc.reserve(5000);
    tokens.reserve(100000);
    load_dictionary(dictionary_fname);
    load_counts(counts_fname);
  }
  void load_dictionary(const std::string& fname) {
    std::ifstream fin(fname.c_str());
    std::string str;
    while(fin.good()) {
      std::getline(fin, str);
      if(fin.good()) { dictionary.push_back(str); nwords++; }
    }
    fin.close();
  }

  void load_counts(const std::string& fname)  {
    std::ifstream fin(fname.c_str());    
    while(fin.good()) {
      // Read a collection of tokens
      const size_t NULL_VALUE(-1);
      size_t word = NULL_VALUE, doc = NULL_VALUE, count = NULL_VALUE;
      fin >> doc >> word >> count;
      if(fin.good()) {
        assert(word != NULL_VALUE && doc != NULL_VALUE && count != NULL_VALUE);
        // update the doc counter
        ndocs = std::max(ndocs, doc + 1);
        // Assert valid word
        assert(word < nwords);
        // Update the words in document counter
        if(doc >= ntokens_in_doc.size())
          ntokens_in_doc.resize(doc+1, 0);
        ntokens_in_doc[doc] += count;
        // Add all the tokens
        token_type tok; tok.word = word; tok.doc = doc;
        for(size_t i = 0; i < count; ++i) tokens.push_back(tok);
      }
    }
    fin.close();
    ntokens = tokens.size();
  } // end of load counts

  void shuffle_tokens() { graphlab::random::shuffle(tokens); }
}; // end of corpus



template<typename T>
class matrix {
private:
  size_t _rows, _cols;
  std::vector<T> data;

  const size_t linear_index(const size_t& i, const size_t& j) const {
    assert(i < _rows && j < _cols);
    return i + j * _rows;
  }

public:
  matrix(const size_t& rows, const size_t& cols, const T& zero = T(0)) :
    _rows(rows), _cols(cols), data(rows*cols, zero) { };
  const T& operator()(const size_t& i, const size_t& j) const {
    return data[linear_index(i,j)];
  }
  size_t rows() const { return _rows; }
  size_t cols() const { return _cols; }
  T& operator()(const size_t& i, const size_t& j) {
    return data[linear_index(i,j)];
  }
  const T& operator()(const size_t& i) const {
    assert(i < data.size());
    return data[i];
  }
  T& operator()(const size_t& i) {
    assert(i < data.size());
    return data[i];
  }
  void zeros() {
    std::fill(data.begin(), data.end(), T(0));
  }
  void operator+=(const matrix& other) {
    assert(_rows == other._rows);
    assert(_cols == other._cols);
    for(size_t i = 0; i < data.size(); ++i) data[i] += other.data[i];
  }
  T sum() const {
    T z(0);
    for(size_t i = 0; i < data.size(); ++i) z += data[i];
    return z;
  }
}; // end of matrix 
typedef matrix<count_type> mat_type;


class collapsed_gibbs {
public:
  const corpus_type* corpus_ptr;
  const size_t ntopics;
  const double alpha, beta;
 
  std::vector< topic_id_type > topics;
  //! n_td(t,d) Number of occurences of topic t in document d
  mat_type n_td;
  //! n_wt(w,t) Number of occurences of word w in topic t
  mat_type n_wt;
  //! n_t(t) The total number of words assigned to topic t
  mat_type n_t;
  //! number of times a token was assigned to a new topic
  size_t nchanges;

  collapsed_gibbs(const corpus_type& corpus, 
                  const size_t& ntopics,
                  const double& alpha,
                  const double& beta) : 
    corpus_ptr(&corpus), ntopics(ntopics), alpha(alpha),
    beta(beta), topics(corpus.ntokens, NULL_TOPIC),
    n_td(ntopics, corpus.ndocs, 0),
    n_wt(corpus.nwords, ntopics, 0),
    n_t(ntopics, 1, 0),
    nchanges(0) { }
  
  
  void iterate() {
    assert(corpus_ptr != NULL);
    const corpus_type& corpus = *corpus_ptr;
    // Reset the number of changes
    nchanges = 0;
    std::vector<double> conditional(ntopics);

    // Loop over all the tokens
    for(size_t i = 0; i < corpus.ntokens; ++i) {
      // Get the word and document for the ith token
      const word_id_type w = corpus.tokens[i].word;
      const doc_id_type d = corpus.tokens[i].doc;
      const topic_id_type old_topic = topics[i];

      // Remove the word from the current counters
      if(old_topic != NULL_TOPIC) {
        --n_td(old_topic, d); --n_wt(w, old_topic), --n_t(old_topic);
      }

      // Construct the conditional
      double normalizer = 0;
      for(size_t t = 0; t < ntopics; ++t) {
        conditional[t] = (alpha + n_td(t,d)) * (beta + n_wt(w,t)) /
          (beta * corpus.nwords + n_t(t)); 
        normalizer += conditional[t];
      }
      assert(normalizer > 0);

      // Draw a new value
      topic_id_type new_topic = 0;
      // normalize and then sample
      for(size_t t = 0; t < ntopics; ++t) conditional[t] /= normalizer;
      new_topic = graphlab::random::multinomial(conditional);

      // Update the topic assignment and counters
      topics[i] = new_topic;
      if(new_topic != old_topic) nchanges++;
      ++n_td(new_topic, d); ++n_wt(w, new_topic), ++n_t(new_topic);
    } // end of loop over tokens

    const size_t n_td_sum = n_td.sum();
    const size_t n_wt_sum = n_wt.sum();
    const size_t n_t_sum = n_t.sum();
    assert(n_td_sum == corpus.ntokens);
    assert(n_wt_sum == corpus.ntokens);
    assert(n_t_sum == corpus.ntokens);
  }
};


double log_likelihood(const double& alpha, const double& beta,
                      const mat_type& n_td, const mat_type& n_wt) {
  using boost::math::lgamma;
  const size_t ndocs  = n_td.cols();
  const size_t ntopics = n_td.rows();
  const size_t nwords = n_wt.rows();

  mat_type n_t(ntopics,1, 0);
  for(size_t t = 0; t < n_wt.cols(); ++t) 
    for(size_t w = 0; w < n_wt.rows(); ++w) n_t(t) += n_wt(w,t);
      
  
  // Matlab Functions:
  //
  //  llik_w_given_z = ...
  //    ntopics * (gammaln(nwords * beta) - nwords * gammaln(beta)) + ...
  //    sum((sum(gammaln(n_wt + beta)) - gammaln( sum(n_wt) + nwords*beta)));
  //
  //  llik_z = ...
  //    ndocs * (gammaln(ntopics * alpha) - ntopics * gammaln(alpha)) + ...
  //    sum(sum(gammaln(n_td + alpha)) - gammaln(sum(n_td) + ntopics * alpha));

  double llik_words_given_topics = 
    ntopics * (lgamma(nwords * beta) - nwords * lgamma(beta));
  for(size_t t = 0; t < ntopics; ++t) {
    for(size_t w = 0; w < nwords; ++w) {
      llik_words_given_topics += lgamma(n_wt(w,t) + beta);
    }
    llik_words_given_topics -= lgamma(n_t(t) + nwords * beta);
  }
  double llik_topics = ndocs * (lgamma(ntopics * alpha) - ntopics * lgamma(alpha));
  for(size_t d = 0; d < ndocs; ++d) {
    size_t ntokens_in_doc = 0;
    for(size_t t = 0; t < ntopics; ++t) {
      llik_topics += lgamma(n_td(t,d) + alpha);
      ntokens_in_doc += n_td(t,d); 
    }
    llik_topics -= lgamma(ntokens_in_doc + ntopics * alpha);
  }
  return llik_words_given_topics + llik_topics;
} // end of log_likelihood


void display_top(const corpus_type& corpus,
                 const mat_type& n_wt,
                 const size_t& ntop) {
  assert(ntop > 0);
  const size_t nwords = n_wt.rows();
  const size_t ntopics = n_wt.cols();
  typedef std::pair<size_t, word_id_type> cw_pair_type;
  for(size_t t = 0; t < ntopics; ++t) {
    std::set< cw_pair_type > top_words;
    for(size_t w = 0; w < nwords; ++w) {
      if(top_words.size() < ntop || n_wt(w,t) > top_words.begin()->first) {
        top_words.insert(std::make_pair(n_wt(w,t), w));
        if(top_words.size() > ntop) top_words.erase(top_words.begin());
      }
    }
    std::cout << std::endl;
    rev_foreach(const cw_pair_type& pair, top_words) {
      std::cout << corpus.dictionary.at(pair.second) << ", ";
    }
    std::cout << std::endl;
  }
} // end of display top








int main(int argc, char** argv) {

  std::string dictionary_fname("dictionary.txt");
  std::string counts_fname("counts.tsv");
  size_t ntopics(50);
  size_t nburnin(50);
  size_t nsamples(10);
  double alpha(50.0/double(ntopics));
  double beta(0.1);
  size_t topk(20);
  std::string llik_fname("llik.txt");
  std::string doctop_fname("doctop.txt");
  std::string wordtop_fname("wordtop.txt");


  // Parse command line options
  namespace po = boost::program_options;
  po::options_description desc("LDA sampler code");
  desc.add_options()
    ("help", "produce help message")
    ("dictionary", po::value<std::string>(&dictionary_fname)->
     default_value(dictionary_fname), "Dictionary file")
    ("counts", po::value<std::string>(&counts_fname)->
     default_value(counts_fname), "Counts file")
    ("ntopics", po::value<size_t>(&ntopics)->
     default_value(ntopics), "Number of topics")
    ("nburnin", po::value<size_t>(&nburnin)->
     default_value(nburnin), "Number of iterations")
    ("nsamples", po::value<size_t>(&nsamples)->
     default_value(nsamples), "Number of iterations")
    ("alpha", po::value<double>(&alpha)->
     default_value(alpha), "Alpha prior")
    ("beta", po::value<double>(&beta)->
     default_value(beta), "Beta prior")
    ("doctop_fname", po::value<std::string>(&doctop_fname)->
     default_value(doctop_fname), "doctop_fname")
    ("wordtop_fname", po::value<std::string>(&wordtop_fname)->
     default_value(wordtop_fname), "wordtop_fname")
    ("topk", po::value<size_t>(&topk)->
     default_value(topk), "number of top k to show");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  if (vm.count("help")) {
    std::cout << desc << "\n";
    return EXIT_FAILURE;
  }

  if (dictionary_fname.length() == 0 || counts_fname.length() == 0) {
    std::cout << "Both counts and dictionary must be specified" << std::endl;
    std::cout << desc << "\n";
    return EXIT_FAILURE; 
  }
  
  std::cout << "Loading the corpus." << std::endl;
  corpus_type corpus(dictionary_fname, counts_fname);

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


  std::cout << "Constructing Gibbs Sampler: " << std::endl;
  collapsed_gibbs gibbs(corpus, ntopics, alpha, beta);
  
  
  std::ofstream llik_fout(llik_fname.c_str());
  llik_fout.precision(16);

  std::cout << "Starting Burnin" << std::endl;
  for(size_t i = 0; i < nburnin; ++i) {
    std::cout << "Burnin iteration: " << i << std::endl;
    gibbs.iterate();
    std::cout << "Computing top " << topk << " of each topic" << std::endl;
    display_top(corpus, gibbs.n_wt, topk);
    std::cout << "Number of changes: " << gibbs.nchanges << std::endl
              << "Prop. Changes:     " 
              << double(gibbs.nchanges)/ corpus.ntokens << std::endl;
    double llik = log_likelihood(gibbs.alpha, gibbs.beta, gibbs.n_td, gibbs.n_wt);
    std::cout << "Log-likelihood:    " // std::setprecision(8) <<
              <<  llik << std::endl;
    llik_fout << llik << '\t' << gibbs.nchanges << std::endl;

  }
  std::cout << "Finished burnin.  Preparing final sample set." << std::endl;
 
  mat_type n_td(ntopics, corpus.ndocs, 0);
  mat_type n_wt(corpus.nwords, ntopics, 0);
  mat_type n_t(ntopics, 1, 0);

  for(size_t i = 0; i < nsamples; ++i) {
    std::cout << "Sampling iteration: " << i << std::endl;
    gibbs.iterate();
    std::cout << "Number of changes: " << gibbs.nchanges << std::endl
              << "Prop. Changes:     " 
              << double(gibbs.nchanges)/ corpus.ntokens << std::endl;
    std::cout << "Accumulating sample" << std::endl;
    n_td += gibbs.n_td;
    n_wt += gibbs.n_wt;
    n_t  += gibbs.n_t;

    std::cout << "Computing top " << topk << " of each topic" << std::endl;
    display_top(corpus, n_wt, topk);    
    std::cout << "Number of changes: " << gibbs.nchanges << std::endl
              << "Prop. Changes:     " 
              << double(gibbs.nchanges)/ corpus.ntokens << std::endl;
    double llik = log_likelihood(gibbs.alpha, gibbs.beta, gibbs.n_td, gibbs.n_wt);
    std::cout << "Log-likelihood:    " // std::setprecision(8) <<
              <<  llik << std::endl;
    llik_fout << llik << '\t' << gibbs.nchanges << std::endl;

  }
  llik_fout.close();

  std::cout << "Saving doctop: " << doctop_fname << std::endl;
  std::ofstream doctop_fout(doctop_fname.c_str());
  for(size_t d = 0; d < corpus.ndocs; ++d) {
    double normalizer = ntopics * alpha; 
    for(size_t t = 0; t < ntopics; ++t) 
      normalizer += double(n_td(t,d)) / double(nsamples);
    for(size_t t = 0; t < ntopics; ++t) {
      const double value = 
        (double(n_td(t,d))/double(nsamples) + alpha) / normalizer;
      doctop_fout << value << ((t+1 < ntopics)? '\t' : '\n');
    }
  }
  doctop_fout.close();
  
  std::cout << "Saving wordtop: " << wordtop_fname << std::endl;



  return EXIT_SUCCESS;
}
