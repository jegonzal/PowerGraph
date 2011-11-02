#include <set>

#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH
#define rev_foreach BOOST_REVERSE_FOREACH

#include <boost/math/special_functions/gamma.hpp>

#include "collapsed_gibbs.hpp"

#include "random.hpp"


inference_type string_to_inference_type(const std::string& str) {
  if(str == "gibbs") return GIBBS;
  else if(str == "annealing") return ANNEALING;
  else if(str == "opt") return OPT;
  else {
    std::cout << "Invalid algorithm type " << str 
              << ", defaulting to GIBBS." << std::endl;
    return GIBBS;
  }
}


std::string inference_type_to_string(const inference_type& alg) {
  switch(alg) {
  case GIBBS: return "gibbs"; break;
  case ANNEALING: return "annealing"; break;
  case OPT: return "opt"; break;
  }
  return "ERROR!";
}



collapsed_gibbs::
collapsed_gibbs(const corpus_type& corpus, 
                const size_t& ntopics,
                const double& alpha,
                const double& beta) : 
  corpus_ptr(&corpus),
  ntopics(ntopics),
  alpha(alpha),
  beta(beta),
  topics(corpus.ntokens, NULL_TOPIC),
  n_td(ntopics, corpus.ndocs, 0),
  n_wt(corpus.nwords, ntopics, 0),
  n_t(ntopics, 1, 0),
  nchanges(0), 
  temperature_(1.1),
  annealing_constant_(0.1),
  alg_(GIBBS) {
  // Initialize the topics table
  // topics.resize(corpus.ntokens, NULL_TOPIC);
}

inference_type& collapsed_gibbs::alg() { return alg_; }
double& collapsed_gibbs::temperature() { return temperature_; }
double& collapsed_gibbs::annealing_constant() {
  return annealing_constant_;
}


void collapsed_gibbs::iterate() {
  // Reset the number of changes
  nchanges = 0;
  const corpus_type& corpus(*corpus_ptr);
  std::vector<double> conditional(ntopics);
  
  std::vector<double> pow_table(10000);
  if(alg_ == ANNEALING) {
    for (size_t i = 0; i < pow_table.size(); ++i)
      pow_table[i] = pow(double(1.0)/(pow_table.size()-1) * i, temperature_);
  }


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
    switch (alg_) { 
    case GIBBS:
      // normalize and then sample
      for(size_t t = 0; t < ntopics; ++t) conditional[t] /= normalizer;
      new_topic = multinomial(conditional);      
      break;
    case ANNEALING:
      for(size_t t = 0; t < ntopics; ++t) conditional[t] /= normalizer;
      normalizer = 0;
      for(size_t t = 0; t < ntopics; ++t)  {
        assert(conditional[t] > 0);
        //normalizer += conditional[t] = pow(conditional[t], temperature_);
        const size_t index(conditional[t]*(pow_table.size()-1.0));
        assert(index < pow_table.size());
        normalizer += conditional[t] = pow_table[index];
      }
      assert(normalizer > 0);
      for(size_t t = 0; t < ntopics; ++t) conditional[t] /= normalizer;
      new_topic = multinomial(conditional);
      break;
    case OPT:
      new_topic = 
        std::max_element(conditional.begin(), conditional.end()) - 
        conditional.begin();
      break;
    } // end of switch
    assert(new_topic < ntopics);

    // Update the topic assignment and counters
    topics[i] = new_topic;
    if(new_topic != old_topic) nchanges++;
    ++n_td(new_topic, d); ++n_wt(w, new_topic), ++n_t(new_topic);
  } // end of loop over tokens

  // If annealing was used advance the temperature
  if(alg_ == ANNEALING) {
    temperature_ += annealing_constant_;
  }
 
  const size_t n_td_sum = n_td.sum();
  const size_t n_wt_sum = n_wt.sum();
  const size_t n_t_sum = n_t.sum();
  //std::cout << n_td_sum << ", " << corpus.ntokens << std::endl;
  //std::cout << n_wt_sum << ", " << corpus.ntokens << std::endl;
  //std::cout << n_t_sum << ", " << corpus.ntokens << std::endl;

  assert(n_td_sum == corpus.ntokens);
  assert(n_wt_sum == corpus.ntokens);
  assert(n_t_sum == corpus.ntokens);

}

size_t collapsed_gibbs::get_nchanges() const { return nchanges; }



double collapsed_gibbs::get_log_likelihood() const {
  using boost::math::lgamma;
  assert(corpus_ptr != NULL);
  const corpus_type& corpus(*corpus_ptr);
  const size_t nwords = corpus.nwords;
  const size_t ndocs  = corpus.ndocs;

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
    for(size_t t = 0; t < ntopics; ++t) {
      llik_topics += lgamma(n_td(t,d) + alpha);
    }
    llik_topics -= lgamma(corpus.ntokens_in_doc[d] + ntopics * alpha);
  }
  return llik_words_given_topics + llik_topics;
}


void collapsed_gibbs::display_top(const size_t& ntop) const {
  assert(ntop > 0);
  assert(corpus_ptr != NULL);
  const corpus_type& corpus(*corpus_ptr);

  typedef std::pair<size_t, word_id_type> cw_pair_type;
  for(size_t t = 0; t < ntopics; ++t) {
    std::set< cw_pair_type > top_words;
    for(size_t w = 0; w < corpus.nwords; ++w) {
      if(top_words.size() < ntop || n_wt(w,t) > top_words.begin()->first) {
        top_words.insert(std::make_pair(n_wt(w,t), w));
        if(top_words.size() > ntop) top_words.erase(top_words.begin());
      }
    }
    std::cout << "----------------" << std::endl;
    rev_foreach(const cw_pair_type& pair, top_words) {
      std::cout << pair.first << "\t\t" << corpus.dictionary.at(pair.second) << std::endl;
    }
  }
    
}






