#ifndef COLLAPSED_GIBBS
#define COLLAPSED_GIBBS

#include <vector>
#include <boost/random.hpp>
// #include <itpp/itbase.h>

#include "corpus.hpp"
#include "matrix.hpp"

enum inference_type {GIBBS, ANNEALING, OPT};


inference_type string_to_inference_type(const std::string& str);
std::string inference_type_to_string(const inference_type& alg);


class collapsed_gibbs {
private:
  const corpus_type* corpus_ptr;
  const size_t ntopics;
  const double alpha, beta;

  typedef uint32_t count_type;
  //  typedef itpp::Mat<uint32_t> mat_type;
  typedef matrix<count_type> mat_type;
  //! topic assignments for each word
  std::vector< topic_id_type > topics;

  //! n_td(t,d) Number of occurences of topic t in document d
  mat_type n_td;

  //! n_wt(w,t) Number of occurences of word w in topic t
  mat_type n_wt;

  //! n_t(t) The total number of words assigned to topic t
  mat_type n_t;

  //! number of times a token was assigned to a new topic
  size_t nchanges;
    
  //! The current temperature
  double temperature_;

  //! the constant used in linear cooling
  double annealing_constant_;

  //! Inference type
  inference_type alg_;

public:

  collapsed_gibbs(const corpus_type& corpus, 
                  const size_t& ntopics,
                  const double& alpha,
                  const double& beta);

  inference_type& alg();
  double& temperature();
  double& annealing_constant();
  
  
  void iterate();

  size_t get_nchanges() const;

  double get_log_likelihood() const;

  void display_top(const size_t& ntop) const;

};





#endif
