#ifndef RANDOM
#define RANDOM

#include <algorithm>
#include <vector>
#include <boost/random.hpp>

typedef boost::mt11213b discrete_rng_type;
typedef boost::uniform_real<double> uniform_dist_type;

// Global variables
extern discrete_rng_type discrete_rng;
extern boost::variate_generator<discrete_rng_type&, uniform_dist_type> rand01;


inline std::ptrdiff_t rand_range(std::ptrdiff_t end) {
  return boost::uniform_int<ptrdiff_t>(0, end-1)(discrete_rng);
}

template<typename T>
void shuffle(std::vector<T>& vec) { 
  std::random_shuffle(vec.begin(), vec.end(), rand_range);
}

size_t multinomial(const std::vector<double>& prb); 

void seed_nondet();







#endif
