#include <fstream>

#include <boost/random.hpp>
#include <boost/integer_traits.hpp>


#include "random.hpp"


discrete_rng_type discrete_rng;
uniform_dist_type uniform_dist(0,1);
boost::variate_generator<discrete_rng_type&, uniform_dist_type> 
rand01(discrete_rng, uniform_dist);


size_t multinomial(const std::vector<double>& prb) {
  const double rnd(rand01());
  size_t ind = 0;
  for(double cumsum = prb[ind]; 
      (cumsum < rnd) && (ind+1 < prb.size());
      cumsum += prb[++ind]);
  assert(ind < prb.size());
  return ind;  
}



/**
 * A truely nondeterministic generator
 */
class nondet_generator {
public:
  typedef size_t result_type;
  BOOST_STATIC_CONSTANT(result_type, min_value = 
                        boost::integer_traits<result_type>::const_min);
  BOOST_STATIC_CONSTANT(result_type, max_value = 
                        boost::integer_traits<result_type>::const_max);
  result_type min BOOST_PREVENT_MACRO_SUBSTITUTION () const { return min_value; }
  result_type max BOOST_PREVENT_MACRO_SUBSTITUTION () const { return max_value; }
  
  nondet_generator() {
    rnd_dev.open("/dev/urandom", std::ios::binary | std::ios::in);
    assert(rnd_dev.good());
  }
  // Close the random number generator
  ~nondet_generator() { rnd_dev.close(); }
  // read a size_t from the source
  result_type operator()() {
    // read a machine word into result
    result_type result(0);   
    assert(rnd_dev.good());
    rnd_dev.read(reinterpret_cast<char*>(&result), sizeof(result_type));
    assert(rnd_dev.good());
    //        std::cout << result << std::endl;
    return result;
  }      
private:
  std::ifstream rnd_dev;
};
nondet_generator global_nondet_rng;




void seed_nondet() {
  discrete_rng.seed(global_nondet_rng);
}



