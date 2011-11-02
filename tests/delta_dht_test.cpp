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



/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


// Test the graph class

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/delta_dht.hpp>

#include <graphlab.hpp>



#include <graphlab/macros_def.hpp>


void basic_test(graphlab::distributed_control& dc) {
  ///! Create the dht 
  graphlab::delta_dht<std::string, int> counts(dc);
  counts.apply_delta("cat", 1);
  counts.apply_delta("dog", 1);
  counts.flush();
  dc.barrier();
  if(dc.procid() == 0) {
    std::cout << "cat: " << counts["cat"] << std::endl;
    std::cout << "dog: " << counts["dog"] << std::endl;
  }
}


///! Create the dht 
typedef size_t word_id_type;

struct topic_vector {
  std::vector<int> n_t;
  topic_vector(size_t ntopics = 0) : n_t(ntopics) { }
  int& operator[](const size_t& index) { return n_t[index]; }
  const int& operator[](const size_t& index) const { return n_t[index]; }
  size_t size() const { return n_t.size(); }
  void resize(const size_t& ntopics) { n_t.resize(ntopics); }


  topic_vector& operator+=(const topic_vector& other) {
    const size_t max_size = std::max(size(), other.size());
    n_t.resize(max_size);
    for(size_t i = 0; i < max_size; ++i) 
      n_t[i] += (i < other.size()? other[i] : 0);
    return *this;
  }

  topic_vector& operator-=(const topic_vector& other) {
    const size_t max_size = std::max(size(), other.size());
    n_t.resize(max_size);
    for(size_t i = 0; i < max_size; ++i) 
      n_t[i] -= (i < other.size()? other[i] : 0);
    return *this;
  } 

  topic_vector operator+(const topic_vector& other) const {
    topic_vector result(*this); result += other;
    return result;
  } 
  topic_vector operator-(const topic_vector& other) const {
    topic_vector result(*this); result -= other;
    return result;
  } 

  void load(graphlab::iarchive& arc) { arc >> n_t; }
  void save(graphlab::oarchive& arc) const { arc << n_t; }
};

typedef graphlab::delta_dht<word_id_type, topic_vector> dictionary_type;

struct functor {
  dictionary_type& dictionary;
  functor(dictionary_type& dictionary) : dictionary(dictionary) { }
  void operator()() {
    std::cout << "Running sampler on machine " << dictionary.procid() << std::endl;
    const size_t n_word_draws = 10000;
    const size_t n_topic_draws = 2;
    const size_t n_topics = 1;
    size_t draws = 0;
    for(size_t i = 0; i < n_word_draws; ++i) {
      const word_id_type word(graphlab::random::gamma(10) * 1000.0);
      // const word_id_type word = 
      //   graphlab::random::fast_uniform<word_id_type>(0, 10000000);
      topic_vector original_vec, vec = dictionary[word];
      if(vec.size() != n_topics) vec.resize(n_topics);
      for(size_t j = 0; j < n_topic_draws; ++j) 
        vec[graphlab::random::fast_uniform<size_t>(0, n_topics-1)]++;
      dictionary.apply_delta(word, vec - original_vec);
      draws++;
    }
    const word_id_type lastword = graphlab::random::fast_uniform<word_id_type>(0, 100); 
    std::cout << "(" << dictionary.procid() << ", " << draws << ", "
              << lastword <<   ")" << std::endl;
    std::cout << "Finished running sampler." << std::endl;
    dictionary.flush();
    std::cout << "Finished flush." << std::endl;
    std::cout << "Dictionary size: " << dictionary.local_size() << std::endl;
  }
};


void large_scale_test(graphlab::distributed_control& dc) {
  ///! Create the dht 
  typedef size_t word_id_type;
  typedef std::vector<int> topic_count_type;
  dictionary_type dictionary(dc);

  graphlab::thread_group threads;
  functor fun(dictionary);
  threads.launch(fun);
  threads.launch(fun);
  threads.join();

  dictionary.flush();


  std::cout << "reached barrier" << std::endl;
  dc.full_barrier();
  


  if(dc.procid() == 0) {
    std::cout << "Finished! " << std::endl
              << "Local Size: " 
              << dictionary.local_size() << std::endl
              << "Dictionary size: " 
              << dictionary.size() << std::endl;
  } 

}


void test_random_gamma(graphlab::distributed_control& dc) {
  graphlab::random::seed(dc.procid());
  word_id_type sum_wid = 0;
  for(size_t i = 0; i < 10000; ++i) {
    const word_id_type word(graphlab::random::gamma(10) * 1000.0);
    sum_wid += word;
  }
  std::cout << "(" << dc.procid() << ", " << sum_wid << ")" << std::endl;
}

int main(int argc, char** argv) {
  std::cout << "Running distributed test" << std::endl;
  
  
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
 
  
  graphlab::random::nondet_seed();
  graphlab::random::seed(dc.procid());
  
  
  

  large_scale_test(dc);
  // test_random_gamma(dc);

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
}


#include <graphlab/macros_undef.hpp>

