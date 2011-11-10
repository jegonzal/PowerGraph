
#include <iostream>
#include <vector>
#include <algorithm>
#include <cxxtest/TestSuite.h>
#include <graphlab/util/fast_multinomial.hpp>
#include <graphlab/parallel/pthread_tools.hpp>


#include <graphlab/macros_def.hpp>


class FastMultinomialTestSuite: public CxxTest::TestSuite {
  public:
  void test_fast_multinomial_uniform(void) {
    for(size_t num_asg = 1; num_asg < 50; ++num_asg) {
      std::cout << "=================================================="
                << std::endl
                << "Num asg: " << num_asg << std::endl;
      graphlab::fast_multinomial multi(num_asg, 1);
      //         multi.print_tree();
      // push 50 uniform elements
      for (size_t i = 0; i < num_asg; ++i) {
        multi.set(i, 1.0);
      }
      //        multi.print_tree();
      std::vector<double> count(num_asg, 0);
      size_t total = 0;
      size_t nsamples = 1000;
      std::cout << "Sampling" << std::endl;
      for(size_t i = 0; i < nsamples; ++i) {
        size_t value = -1;
        if(multi.sample(value, 0)) {
          count.at(value)++;
          total++;
        } else {
          std::cout << "FAIL!!" << std::endl;
           assert(false);
        }
      }
      // normalize the counts
      for(size_t i = 0; i < count.size(); ++i) {
        count[i] /= total;
      }        
      double max = *std::max_element(count.begin(), count.end());
      double min = *std::min_element(count.begin(), count.end());
      std::cout << "max : min -- " << max << " : " << min << std::endl;
      
      std::cout << "====================================="
                << std::endl
                << "Pop Test" << std::endl;
      
      // Try again using a popping strategy
      std::fill(count.begin(), count.end(), 0);
      total = 0;
      nsamples = 10000;
      for(size_t i = 0; i < nsamples; ++i) {
        size_t value = -1;
        if(multi.pop(value, 0)) {
          assert(value < num_asg);
          count.at(value)++;
          total++;
          multi.set(value, 1.0);
        } else {
          std::cout << "FAIL!!" << std::endl;
          assert(false);
        }
      }
      // normalize the counts
      for(size_t i = 0; i < count.size(); ++i) {
        count[i] /= total;
      }
      max = *std::max_element(count.begin(), count.end());
      min = *std::min_element(count.begin(), count.end());
      std::cout << "max : min -- " << max << " : " << min << std::endl;
      
    }
    
  }

  void test_fast_multinomial_nonuniform(void) {
    for(size_t num_asg = 1; num_asg < 50; ++num_asg) {
      std::cout << "=================================================="
                << std::endl
                << "Num asg: " << num_asg << std::endl;
      graphlab::fast_multinomial multi(num_asg,1);
      //         multi.print_tree();
      // push 50 uniform elements
      for (size_t i = 0; i < num_asg; ++i) {
        multi.set(i, double(i + 1));
      }
      //        multi.print_tree();
      std::vector<double> count(num_asg, 0);
      size_t total = 0;
      size_t nsamples = 10000;
      std::cout << "Sampling" << std::endl;
      for(size_t i = 0; i < nsamples; ++i) {
        size_t value = -1;
        if(multi.sample(value, 0)) {
          count.at(value)++;
          total++;
        } else {
          std::cout << "FAIL!!" << std::endl;
           assert(false);
        }
      }
      // normalize the counts
      for(size_t i = 0; i < count.size(); ++i) {
        count[i] /= total;
        std::cout << count[i] << "  ";
      }
      std::cout << std::endl;
      double max = *std::max_element(count.begin(), count.end());
      double min = *std::min_element(count.begin(), count.end());
      std::cout << "max : min -- " << max << " : " << min << std::endl;      




      std::cout << "====================================="
                << std::endl
                << "Pop Test" << std::endl;
      
      // Try again using a popping strategy
      std::fill(count.begin(), count.end(), 0);
      total = 0;
      nsamples = 10000;
      for(size_t i = 0; i < nsamples; ++i) {
        size_t value = -1;
        if(multi.pop(value, 0)) {
          assert(value < num_asg);
          count.at(value)++;
          total++;
          multi.set(value, double(value + 1.0));
        } else {
          std::cout << "FAIL!!" << std::endl;
          assert(false);
        }
      }
      // normalize the counts
      for(size_t i = 0; i < count.size(); ++i) {
        count[i] /= total;
      }
      max = *std::max_element(count.begin(), count.end());
      min = *std::min_element(count.begin(), count.end());
      std::cout << "max : min -- " << max << " : " << min << std::endl;
      

    }
  }


  class worker {
  public:
    size_t id;
    size_t nsamples;
    size_t num_asg;
    std::vector<double> count;
    graphlab::fast_multinomial* multi;
  
    void run() {
      for(size_t i = 0; i < nsamples; ++i) {
        size_t value = -1;
        if(multi->pop(value, 0)) {
          assert(value < num_asg);
          count.at(value)++;
          multi->set(value, 1.0);
        } else {
          std::cout << "FAIL!!" << std::endl;
          multi->print_tree();
          assert(false);
        }
      }
    }

  }; // end of worker

  
  class biased_worker  {
  public:
    size_t id;
    size_t nsamples;
    size_t num_asg;
    std::vector<double> count;
    graphlab::fast_multinomial* multi;
  
    void run() {
      for(size_t i = 0; i < nsamples; ++i) {
        size_t value = -1;
        if(multi->pop(value, 0)) {
          assert(value < num_asg);
          count.at(value)++;
          multi->set(value, double(value + 1.0));
        } else {
          std::cout << "FAIL!!" << std::endl;
          multi->print_tree();
          assert(false);
        }
      }
    }

  }; // end of worker

  
  void test_fast_multinomial_threaded(void) {
    for(size_t num_asg = 1; num_asg < 50; ++num_asg) {
      std::cout << "=================================================="
                << std::endl
                << "Num asg: " << num_asg << std::endl;
     size_t num_threads = num_asg < 5 ? 1 : num_asg - 3; 
     graphlab::fast_multinomial multi(num_asg, num_threads);
      for (size_t i = 0; i < num_asg; ++i) {
        multi.set(i, 1.0);
      }

      size_t nsamples = 10000;
 

      std::cout << "Threads: " << num_threads << std::endl;
      std::vector<worker> workers(num_threads);
      graphlab::thread_group threads;
      for(size_t i = 0; i < workers.size(); ++i) {
        workers[i].nsamples = nsamples;
        workers[i].num_asg = num_asg;
        workers[i].count.resize(num_asg, 0);
        workers[i].multi = &multi;
        workers[i].id = i;
        threads.launch(boost::bind(&worker::run, &workers[i]));
      }
      threads.join();

      // Add up all the threads
      std::vector<double> count(num_asg, 0);
      for(size_t i = 0; i < workers.size(); ++i) {
        for(size_t asg = 0; asg < num_asg; ++asg) {
          count[asg] += workers[i].count[asg];
        }
      }
      
      // normalize the counts
      double total = 0;
      for(size_t i = 0; i < count.size(); ++i) total += count[i];
      for(size_t i = 0; i < count.size(); ++i) count[i] /= total;
      double max = *std::max_element(count.begin(), count.end());
      double min = *std::min_element(count.begin(), count.end());
      std::cout << "max : min -- " << max << " : " << min << std::endl;
      
    }
    
  }

  void test_fast_multinomial_threaded_biased(void) {
    for(size_t num_asg = 1; num_asg < 50; ++num_asg) {
      std::cout << "=================================================="
                << std::endl
                << "Num asg: " << num_asg << std::endl;
      size_t num_threads = num_asg < 5 ? 1 : num_asg - 3;
      graphlab::fast_multinomial multi(num_asg, num_threads);
      for (size_t i = 0; i < num_asg; ++i) {
        multi.set(i, 1.0);
      }

      size_t nsamples = 10000;


      std::cout << "Threads: " << num_threads << std::endl;
      std::vector<biased_worker> workers(num_threads);
      graphlab::thread_group threads;
      for(size_t i = 0; i < workers.size(); ++i) {
        workers[i].nsamples = nsamples;
        workers[i].num_asg = num_asg;
        workers[i].count.resize(num_asg, 0);
        workers[i].multi = &multi;
        workers[i].id = i;
        threads.launch(boost::bind(&biased_worker::run, &workers[i]));
      }
      threads.join();

      // Add up all the threads
      std::vector<double> count(num_asg, 0);
      for(size_t i = 0; i < workers.size(); ++i) {
        for(size_t asg = 0; asg < num_asg; ++asg) {
          count[asg] += workers[i].count[asg];
        }
      }
      
      // normalize the counts
      double total = 0;
      for(size_t i = 0; i < count.size(); ++i) total += count[i];
      for(size_t i = 0; i < count.size(); ++i) count[i] /= total;
      double max = *std::max_element(count.begin(), count.end());
      double min = *std::min_element(count.begin(), count.end());
      std::cout << "max : min -- " << max << " : " << min << std::endl;
      
    }
    
  }


    
};
