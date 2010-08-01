#include <cxxtest/TestSuite.h>
#include <graphlab/util/sampling_tree.hpp>

class SamplingTreeTestSuite: public CxxTest::TestSuite {
  public:
    void test_sampling_tree_uniform(void) {
      graphlab::sampling_tree<size_t, double> stree(50);
      // push 50 uniform elements
      for (size_t i = 0;i < 50;++i) {
        stree.push(i, i, 1.0);
      }
      TS_ASSERT_EQUALS(stree.size(), size_t(50));
      
      
      size_t counts[50] = {0};
      // sample a bunch of times...
      for (size_t i = 0;i < 1000000 ;++i) {
        size_t s=0;
        TS_ASSERT(stree.sample(s));
        counts[s]++;
      }
      for (size_t i = 0;i < 50; ++i) {
        std::cout << i << ": " << counts[i] << "\n";
        counts[i] = 0;
      }
      
      
      // now we pop
      for (size_t i = 0;i < 50; ++i) {
        size_t key=0;
        size_t d=0; 
        double priority=0;
        TS_ASSERT(stree.pop(key, d, priority));
        counts[key]++;
        std::cout << key << "\t";
        TS_ASSERT_EQUALS(counts[key], size_t(1));
      }
      std::cout << "\n";
      TS_ASSERT_EQUALS(stree.size(), size_t(0));
      
      size_t key;
      size_t d; 
      double priority;
      TS_ASSERT(! stree.sample(key));
      TS_ASSERT(! stree.pop(key, d, priority));
    }

    void test_sampling_tree_exponential(void) {
      graphlab::sampling_tree<size_t, double> stree(50);
      // push 50 uniform elements
      for (size_t i = 0;i < 50;++i) {
        stree.push(i, i, std::exp(-double(i)));
      }
      TS_ASSERT_EQUALS(stree.size(), size_t(50));
      
      
      size_t counts[50] = {0};
      // sample a bunch of times...
      for (size_t i = 0;i < 1000000 ;++i) {
        size_t s = 0;
        TS_ASSERT(stree.sample(s));
        counts[s]++;
      }
      for (size_t i = 0;i < 50; ++i) {
        std::cout << i << ": " << counts[i] << "\n";
        counts[i] = 0;
      }
      
      
      // now we pop
      for (size_t i = 0;i < 50; ++i) {
        size_t key = 0;
        size_t d =0; 
        double priority=0;
        TS_ASSERT(stree.pop(key, d, priority));
        counts[key]++;
        std::cout << key << "\t";
        TS_ASSERT_EQUALS(counts[key], size_t(1));
      }
      std::cout << "\n";
      TS_ASSERT_EQUALS(stree.size(), size_t(0));
      
      size_t key;
      size_t d; 
      double priority;
      stree.print_tree();
      TS_ASSERT(! stree.sample(key));
      TS_ASSERT(! stree.pop(key, d, priority));
    }
    
    void test_sampling_tree_order(void) {
        graphlab::sampling_tree<size_t, double> stree(50);
      	// push 50 uniform elements
      	for (size_t i = 0;i < 50;++i) {
        	stree.push(i, i, std::exp(-double(i)));
      	}
    	 for (size_t i = 0;i < 50;++i) {
    	     size_t key = 0;
        	size_t d =0; 
        	double priority=0;
      		stree.pop(key, d, priority);
      		std::cout << "Popped " << key << ": " << priority << "\n";
      	}
    }
    
};
