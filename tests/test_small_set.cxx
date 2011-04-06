#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/unordered_set.hpp>

#include <cxxtest/TestSuite.h>

#include <graphlab.hpp>
#include <graphlab/util/small_set.hpp>

using namespace graphlab;

#include <graphlab/macros_def.hpp>
class test_small_set : public CxxTest::TestSuite {
public:

  
  void test_union() {
    std::cout << std::endl;
    std::cout << "Testing set union" << std::endl;

    typedef small_set<10, int> set_type;
    typedef small_set<5, int> small_set_type;
    small_set<0, int> empty_set;
    small_set<10, int> set1;
    small_set<10, int> set2;
    set1 += set_type(1) + small_set_type(3) + set_type(2) + empty_set;
    set1 += 1;
    set1 += 3;
    set1 += 2;
    set1 += empty_set;
    set1 += 1;
    std::set<int> true_set1;
    true_set1.insert(1);
    true_set1.insert(2);
    true_set1.insert(3);
    ASSERT_EQ(set_type(true_set1), set1);
    std::cout << "set1: " << set1 << std::endl;
    set2 += set_type(2) + small_set_type(5) + small_set_type(3) + set_type(7);
    set2.insert(0);
    set2 += 7;
    set2 += 0;

    std::set<int> true_set2;
    true_set2.insert(0);
    true_set2.insert(2);
    true_set2.insert(5);
    true_set2.insert(3);
    true_set2.insert(7);
    ASSERT_EQ(set_type(true_set2), set2);    
    std::cout << "set2: " << set2 << std::endl;

    small_set<7, int> set3 = set1 + set2;
    std::set<int> true_set3 = set_union(true_set1, true_set2);
    ASSERT_EQ(set_type(true_set3), set3);
    std::cout << "set3 = set1 + set2: " << set3 << std::endl;
    std::cout << "set3 + set3: " << (set3  + set3) << std::endl;
    ASSERT_EQ(set_type(true_set3), (set3 + set3));    
  }


  void test_union_speed() {
    typedef small_set<20, int> set_type;
    typedef std::set<int> true_set_type;
    typedef boost::unordered_set<int> boost_set_type;
    size_t max_iter = 1000000;
    true_set_type true_set1;
    set_type set1;
    boost_set_type boost_set1;
    for(size_t i = 0; i < 15; ++i) {
      true_set1.insert(i);
      set1.insert(i);
      boost_set1.insert(i);
    }
    ASSERT_EQ(set_type(true_set1), set1);
    true_set_type true_set2;
    set_type set2;
    boost_set_type boost_set2;
    for(size_t i = 5; i < 20; ++i) {
      true_set2.insert(i);
      set2.insert(i);
      boost_set2.insert(i);
    }
    ASSERT_EQ(set_type(true_set2), set2);
    

    timer time;



    set_type set3;
    time.start();
    for(size_t i = 0; i < max_iter; ++i) {
      set3 = set1 + set2;
    }
    std::cout << "small set time: " << time.current_time() << std::endl;


    boost_set_type boost_set3;
    time.start();
    for(size_t i = 0; i < max_iter; ++i) {
      boost_set3.clear();
      boost_set3.insert(boost_set1.begin(), boost_set1.end());
      boost_set3.insert(boost_set2.begin(), boost_set2.end());
    }
    std::cout << "boost set time: " << time.current_time() << std::endl;


    true_set_type true_set3;
    time.start();
    for(size_t i = 0; i < max_iter; ++i) {
      true_set3 = set_union(true_set1, true_set2);
    }
    std::cout << "Std set time: " << time.current_time() << std::endl;

    ASSERT_EQ(set_type(true_set3), set3);


  }
  

};

#include <graphlab/macros_undef.hpp>
