#include <cxxtest/TestSuite.h>
#include <iostream>
#include <vector>
#include <map>
#include <iterator>


#include <graphlab/macros_def.hpp>

class FTestSuite : public CxxTest::TestSuite {
public:
  void test_foreach(void) {
    // iterate over a vector
    // create the vector
    std::vector<int> int_vector;
    for (int i = 1;i < 10; ++i) {
      int_vector.push_back(i);
    }
    // use foreach to print every element in the vector
    foreach(int& i, int_vector) {
      std::cout << i << "\t";
    }


    std::cout << "\n";

    // create a map
    std::map<int,int> int_map;
    for (int i = 1;i < 10; ++i) {
      int_map[i] = i+1;
    }

    // now, the map's 'value type' is a pair. but the comma in the
    // map definition, messes up the foreach macro expansion.
    typedef std::map<int,int>::value_type int_mapentry;
    foreach(int_mapentry& i, int_map) {
      std::cout << i.first << "->" << i.second << "\n";
    }

    // but the typedef is annoying, so an alternative method,
    // which is gcc specific is to use the typeof command
    // (I believe it might work in MSVC also, but not tested)
    // the former is prefered
    foreach(typeof(*int_map.begin()) &i, int_map) {
      std::cout << i.first << "->" << i.second << "\n";
    }
  }
};
