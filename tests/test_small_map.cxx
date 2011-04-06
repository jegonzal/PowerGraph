#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/unordered_set.hpp>

#include <cxxtest/TestSuite.h>

#include <graphlab.hpp>
#include <graphlab/util/small_map.hpp>

using namespace graphlab;

#include <graphlab/macros_def.hpp>
class test_small_map : public CxxTest::TestSuite {
public:

  void test_lookup() {
    typedef small_map<10, size_t, std::string> map_type;
    map_type map;
    map[1] = "hello";
    map[0] = "A)";
    map[5] = "!";
    map[2] = " world";
    std::cout << map << std::endl;
    
  }



};
#include <graphlab/macros_undef.hpp>
