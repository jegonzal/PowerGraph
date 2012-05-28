#include <vector>
#include <algorithm>
#include <iostream>

#include <graphlab/util/empty.hpp>

#include <cxxtest/TestSuite.h>

#include <graphlab/macros_def.hpp>

using namespace graphlab;

class empty_test : public CxxTest::TestSuite {
public:

  
  void test_empty() {
    std::vector<empty> v;
    v.resize(100);
    TS_ASSERT_EQUALS(v.size(), (size_t)100);
    size_t c = 0;
    foreach(empty e, v) {
      ++c;
    }
    TS_ASSERT_EQUALS(c, (size_t)100);

    TS_ASSERT_EQUALS(v.end() - v.begin(), (int)100);
    TS_ASSERT_EQUALS(v.rend() - v.rbegin(), (int)100);
    
    std::vector<empty> v2;
    v2.assign(v.begin(), v.end());
    TS_ASSERT_EQUALS(v2.size(), (size_t)100);
    v2.assign(v.rbegin(), v.rend());
    TS_ASSERT_EQUALS(v2.size(), (size_t)100);

    v.insert(v.begin(), empty());
    TS_ASSERT_EQUALS(v.size(), (size_t)101);
     v.insert(v.end(), empty());
    TS_ASSERT_EQUALS(v.size(), (size_t)102);

    std::vector<empty>::const_iterator iter = v.begin();
    (*iter);
    ++iter; TS_ASSERT_EQUALS(v.end() - iter, (int)101);
    --iter; TS_ASSERT_EQUALS(v.end() - iter, (int)102);
    iter+=10; TS_ASSERT_EQUALS(v.end() - iter, (int)92);
    iter-=10; TS_ASSERT_EQUALS(v.end() - iter, (int)102);
    std::vector<empty>::const_iterator iter2 = iter;
    iter2 += 10;
    TS_ASSERT_EQUALS(iter2 - iter, (int)10);
  }

};

