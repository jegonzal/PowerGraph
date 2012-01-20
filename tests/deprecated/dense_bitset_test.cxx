#include <cxxtest/TestSuite.h>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/macros_def.hpp>
using namespace graphlab;

class DenseBitsetTestSuite : public CxxTest::TestSuite {
public:
  void test_densebitset(void) {
    dense_bitset d;
    d.resize(100);
    uint32_t probelocations[7] = {0, 10, 12, 50, 66, 81, 99};
    // test setting
    for (size_t i= 0;i < 7; ++i) {
      d.set_bit(probelocations[i]);
    }
    
    for (size_t i = 0;i< 100; ++i) {
      bool inprobe=false;
      for (size_t j = 0;j <7; ++j) inprobe |= (probelocations[j] == i);
      TS_ASSERT_EQUALS(d.get(i), inprobe);
    }

    // test iteration
    uint32_t iter = (uint32_t)(-1);
    TS_ASSERT_EQUALS(d.first_bit(iter), true)
    for (size_t i= 0;i < 7; ++i) {
      TS_ASSERT_EQUALS(iter, probelocations[i]);
      bool ret = d.next_bit(iter);
      TS_ASSERT_EQUALS(ret, i < 6);
    }
    size_t ctr = 0;
    foreach(iter, d) {
      TS_ASSERT(ctr < 7);
      TS_ASSERT_EQUALS(iter, probelocations[ctr]);
      ++ctr;
    }

    // testclearing
    for (size_t i= 0;i < 7; ++i) {
      d.clear_bit(probelocations[i]);
    }
    for (size_t i = 0;i< 100; ++i) {
      TS_ASSERT_EQUALS(d.get(i), false);
    }
  }


  void test_fixeddensebitset(void) {
    fixed_dense_bitset<100> d;
    uint32_t probelocations[7] = {0, 10, 12, 50, 66, 81, 99};
    // test setting
    for (size_t i= 0;i < 7; ++i) {
      d.set_bit(probelocations[i]);
    }
    
    for (size_t i = 0;i< 100; ++i) {
      bool inprobe=false;
      for (size_t j = 0;j <7; ++j) inprobe |= (probelocations[j] == i);
      TS_ASSERT_EQUALS(d.get(i), inprobe);
    }

    // test iteration
    uint32_t iter = (uint32_t)(-1);
    TS_ASSERT_EQUALS(d.first_bit(iter), true)
    for (size_t i= 0;i < 7; ++i) {
      TS_ASSERT_EQUALS(iter, probelocations[i]);
      bool ret = d.next_bit(iter);
      TS_ASSERT_EQUALS(ret, i < 6);
    }
    
    size_t ctr = 0;
    foreach(iter, d) {
      TS_ASSERT(ctr < 7);
      TS_ASSERT_EQUALS(iter, probelocations[ctr]);
      ++ctr;
    }

    // testclearing
    for (size_t i= 0;i < 7; ++i) {
      d.clear_bit(probelocations[i]);
    }
    for (size_t i = 0;i< 100; ++i) {
      TS_ASSERT_EQUALS(d.get(i), false);
    }
  }


};

