#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>

#include <graphlab/logger/logger.hpp>

#include <graphlab/serialization/serialize.hpp>
#include <graphlab/serialization/vector.hpp>
#include <graphlab/serialization/map.hpp>
#include <graphlab/serialization/list.hpp>
#include <graphlab/serialization/set.hpp>
#include <graphlab/serialization/serializable_pod.hpp>
#include <graphlab/serialization/podify.hpp>

#include <graphlab/util/generics/any.hpp>

using namespace graphlab;



struct TestClass1{
  int z;
};
SERIALIZABLE_POD(TestClass1);

class TestClass2{
public:
  int i;
  int j;
  std::vector<int> k;
  TestClass1 l;
  void save(oarchive &a) const {
    a << i << j << k << l;
  }
  void load(iarchive &a) {
    a >> i >> j >> k >> l;
  }
};


struct TestClass3{
  int z;
};


void testclass3test(TestClass3 tc) {
  ASSERT_EQ(tc.z, 100);
}

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  std::ifstream f;
  f.open("test.bin",std::fstream::binary);
  iarchive iarc(f);
  graphlab::any variant;

  // store a bunch of stuff
  iarc >> variant;
  ASSERT_EQ(variant.as<int>(), 10);

  iarc >> variant;
  ASSERT_GT(variant.as<double>(), 3.14158);
  ASSERT_LE(variant.as<double>(), 3.1416);

  iarc >> variant;
  ASSERT_EQ(variant.as<TestClass1>().z, 4321);


  iarc >> variant;
  TestClass2 tc = variant.as<TestClass2>();
  ASSERT_EQ(tc.i, 1);
  ASSERT_EQ(tc.j, 2);
  ASSERT_EQ(tc.l.z, 3);
  ASSERT_EQ(tc.k.size(), 10);
  for (size_t i = 0;i < 10; ++i) {
    ASSERT_EQ(tc.k[i], i);
  }
  
  TestClass3 t3;
  t3.z = 100;
  iarc >> PODIFY(t3);
  
  char *c = NULL;
  iarc >> c;
  std::cout << c;
  testclass3test(PODIFY(t3));
}
