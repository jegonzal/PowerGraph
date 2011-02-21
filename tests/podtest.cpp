#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>

#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/serialize.hpp>
#include <graphlab/serialization/vector.hpp>
#include <graphlab/serialization/map.hpp>
#include <graphlab/serialization/list.hpp>
#include <graphlab/serialization/set.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/serialization/serializable_pod.hpp>
#include <graphlab/serialization/podify.hpp>

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

void is10(const graphlab::any a) {
  ASSERT_EQ(a.as<int>(), 10);
}

void testclass3test(TestClass3 &tc) {
  ASSERT_EQ(tc.z, 100);
}
int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  std::ofstream f;
  f.open("test.bin",std::fstream::binary);
  oarchive oarc(f);
  graphlab::any variant;

  // store a bunch of stuff
  int i = 10;
  variant = i;
  ASSERT_EQ(variant.as<int>(), 10);
  is10(variant);
  oarc << variant;
  
  
  double d = 3.14159;
  variant = d;
  ASSERT_GT(variant.as<double>(), 3.14158);
  ASSERT_LE(variant.as<double>(), 3.1416);
  oarc << variant;

  TestClass1 t;
  t.z = 4321;
  variant = t;
  ASSERT_EQ(variant.as<TestClass1>().z, 4321);
  oarc << variant;

  TestClass2 t2;
  t2.i = 1;
  t2.j = 2;
  t2.l.z = 3;
  for (size_t i = 0;i < 10; ++i) {
    t2.k.push_back(i);
  }
  variant = t2;
  oarc << variant;
  
  TestClass3 t3;
  t3.z = 100;
  oarc << PODIFY(t3);
  testclass3test(PODIFY(t3));

  const char* c = "hello world";
  oarc << c;
  f.close();
}
