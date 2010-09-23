#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>

#include <cxxtest/TestSuite.h>

#include <graphlab/serialization/serialize.hpp>
#include <graphlab/serialization/vector.hpp>
#include <graphlab/serialization/map.hpp>
#include <graphlab/serialization/list.hpp>
#include <graphlab/serialization/set.hpp>


using namespace graphlab;



struct A{
  int z;
  void save(oarchive &a) const {
    a << z;
  }
  void load(iarchive a) {
    a >> z;
  }
};

class TestClass{
public:
  int i;
  int j;
  std::vector<int> k;
  A l;
  void save(oarchive &a) const {
    a << i << j << k << l;
  }
  void load(iarchive &a) {
    a >> i >> j >> k >> l;
  }
};


class SerializeTestSuite : public CxxTest::TestSuite {
  public:

    // Look for the class TestClass() to see the most interesting tutorial on how to
    // use the serializer
    void test_basic_datatype(void) {
      char t1 = 'z';
      bool t2 = true;
      int t3 = 10;
      int t4 = 18345;
      long t5 = 30921233;
      long long t6 = (long long)(t5)*100;
      float t7 = 10.35;
      double t8 = 3.14156;
      const char *t9 = "hello world";
      const char * t10 = "blue";

      char r1;
      bool r2;
      int r3;
      int r4;
      long r5;
      long long r6;
      float r7;
      double r8;
      char r9[100];
      char r10[10];

      // serialize t1-10
      std::ofstream f;
      f.open("test.bin",std::fstream::binary);
      oarchive a(f);
      a << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8;
      serialize(a, t9, strlen(t9) + 1);
      serialize(a, t10, strlen(t10) + 1);
      f.close();

      // deserialize into r1-10
      std::ifstream g;
      g.open("test.bin",std::fstream::binary);
      iarchive b(g);
      b >> r1 >> r2 >> r3 >> r4 >> r5 >> r6 >> r7 >> r8;
      deserialize(b, &r9, strlen(t9) + 1);
      deserialize(b, r10, strlen(t10) + 1);
      g.close();

      TS_ASSERT_EQUALS(t1, r1);
      TS_ASSERT_EQUALS(t2, r2);
      TS_ASSERT_EQUALS(t3, r3);
      TS_ASSERT_EQUALS(t4, r4);
      TS_ASSERT_EQUALS(t5, r5);
      TS_ASSERT_EQUALS(t6, r6);
      TS_ASSERT_EQUALS(t7, r7);
      TS_ASSERT_EQUALS(t8, r8);
      TS_ASSERT_SAME_DATA(t9, r9, strlen(t9) + 1);
      TS_ASSERT_SAME_DATA(t10, r10, strlen(t10) + 1);
    }

    void test_vector_serialization(void) {
      std::vector<int> v;
      for (int i = 0;i< 10; ++i) {
        v.push_back(i);
      }
      std::ofstream f;
      f.open("test.bin",std::fstream::binary);
      oarchive a(f);
      a << v;
      f.close();

      std::vector<int> w;
      std::ifstream g;
      iarchive b(g);
      g.open("test.bin",std::fstream::binary);
      b >> w;
      g.close();

      for (int i = 0;i< 10; ++i) {
        TS_ASSERT_EQUALS(v[i], w[i]);
      }
    }


    void test_class_serialization(void) {
      // create a test class
      TestClass t;
      t.i=10;
      t.j=20;
      t.k.push_back(30);

      //serialize
      std::ofstream f;
      f.open("test.bin",std::fstream::binary);
      oarchive a(f);
      a << t;
      f.close();
      //deserialize into t2
      TestClass t2;
      std::ifstream g;
      g.open("test.bin",std::fstream::binary);
      iarchive b(g);
      b >> t2;
      g.close();
      // check
      TS_ASSERT_EQUALS(t.i, t2.i);
      TS_ASSERT_EQUALS(t.j, t2.j);
      TS_ASSERT_EQUALS(t.k.size(), t2.k.size());
      TS_ASSERT_EQUALS(t.k[0], t2.k[0]);
    }

    void test_vector_of_classes(void) {
      // create a vector of test classes
      std::vector<TestClass> vt;
      vt.resize(10);
      for (int i=0;i<10;i++) {
        vt[i].i=i;
        vt[i].j=i*21;
        vt[i].k.resize(10);
        vt[i].k[i]=i*51;
      }

      //serialize
      std::ofstream f;
      f.open("test.bin",std::fstream::binary);
      oarchive a(f);
      a << vt;
      f.close();

      //deserialize into vt2
      std::vector<TestClass> vt2;
      std::ifstream g;
      g.open("test.bin",std::fstream::binary);
      iarchive b(g);
      b >> vt2;
      g.close();
      // check
      TS_ASSERT_EQUALS(vt.size(), vt2.size());
      for (size_t i=0;i<10;i++) {
        TS_ASSERT_EQUALS(vt[i].i, vt2[i].i);
        TS_ASSERT_EQUALS(vt[i].j, vt2[i].j);
        TS_ASSERT_EQUALS(vt[i].k.size(), vt2[i].k.size());
        for (size_t j = 0; j < vt[i].k.size(); ++j) {
          TS_ASSERT_EQUALS(vt[i].k[j], vt2[i].k[j]);
        }
      }
    }

    void test_vector_of_strings(void) {
      std::string x = "Hello world";
      std::string y = "This is a test";
      std::vector<std::string> v;
      v.push_back(x); v.push_back(y);

      std::ofstream f;
      f.open("test.bin",std::fstream::binary);
      oarchive a(f);
      a << v;
      f.close();

      //deserialize into vt2
      std::vector<std::string> v2;
      std::ifstream g;
      g.open("test.bin",std::fstream::binary);
      iarchive b(g);
      b >> v2;
      g.close();
      TS_ASSERT_EQUALS(v[0], v2[0]);
      TS_ASSERT_EQUALS(v[1], v2[1]);
    }

    void test_map_serialization(void) {
      std::map<std::string,int> v;
      v["one"] = 1;
      v["two"] = 2;
      v["three"] = 3;

      std::ofstream f;
      f.open("test.bin",std::fstream::binary);
      oarchive a(f);
      a << v;
      f.close();

      //deserialize into vt2
      std::map<std::string,int> v2;
      std::ifstream g;
      g.open("test.bin",std::fstream::binary);
      iarchive b(g);
      b >> v2;
      g.close();
      TS_ASSERT_EQUALS(v["one"], v2["one"]);
      TS_ASSERT_EQUALS(v["two"], v2["two"]);
      TS_ASSERT_EQUALS(v["three"], v2["three"]);
    }

};

