/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


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
#include <graphlab/util/generics/any.hpp>

using namespace graphlab;



struct TestClass1{
  int z;
  void save(oarchive &a) const {
    a << z;
  }
  void load(iarchive &a) {
    a >> z;
  }
};

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
}

