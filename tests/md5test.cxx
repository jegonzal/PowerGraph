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


#include <graphlab/util/md5.hpp>
#include <fstream>

void hash_verify(std::string teststring, std::string realhash) {
  std::stringstream strm;
  strm << graphlab::MD5(teststring);
  TS_ASSERT_EQUALS(realhash, strm.str());
}

class RandomTestSuite: public CxxTest::TestSuite {
 public:  
  void test_md5_verification() {
   // test hashes as according to RFC 1321
   hash_verify("", "d41d8cd98f00b204e9800998ecf8427e");
   hash_verify("a", "0cc175b9c0f1b6a831c399e269772661");
   hash_verify("abc", "900150983cd24fb0d6963f7d28e17f72");
   hash_verify("message digest", 
               "f96b697d7cb7938d525a2f31aaf161d0");
   hash_verify("abcdefghijklmnopqrstuvwxyz", 
               "c3fcd3d76192e4007dfb496cca67e13b");
   hash_verify("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
               "d174ab98d277d9f5a5611c2c9f419d9f");
   hash_verify("12345678901234567890123456789012345678901234567890123456789012345678901234567890",
               "57edf4a22be3c955ac49da2e2107b67a");
  }
};