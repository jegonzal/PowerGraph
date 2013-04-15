/*  
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
#include <iostream>
#include <cxxtest/TestSuite.h>

#include <graphlab/util/generics/csr_storage.hpp>
#include <graphlab/util/generics/dynamic_csr_storage.hpp>
#include <graphlab/util/generics/shuffle.hpp>
#include <graphlab/logger/assertions.hpp>

class csr_storage_test : public CxxTest::TestSuite {  
 public:
  typedef int valuetype;
  typedef int keytype;
  typedef size_t sizetype;

  typedef graphlab::csr_storage<valuetype, sizetype> csr_storage;
  typedef graphlab::dynamic_csr_storage<valuetype, sizetype, 4> dynamic_csr_storage;
  
 public:
  csr_storage_test() {
    keytype keyin_arr[] = {1, 3, 6, 9, 5, 2};
    valuetype valin_arr[] = {3, 2, 1, 4, 4, 4};

    _keyin.assign(keyin_arr, keyin_arr + sizeof(keyin_arr) / sizeof(keytype));
    _valin.assign(valin_arr, valin_arr + sizeof(valin_arr) / sizeof(valuetype));

    keytype keyout_arr[] = {1, 2, 3, 5, 6, 9};
    valuetype valout_arr[] = {3, 4, 2, 4, 1, 4};

    _keyout.assign(keyout_arr, keyout_arr + sizeof(keyout_arr) / sizeof(keytype));
    _valout.assign(valout_arr, valout_arr + sizeof(valout_arr) / sizeof(valuetype));
  }

  void test_csr_storage() {
    std::cout << "Test csr_storage constructor" << std::endl;
    csr_storage csr(get_keyin(), get_valin());
    check(csr, get_keyout(), get_valout());
    csr.print(std::cout);
    printf("+ Pass test: csr_storage constructor :)\n\n");
  }

  void test_csr_storage2() {
    std::cout << "Test csr_storage wrap " << std::endl;
    std::vector<keytype> keys(get_keyin());
    std::vector<valuetype> values(get_valin());

    std::vector<sizetype> permute_index;
    std::vector<sizetype> prefix;

    graphlab::counting_sort(keys, permute_index, &prefix);
    graphlab::outofplace_shuffle(values, permute_index);

    csr_storage csr;
    csr.wrap(prefix, values);
    check(csr, get_keyout(), get_valout());
    printf("+ Pass test: csr_storage wrap :)\n\n");
  }

  void test_dynamic_csr_storage() {
    std::cout << "Test dynamic csr_storage constructor" << std::endl;
    dynamic_csr_storage csr(get_keyin(), get_valin());
    check(csr, get_keyout(), get_valout());
    printf("+ Pass test: dynamic_csr_storage constructor :)\n\n");
  }
  
  void test_dynamic_csr_storage2() {
    std::cout << "Test dynamic csr_storage wrap" << std::endl;
    std::vector<keytype> keys(get_keyin());
    std::vector<valuetype> values(get_valin());

    std::vector<sizetype> permute_index;
    std::vector<sizetype> prefix;

    graphlab::counting_sort(keys, permute_index, &prefix);
    graphlab::outofplace_shuffle(values, permute_index);
    dynamic_csr_storage csr;
    csr.wrap(prefix, values);
    check(csr, get_keyout(), get_valout());
    printf("+ Pass test: dynamic_csr_storage wrap:)\n\n");
  }

  void test_dynamic_csr_storage_insertion() {
    std::vector<keytype> keys(get_keyin());
    std::vector<valuetype> values(get_valin());

    dynamic_csr_storage csr;
    for (size_t i = 0; i < keys.size(); ++i) {
      csr.insert(keys[i], values[i]);
    }
    check(csr, get_keyout(), get_valout());
  }

 private:
  template<typename csr_type>
      void check(csr_type& csr,
                 std::vector<keytype> keyout,
                 std::vector<valuetype> valout) {
        typedef typename csr_type::iterator iterator;
        size_t id = 0;
        for (size_t i = 0; i < csr.num_keys(); ++i) {
          iterator iter = csr.begin(i);
          while (iter != csr.end(i)) {
            ASSERT_EQ(i, keyout[id]);
            ASSERT_EQ(*iter, valout[id]); 
            ++iter;
            ++id;
          }
        }
      }

  std::vector<keytype> get_keyin() { return std::vector<keytype>(_keyin); }

  std::vector<valuetype> get_valin() { return std::vector<valuetype>(_valin); }
  
  std::vector<keytype> get_keyout() { return std::vector<keytype>(_keyout); }

  std::vector<valuetype> get_valout() { return std::vector<valuetype>(_valout); }

  std::vector<keytype> _keyin;
  std::vector<keytype> _keyout;
  std::vector<valuetype> _valin;
  std::vector<valuetype> _valout;
}; // end of test
