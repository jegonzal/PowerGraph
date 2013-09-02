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

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/thread_pool.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/timer.hpp>
#include <boost/bind.hpp>

using namespace graphlab;

atomic<int> testval;

void test_inc() {
  usleep(100000);
  testval.inc();
}

void test_dec() {
  usleep(100000);
  testval.dec();
}



void thread_assert_false() {
  ASSERT_TRUE(false);
}


void test_group_exception_forwarding(){
  std::cout << "\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << "This test will print a  large number of assertional failures\n";
  std::cout << "and back traces. This is intentional as we are testing the\n" ;
  std::cout << "exception forwarding scheme\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << std::endl;

  thread_group group;

  
  thread thr3;
  thr3.launch(thread_assert_false);
  try {
    thr3.join();
  }
  catch(const char* c) {
    std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
  }
  
  
  for (size_t i = 0;i < 10; ++i) {
    group.launch(thread_assert_false);
  }
  
  size_t numcaught = 0;
  while (group.running_threads() > 0) {
    try {
      group.join();
    }
    catch (const char* c){
      std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
      numcaught++;
    }
  }
  std::cout << "Caught " << numcaught << " exceptions!" << std::endl;
  TS_ASSERT_EQUALS(numcaught, (size_t)10);
}

void test_pool(){
  testval.value = 0;
  thread_pool pool(4);
  for (size_t j = 0;j < 10; ++j) {
    for (size_t i = 0;i < 10; ++i) {
      pool.launch(test_inc);
    }
    for (size_t i = 0;i < 10; ++i) {
      pool.launch(test_dec);
    }
    pool.set_cpu_affinity(j % 2);
  }
  
  pool.join();
  TS_ASSERT_EQUALS(testval.value, 0);
}

void test_pool_exception_forwarding(){
  std::cout << "\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << "This test will print a  large number of assertional failures\n";
  std::cout << "and back traces. This is intentional as we are testing the\n" ;
  std::cout << "exception forwarding scheme\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << std::endl;
  thread_pool pool(10);

  
  thread thr3;
  thr3.launch(thread_assert_false);
  try {
    thr3.join();
  }
  catch(const char* c) {
    std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
  }
  
  
  for (size_t i = 0;i < 10; ++i) {
    pool.launch(thread_assert_false);
    if (i == 50) {
      pool.set_cpu_affinity(true);
    }
  }
  
  size_t numcaught = 0;
  while (1) {
    try {
      pool.join();
      break;
    }
    catch (const char* c){
      std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
      numcaught++;
    }
  }
  std::cout << "Caught " << numcaught << " exceptions!" << std::endl;
  TS_ASSERT_EQUALS(numcaught, (size_t)10);
}






class ThreadToolsTestSuite : public CxxTest::TestSuite {
public:
  void test_thread_group_exception(void) {
   test_group_exception_forwarding();
  }

  void test_thread_pool(void) {
   test_pool();
  }
   
  void test_thread_pool_exception(void) {
    test_pool_exception_forwarding();
  }

};
