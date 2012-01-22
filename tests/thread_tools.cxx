#include <iostream>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/thread_pool.hpp>
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





template<typename Mutex>
void test_adaptive_mutex_helper(Mutex* mut, size_t* val) {
  for(size_t i = 0; i < 20000000; ++i) {
    mut->lock(); 
    *val += size_t(log(i+1.0) + log(exp(i+1.0) + 1.0) + 1) ; 
    mut->unlock();
  }
}

size_t counter = 0;

void adaptive_mutex_test() {
  const size_t nthreads = 4;
  std::cout << std::endl;
  thread_pool pool(nthreads);
  {
    timer ti;
    ti.start();
    typedef adaptive_mutex<0> mutex_type;
    mutex_type mut;

    for (size_t i = 0; i < nthreads; ++i) {
      pool.launch(boost::bind(test_adaptive_mutex_helper<mutex_type>, &mut, &counter) );
    }
    pool.join();
    std::cout << counter << std::endl;   
    std::cout << ti.current_time() << std::endl;
  }
 
  {
    timer ti;
    ti.start();
    typedef adaptive_mutex<100> mutex_type;
    mutex_type mut;
    for (size_t i = 0; i < nthreads; ++i) {
      pool.launch(boost::bind(test_adaptive_mutex_helper<mutex_type>, &mut, &counter) );
    }
    pool.join();
    std::cout << counter << std::endl;
    std::cout << ti.current_time() << std::endl;
  }

}



class ThreadToolsTestSuite : public CxxTest::TestSuite {
public:
  // void test_thread_group_exception(void) {
  //  test_group_exception_forwarding();
  // }

  // void test_thread_pool(void) {
  //  test_pool();
  // }
   
  // void test_thread_pool_exception(void) {
  //   test_pool_exception_forwarding();
  // }

  void test_adaptive_mutex() {
    adaptive_mutex_test();
  }

};
