#include <iostream>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/assertions.hpp>
#include <boost/bind.hpp>

using namespace graphlab;
void test() {
  std::cout << "test" << std::endl;
}

void test_print(size_t s) {
  std::cout << s << std::endl;
}


void thread_assert_false() {
  ASSERT_TRUE(false);
}

int main(int argc, char** argv) {
  // use the helper function
  thread thr1 = launch_in_new_thread(test);
  thr1.join();
  
  // use the built in function. and set to run on cpu 0
  thread thr2;
  thr2.launch(test, 0 /* cpu affinity */ );
  thr2.join();

  
  thread_group thrgroup;
  for (size_t i = 0;i < 10; ++i) {
    // binds the argument (i) to the function test_print
    // so the first thread will call test_print(0)
    // the next thread will call test_print(1), etc
    thrgroup.launch(boost::bind(test_print, i));
  }
  thrgroup.join();
  std::cout << "threads test ok" << std::endl;
  
  thread thr3;
  thr3.launch(thread_assert_false);
  try {
    thr3.join();
  }
  catch(const char* c) {
    std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
  }
  
  
  for (size_t i = 0;i < 10; ++i) {
    thrgroup.launch(thread_assert_false);
  }
  
  size_t numcaught = 0;
  while (thrgroup.running_threads() > 0) {
    try {
      thrgroup.join();
    }
    catch (const char* c){
      std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
      numcaught++;
    }
  }
  std::cout << "Caught " << numcaught << " exceptions!" << std::endl;
  ASSERT_EQ(numcaught, 10);
}
