#include <iostream>
#include <graphlab/parallel/pthread_tools.hpp>
#include <boost/bind.hpp>

using namespace graphlab;
void test() {
  std::cout << "test" << std::endl;
}

void test_print(size_t s) {
  std::cout << s << std::endl;
}

int main(int argc, char** argv) {
  thread_group thrgroup;
  thread thr1 = launch_in_new_thread(test);
  thread thr2 = launch_in_new_thread(test);
  thr1.join();
  thr2.join();
  
  for (size_t i = 0;i < 10; ++i) {
    launch_in_new_thread(thrgroup, boost::bind(test_print, i));
  }
  thrgroup.join();
}
