#include <iostream>
#include <graphlab/parallel/fiber_control.hpp>
#include <graphlab/util/timer.hpp>
using namespace graphlab;
int numticks = 0;
void threadfn() {

  timer ti; ti.start();
  while(1) {
    if (ti.current_time() >= 1) break;
    fiber_control::yield();
    __sync_fetch_and_add(&numticks, 1);
  }
}

int main(int argc, char** argv) {
  timer ti; ti.start();
  for (int i = 0;i < 100000; ++i) {
    fiber_control::get_instance().launch(threadfn);
  }
  fiber_control::get_instance().join();
  std::cout << "Completion in " << ti.current_time() << "s\n";
  std::cout << "Context Switches: " << numticks << "\n";
}
