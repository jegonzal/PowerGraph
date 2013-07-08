#include <iostream>
#include <graphlab/parallel/fiber_group.hpp>
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


void threadfn2() {

  timer ti; ti.start();
  while(1) {
    if (ti.current_time() >= 2) break;
    fiber_control::yield();
    __sync_fetch_and_add(&numticks, 2);
  }
}

int main(int argc, char** argv) {
  timer ti; ti.start();
  fiber_group group;
  fiber_group group2;
  for (int i = 0;i < 100000; ++i) {
    group.launch(threadfn);
    group2.launch(threadfn2);
  }
  group.join();
  std::cout << "Completion in " << ti.current_time() << "s\n";
  std::cout << "Context Switches: " << numticks << "\n";
  group2.join();
  std::cout << "Completion in " << ti.current_time() << "s\n";
  std::cout << "Context Switches: " << numticks << "\n";
}
