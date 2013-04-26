#include <iostream>
#include <graphlab/parallel/fiber.hpp>
#include <graphlab/util/timer.hpp>
using namespace graphlab;
int numticks = 0;
void threadfn() {

  timer ti; ti.start();
  while(1) {
    if (ti.current_time() >= 1) break;
    fiber_group::yield();
    __sync_fetch_and_add(&numticks, 1);
  }
}

int main(int argc, char** argv) {
  fiber_group fibers(4, 8192);
  timer ti; ti.start();
  for (int i = 0;i < 100000; ++i) {
    fibers.launch(threadfn);
  }
  fibers.join();
  std::cout << "Completion in " << ti.current_time() << "s\n";
  std::cout << "Context Switches: " << numticks << "\n";
}
