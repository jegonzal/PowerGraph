#include <iostream>
#include <graphlab/parallel/fiber.hpp>
using namespace graphlab;
int numticks = 0;
void threadfn(void*) {
  timeval t;
  gettimeofday(&t, NULL);
  double start_time = t.tv_sec + ((double)t.tv_usec)/1.0E6;

  while(1) {
    gettimeofday(&t, NULL);
    double cur_time = t.tv_sec + ((double)t.tv_usec)/1.0E6;
    if (cur_time - start_time >= 1) break;
    fiber_group::yield();
    __sync_fetch_and_add(&numticks, 1);
  }
}

int main(int argc, char** argv) {
  fiber_group fibers(4, 8192);
  timeval t;
  gettimeofday(&t, NULL);
  double start_time = t.tv_sec + ((double)t.tv_usec)/1.0E6;
  for (int i = 0;i < 100000; ++i) {
    fibers.launch(threadfn, NULL);
  }

  fibers.join();
  gettimeofday(&t, NULL);
  double cur_time = t.tv_sec + ((double)t.tv_usec)/1.0E6;
  std::cout << "Completion in " << cur_time - start_time << "s\n";
  std::cout << "Context Switches: " << numticks << "\n";
}
