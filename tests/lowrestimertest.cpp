#include <graphlab/util/timer.hpp>

int main(int argc, char** argv) {
  graphlab::timer ti;
  ti.start();
  for (size_t i = 0;i < 10; ++i) {
    std::cout << "timer: " << ti << std::endl;
    std::cout << "lowres timer: " << graphlab::lowres_time_millis() << std::endl;
    getchar();
  }
}