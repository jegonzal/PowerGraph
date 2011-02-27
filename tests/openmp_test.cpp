
#include <omp.h>
#include <cmath>
#include <iostream>


template <typename T>
struct S {
  T pika;
  void test(const T& t) {
#pragma omp parallel for
    for(size_t i = 0;i < 100; ++i) {
      pika = 10;
      for (size_t j = 0; j < 100000000; ++j) {
        volatile double g = 1000000.0;
        volatile double t = sqrt(g);
        g = g + t;
      } 
      std::cout << i << ": " << omp_get_thread_num() << std::endl;
    }
  }
};

int main(int argc, char** argv) {
  S<int> a;
  S<double> b;
  a.test(1);
  b.test(1.0);
}


