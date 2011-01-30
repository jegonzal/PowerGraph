#include <cstring>
#include <iostream>
#include <vector>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/mmap_wrapper.hpp>

using namespace graphlab;


const size_t passes = 10;
const size_t msize = 100 * 1024 * 1024;

template <size_t WSIZE>
double random_writetest(char* c) {
  char* v = new char[WSIZE];
  for (size_t i = 0;i < WSIZE; ++i) {
    v[i] = rand() % 256;
  }
  std::vector<size_t> randpos;
  // how many?
  size_t numrand = passes * msize / WSIZE;
  randpos.resize(numrand);
  for (size_t i = 0;i < numrand; ++i) {
    randpos[i] = graphlab::random::rand_int(msize - WSIZE - 1);
  }
  
  
  timer ti;
  ti.start();
  for (size_t i = 0; i < randpos.size(); ++i) {
        memcpy(c + randpos[i], v, WSIZE);
  }
  
  randpos.clear();
  delete [] v;
  return ti.current_time();
}


template <size_t WSIZE>
double writetest(char* c) {
  char* v = new char[WSIZE];
  for (size_t i = 0;i < WSIZE; ++i) {
    v[i] = rand() % 256;
  }

  timer ti;
  ti.start();
  size_t steps = msize / WSIZE;
  for (size_t i = 0; i < passes; ++i) {
    for (size_t j = 0;j < steps; ++j) {
      memcpy(c + j * WSIZE, v, WSIZE);
    }
  }
  delete [] v;
  return ti.current_time();
}

template<size_t WSIZE>
double sustained_writetest(char* c, double time) { 
  char* v = new char[WSIZE];
  for (size_t i = 0;i < WSIZE; ++i) {
    v[i] = rand() % 256;
  }
  timer ti;
  ti.start();
  size_t passes = 0;
  size_t steps = msize / WSIZE;
  while(ti.current_time() < time) {
    for (size_t j = 0; j < steps; ++j) {
      memcpy(c + j * WSIZE, v, WSIZE);
    }
    passes++;
  }
  // compute the number of bytes written per second.
  return (double(passes) * msize) / ti.current_time();
}


template<size_t WSIZE>
double sustained_random_writetest(char* c, double time) { 
  char* v = new char[WSIZE];
  for (size_t i = 0;i < WSIZE; ++i) {
    v[i] = rand() % 256;
  }
  std::vector<size_t> randpos;
  // how many?
  size_t numrand = passes * msize / WSIZE;
  randpos.resize(numrand);
  for (size_t i = 0;i < numrand; ++i) {
    randpos[i] = graphlab::random::rand_int(msize - WSIZE - 1);
  }
  
  
  timer ti;
  ti.start();
  size_t passes = 0;
  while(ti.current_time() < time) {
    for (size_t i = 0; i < randpos.size(); ++i) {
        memcpy(c + randpos[i], v, WSIZE);
    }
    passes++;
  }
  // compute the number of bytes written per second.
  return (double(passes) * msize) / ti.current_time();
}

void bench_mmap(mmap_wrapper &wrap) {
  char* c = (char*)wrap.mapped_ptr();
  std::cout << "Seq Write test. " << std::endl;
  std::cout << 16 << "b block: " << writetest<16>(c) << std::endl;
  wrap.sync_all();
  std::cout << 64 << "b block: " << writetest<64>(c) << std::endl;
  wrap.sync_all();
  std::cout << 256 << "b block: " << writetest<256>(c) << std::endl;
  wrap.sync_all();
  std::cout << 1024 << "b block: " << writetest<1024>(c) << std::endl;
  wrap.sync_all();
  std::cout << 4096 << "b block: " << writetest<4096>(c) << std::endl;
  wrap.sync_all();
  std::cout << 8192 << "b block: " << writetest<8192>(c) << std::endl;
  wrap.sync_all();
  std::cout << "Sustained Sequential (64b): " << sustained_writetest<64>(c, 30)/1024/1024 << " MBps" << std::endl;
  wrap.sync_all();

  std::cout << "Random Write test. " << std::endl;
  std::cout << 16 << "b block: " << random_writetest<16>(c) << std::endl;
  wrap.sync_all();
  std::cout << 64 << "b block: " << random_writetest<64>(c) << std::endl;
  wrap.sync_all();
  std::cout << 256 << "b block: " << random_writetest<256>(c) << std::endl;
  wrap.sync_all();
  std::cout << 1024 << "b block: " << random_writetest<1024>(c) << std::endl;
  wrap.sync_all();
  std::cout << 4096 << "b block: " << random_writetest<4096>(c) << std::endl;
  wrap.sync_all();
  std::cout << 8192 << "b block: " << random_writetest<8192>(c) << std::endl;
  wrap.sync_all();
  std::cout << "Sustained Random (64b): " << sustained_random_writetest<64>(c, 30)/1024/1024 << " MBps" << std::endl;
  wrap.sync_all();
 
}



void bench(char* c) {
  std::cout << "Seq Write test. " << std::endl;
  std::cout << 16 << "b block: " << writetest<16>(c) << std::endl;
  std::cout << 64 << "b block: " << writetest<64>(c) << std::endl;
  std::cout << 256 << "b block: " << writetest<256>(c) << std::endl;
  std::cout << 1024 << "b block: " << writetest<1024>(c) << std::endl;
  std::cout << 4096 << "b block: " << writetest<4096>(c) << std::endl;
  std::cout << 8192 << "b block: " << writetest<8192>(c) << std::endl;
  std::cout << "Sustained Sequential (64b): " << sustained_writetest<64>(c, 30)/1024/1024 << " MBps" << std::endl;

  std::cout << "Random Write test. " << std::endl;
  std::cout << 16 << "b block: " << random_writetest<16>(c) << std::endl;
  std::cout << 64 << "b block: " << random_writetest<64>(c) << std::endl;
  std::cout << 256 << "b block: " << random_writetest<256>(c) << std::endl;
  std::cout << 1024 << "b block: " << random_writetest<1024>(c) << std::endl;
  std::cout << 4096 << "b block: " << random_writetest<4096>(c) << std::endl;
  std::cout << 8192 << "b block: " << random_writetest<8192>(c) << std::endl;
  std::cout << "Sustained Random (64b): " << sustained_random_writetest<64>(c, 30)/1024/1024 << " MBps" << std::endl;
}

int main(int argc, char** argv) { 
  std::cout << "Creating testfile" << std::endl;
  mmap_wrapper wrap("testfile.txt", msize);
  std::cout << "MMap perf" << std::endl;
  bench_mmap(wrap);
  wrap.close();

  mmap_wrapper wrap2("testfile.txt", msize, true);
  std::cout << "MMap perf (uncached)" << std::endl;
  bench_mmap(wrap2);
  wrap2.close();
  
  std::cout << "RAM perf" << std::endl;
  char *c = new char[msize];
  bench(c);
  delete [] c;
  
}
