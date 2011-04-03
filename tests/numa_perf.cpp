#include <iostream>
#include <fstream>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>
#include <boost/bind.hpp>
using namespace graphlab;

std::ofstream fout;
std::vector<void*> mems;
graphlab::barrier *bar;

#define CACHE_LINE_SIZE 64*1024

void cache_flush(void* _start, void* _end) {
  char* start = (char*)_start;
  char* end = (char*)_end;
  while (start < end) {
    asm volatile("clflush (%0)" :: "r"(start));
    start += CACHE_LINE_SIZE;
  }
}

template <typename WordType>
double read_mem_seq(void* m, size_t len) {
  cache_flush(m, (char*)m + len);
  size_t numreads = len / sizeof(WordType);
  // allocate a target
  volatile WordType* readtarget = (volatile WordType*)malloc(len);
  memset((WordType*)readtarget, 1, len);
  WordType* readsource = (WordType*)m;
  timer ti;
  ti.start();
  for (size_t i = 0;i < numreads; ++i) {
    readtarget[i] = readsource[i];
  }
  double rt = ti.current_time();
  free((WordType*)readtarget);
  return rt;
}


// plen should divide len
template <typename WordType>
double read_mem_seq_prefetch(void* m, size_t len, size_t prefetchlen) {
  cache_flush(m, (char*)m + len);
  size_t numprefetch = len / prefetchlen;
  size_t idx_per_prefetch = prefetchlen / sizeof(WordType);
  
  volatile WordType* readtarget = (volatile WordType*)malloc(len);
  memset((WordType*)readtarget, 2, len);
  WordType* readsource = (WordType*)m;

  size_t i = 0;
  timer ti;
  ti.start();

  for (size_t p = 0; p < numprefetch; ++p) {
    prefetch_range(&(readsource[i]), prefetchlen);
    size_t nextidx = i + idx_per_prefetch;
    for (;i < nextidx; ++i) {
      readtarget[i] = readsource[i];
    }
  }
  double rt = ti.current_time();
  free((void*)(WordType*)readtarget);
  return rt;
}


template <typename WordType>
double write_mem_seq(void* m, size_t len) {
  cache_flush(m, (char*)m + len);
  size_t numreads = len / sizeof(WordType);
  // allocate a target
  volatile WordType* readsource = (volatile WordType*)malloc(len);
  memset((WordType*)readsource, 3, len);
  
  WordType* readtarget = (WordType*)m;

  timer ti;
  ti.start();

  for (size_t i = 0;i < numreads; ++i) {
    readtarget[i] = readsource[i];
  }
  double rt = ti.current_time();
  free((void*)(WordType*)readsource);
  return rt;
}


// plen should divide len
template <typename WordType>
double write_mem_seq_prefetch(void* m, size_t len, size_t prefetchlen) {
  cache_flush(m, (char*)m + len);
  size_t numprefetch = len / prefetchlen;
  size_t idx_per_prefetch = prefetchlen / sizeof(WordType);
  
  volatile WordType* readsource = (volatile WordType*)malloc(len);
  memset((WordType*)readsource, 4, len);
  WordType* readtarget = (WordType*)m;
  
  size_t i = 0;
  timer ti;
  ti.start();

  for (size_t p = 0; p < numprefetch; ++p) {
    prefetch_range_write((WordType*)&(readsource[i]), prefetchlen);
    size_t nextidx = i + idx_per_prefetch;
    for (;i < nextidx; ++i) {
      readtarget[i] = readsource[i];
    }
  }
  double rt = ti.current_time();
  free((WordType*)readsource);
  return rt;
}

template <typename WordType>
void runtest(size_t memlen, size_t src, size_t target, bool write, size_t prefetchlen){

  double time = 0.0;
  if (prefetchlen == 0) {
    if (write) time = write_mem_seq<WordType>(mems[target], memlen);
    else time = read_mem_seq<WordType>(mems[target], memlen);
  }
  else {
    if (write) time = write_mem_seq_prefetch<WordType>(mems[target], memlen, prefetchlen);
    else time = read_mem_seq_prefetch<WordType>(mems[target], memlen, prefetchlen);
  }
  size_t wordsize = sizeof(WordType);
  fout << memlen << "\t" << src << "\t" << target << "\t" << (int)write << "\t" << wordsize << "\t" << prefetchlen << "\t" << time << "\n";
  if (write) {
    std::cout << "WRITE: " << "M:" << memlen << "\tS:" << src << "\tT:" << target << "\tW:" << wordsize << "\tP:" << prefetchlen << "\t" << time << "\n";
  }
  else {
    std::cout << "READ : " << "M:" << memlen << "\tS:" << src << "\tT:" << target << "\tW:" << wordsize << "\tP:" << prefetchlen << "\t" << time << "\n";
  }
}

void test_thread(size_t numthreads, size_t thrid, size_t memlen) {
  mems[thrid] = malloc(memlen);
  memset(mems[thrid], 5, memlen);
  bar->wait();
  
  if (thrid == 0) std::cout << "Running all pairs test (Word size)...\n";
  for (size_t p = 0; p < numthreads; ++p) {
    if (thrid == p) {
      std::cout << "Processor " << p << ":" << std::endl;
      for (size_t i = 0;i < numthreads; ++i) {
        runtest<size_t>(memlen, p, i, false, 0);
        runtest<size_t>(memlen, p, i, true, 0);
      }
    }
    bar->wait();
  }

  if (thrid == 0) std::cout << "Running all pairs test (Word size, Prefetch 64KB)...\n";
  for (size_t p = 0; p < numthreads; ++p) {
    if (thrid == p) {
      std::cout << "Processor " << p << ":" << std::endl;
      for (size_t i = 0;i < numthreads; ++i) {
        runtest<size_t>(memlen, p, i, false, 64 * 1024);
        runtest<size_t>(memlen, p, i, true, 64 * 1024);
      }
    }
    bar->wait();
  }

  if (thrid == 0) {
    std::cout << "P0 to P0, varying Word Size, Varying Prefetch...\n";
    size_t pflen[] = {0, 32, 64, 128, 256, 1024, 2 * 1024, 4* 1024, 8 * 1024, 16 * 1024, 32 * 1024, 64 * 1024, 
                      128 * 1024, 256 * 1024, 512 * 1024, 1024 * 1024};
    size_t numpf = sizeof(pflen) / sizeof(size_t);
    for (size_t pf = 0; pf < numpf; ++pf) {
      runtest<char>(memlen, 0, 0, false, pflen[pf]);
      runtest<char>(memlen, 0, 0, true, pflen[pf]);
      runtest<uint16_t>(memlen, 0, 0, false, pflen[pf]);
      runtest<uint16_t>(memlen, 0, 0, true, pflen[pf]);
      runtest<uint32_t>(memlen, 0, 0, false, pflen[pf]);
      runtest<uint32_t>(memlen, 0, 0, true, pflen[pf]);
      runtest<uint64_t>(memlen, 0, 0, false, pflen[pf]);
      runtest<uint64_t>(memlen, 0, 0, true, pflen[pf]);
      runtest<double>(memlen, 0, 0, false, pflen[pf]);
      runtest<double>(memlen, 0, 0, true, pflen[pf]);
    }
    
    std::cout << "P0 to P1, varying Word Size, Varying Prefetch...\n";
    if (numthreads > 1) {
      for (size_t pf = 0; pf < numpf; ++pf) {
      runtest<char>(memlen, 0, 1, false, pflen[pf]);
      runtest<char>(memlen, 0, 1, true, pflen[pf]);
      runtest<uint16_t>(memlen, 0, 1, false, pflen[pf]);
      runtest<uint16_t>(memlen, 0, 1, true, pflen[pf]);
      runtest<uint32_t>(memlen, 0, 1, false, pflen[pf]);
      runtest<uint32_t>(memlen, 0, 1, true, pflen[pf]);
      runtest<uint64_t>(memlen, 0, 1, false, pflen[pf]);
      runtest<uint64_t>(memlen, 0, 1, true, pflen[pf]);
      runtest<double>(memlen, 0, 1, false, pflen[pf]);
      runtest<double>(memlen, 0, 1, true, pflen[pf]);      }
    }
  }

  bar->wait();
}


int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "numa_perf [numthreads] [memtestlen in MB]\n";
    return 0;
  }
  size_t numthreads = atoi(argv[1]);
  size_t testlen = atoi(argv[2]) * 1024 * 1024;
  bar = new graphlab::barrier(numthreads);
  mems.resize(numthreads);
  fout.open("log.txt");
  thread_group grp;
  for (size_t i = 0; i < numthreads; ++i) {
    launch_in_new_thread(grp,
                        boost::bind(test_thread, numthreads, i, testlen), 
                        i);
  }
  grp.join();

}
