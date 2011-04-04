#include <iostream>
#include <fstream>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>
#include <boost/bind.hpp>
using namespace graphlab;

std::ofstream fout;
std::vector<void*> mems;
graphlab::barrier *bar;

#define CACHE_LINE_SIZE 64


inline void prefetch_range_2(void *addr, size_t len) {
  char *cp;
  char *end = (char*)(addr) + len;

  for (cp = (char*)(addr); cp < end; cp += CACHE_LINE_SIZE) {
    __builtin_prefetch(cp, 0); 
  }
}

inline void prefetch_range_write_2(void *addr, size_t len) {
  char *cp;
  char *end = (char*)(addr) + len;

  for (cp = (char*)(addr); cp < end; cp += CACHE_LINE_SIZE) {
    __builtin_prefetch(cp, 1);
  }
}

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
  volatile WordType readtarget;
  WordType* readsource = (WordType*)m;
  timer ti;
  ti.start();
  for (size_t i = 0;i < numreads; ++i) {
    readtarget = readsource[i];
  }
  double rt = ti.current_time();
  return rt;
}



template <typename WordType>
double write_mem_seq(void* m, size_t len) {
  cache_flush(m, (char*)m + len);
  size_t numreads = len / sizeof(WordType);
  // allocate a target
  volatile WordType readsource;
  
  WordType* readtarget = (WordType*)m;

  timer ti;
  ti.start();

  for (size_t i = 0;i < numreads; ++i) {
    readtarget[i] = readsource;
  }
  double rt = ti.current_time();
  return rt;
}

// plen should divide len
template <typename WordType>
double read_mem_seq_prefetch(void* m, size_t len, size_t prefetchgap) {
  cache_flush(m, (char*)m + len);
  size_t numreads = len / sizeof(WordType);
 
  volatile WordType readtarget;
  WordType* readsource = (WordType*)m;

  size_t i = 0;
  timer ti;
  ti.start();
  for (i = 0;i < numreads - prefetchgap; ++i) {
    __builtin_prefetch(&(readsource[i+prefetchgap]), 0); 
    readtarget = readsource[i];
  }
  for (;i < numreads; ++i) {
    readtarget = readsource[i];
  }
 
  double rt = ti.current_time();
  return rt;
}


// plen should divide len
template <typename WordType>
double write_mem_seq_prefetch(void* m, size_t len, size_t prefetchgap) {
  cache_flush(m, (char*)m + len);
  size_t numreads = len / sizeof(WordType);
 
  volatile WordType readsource;
  WordType* readtarget = (WordType*)m;
  
  size_t i = 0;
  timer ti;
  ti.start();
  for (i = 0;i < numreads - prefetchgap; ++i) {
    __builtin_prefetch(&(readtarget[i+prefetchgap]), 1); 
    readtarget[i] = readsource;
  }
  for (;i < numreads; ++i) {
    readtarget[i] = readsource;
  }
  double rt = ti.current_time();
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

  if (thrid == 0) std::cout << "Running all pairs test (Word size, Prefetch 64K)...\n";
  for (size_t p = 0; p < numthreads; ++p) {
    if (thrid == p) {
      std::cout << "Processor " << p << ":" << std::endl;
      for (size_t i = 0;i < numthreads; ++i) {
        runtest<size_t>(memlen, p, i, false, 64);
        runtest<size_t>(memlen, p, i, true, 64);
      }
    }
    bar->wait();
  }

  if (thrid == 0) {
    std::cout << "P0 to P0, varying Word Size, Varying Prefetch Gap...\n";
    for (size_t pf = 0; pf <= 128; pf+=8) {
      runtest<char>(memlen, 0, 0, false, pf);
      runtest<char>(memlen, 0, 0, true, pf);
      runtest<uint16_t>(memlen, 0, 0, false, pf);
      runtest<uint16_t>(memlen, 0, 0, true, pf);
      runtest<uint32_t>(memlen, 0, 0, false, pf);
      runtest<uint32_t>(memlen, 0, 0, true, pf);
      runtest<uint64_t>(memlen, 0, 0, false, pf);
      runtest<uint64_t>(memlen, 0, 0, true, pf);
    }
    
    std::cout << "P0 to P1, varying Word Size, Varying Prefetch...\n";
    if (numthreads > 1) {
      for (size_t pf = 0; pf <= 128; pf+=8) {
      runtest<char>(memlen, 0, 1, false, pf);
      runtest<char>(memlen, 0, 1, true, pf);
      runtest<uint16_t>(memlen, 0, 1, false, pf);
      runtest<uint16_t>(memlen, 0, 1, true, pf);
      runtest<uint32_t>(memlen, 0, 1, false, pf);
      runtest<uint32_t>(memlen, 0, 1, true, pf);
      runtest<uint64_t>(memlen, 0, 1, false, pf);
      runtest<uint64_t>(memlen, 0, 1, true, pf);
      }
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
  fout.close();
}
