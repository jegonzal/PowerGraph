#ifndef GRAPHLAB_UTIL_TRACEPOINT_HPP
#define GRAPHLAB_UTIL_TRACEPOINT_HPP
#include <iostream>
#include <vector>
#include <string>
#include <graphlab/util/timer.hpp>
#include <graphlab/parallel/atomic.hpp>

#define USE_TRACEPOINT

namespace graphlab{

struct trace_count{
  std::string name;
  std::string description;
  graphlab::atomic<unsigned long long> count;
  graphlab::atomic<unsigned long long> total;
  unsigned long long minimum;
  unsigned long long maximum;
  trace_count():count(0),
                total(0),
                minimum(std::numeric_limits<unsigned long long>::max()),
                maximum(0) { }
                
  inline void incorporate_unsafe(unsigned long long val) {
    ++count.value;
    total.value += val;
    if (val < minimum) minimum = val; 
    if (val > maximum) maximum = val; 
  }
  
  inline void incorporate(unsigned long long val) {
    count.inc();
    total.inc(val);
    while(1) {
      unsigned long long m = minimum;
      if (val > m || graphlab::atomic_compare_and_swap(minimum, m, val)) break;
    }
    while(1) {
      unsigned long long m = maximum;
      if (val < m || graphlab::atomic_compare_and_swap(maximum, m, val)) break;
    }
  }
  
  
  inline void incorporate(const trace_count &val) {
    count.inc(val.count.value);
    total.inc(val.total.value);
    while(1) {
      unsigned long long m = minimum;
      if (val.minimum > m || graphlab::atomic_compare_and_swap(minimum, m, val.minimum)) break;
    }
    while(1) {
      unsigned long long m = maximum;
      if (val.maximum < m || graphlab::atomic_compare_and_swap(maximum, m, val.maximum)) break;
    }
  }
  
  void print(std::ostream& out, unsigned long long tpersec = 0) const;
};

} // namespace

/**
REGISTER_TRACEPOINT(name, desc)
  registers a tracepoint name, and associates it with a longer description.
  name is preferably short. Note that name is not a string, but desc is a string.
  REGISTER_TRACEPOINT should be called serially before any 
  BEGIN_TRACEPOINT/END_TRACEPOINT calls are made. 

BEGIN_TRACEPOINT(name)
END_TRACEPOINT(name)
  Times a block of code. Every END_TRACEPOINT must be matched with a
  BEGIN_TRACEPOINT within the same scope. Tracepoints are parallel.
  Note that name is not a string, but desc is a string.
  
Example Usage:
  REGISTER_TRACEPOINT(classname_someevent, "hello world");
  Then later on...
  BEGIN_TRACEPOINT(classname_someevent)
  ...
  END_TRACEPOINT(classname_someevent)
*/

#ifdef USE_TRACEPOINT
#define REGISTER_TRACEPOINT(name, desc) register_trace(#name, desc); 
#define BEGIN_TRACEPOINT(name) unsigned long long __ ## name ## _trace_ = rdtsc();
#define END_TRACEPOINT(name) store_trace(#name, rdtsc() - __ ## name ## _trace_);
#define END_AND_BEGIN_TRACEPOINT(endname, beginname) unsigned long long __ ## beginname ## _trace_ = rdtsc(); \
                                                     store_trace(#endname, __ ## beginname ## _trace_ - __ ## endname ## _trace_);

#define CREATE_ACCUMULATING_TRACEPOINT(name) graphlab::trace_count __ ## name ## _acc_trace_; \
                                             unsigned long long __ ## name ## _acc_trace_elem_;
#define BEGIN_ACCUMULATING_TRACEPOINT(name) __ ## name ## _acc_trace_elem_ = rdtsc();
#define END_ACCUMULATING_TRACEPOINT(name) __ ## name ## _acc_trace_.incorporate_unsafe(rdtsc() - __ ## name ## _acc_trace_elem_);

#define END_AND_BEGIN_ACCUMULATING_TRACEPOINT(endname, beginname) __ ## beginname ## _acc_trace_elem_ = rdtsc(); \
                                                                  __ ## endname ## _acc_trace_.incorporate_unsafe(__ ## beginname ## _acc_trace_elem_ - __ ## endname ## _acc_trace_elem_)

#define STORE_ACCUMULATING_TRACEPOINT(name) store_trace(#name, __ ## name ## _acc_trace_);
#else
#define REGISTER_TRACEPOINT(name, desc)
#define BEGIN_TRACEPOINT(name) 
#define END_TRACEPOINT(name) 

#define CREATE_ACCUMULATING_TRACEPOINT(name) 
#define BEGIN_ACCUMULATING_TRACEPOINT(name) 
#define END_ACCUMULATING_TRACEPOINT(name)
#define STORE_ACCUMULATING_TRACEPOINT(name)

#define END_AND_BEGIN_ACCUMULATING_TRACEPOINT(endname, beginname)
#define END_AND_BEGIN_TRACEPOINT(endname, beginname) 
#endif           

void register_trace(const char* name, const char* description);
void store_trace(const char* name, unsigned long long val);
void store_trace(const char* name, const graphlab::trace_count& val);
void dump_trace(std::ostream& out, unsigned long long tpersec = 0);

#endif
