#ifndef GRAPHLAB_UTIL_TRACEPOINT_HPP
#define GRAPHLAB_UTIL_TRACEPOINT_HPP
#include <graphlab/util/timer.hpp>
#include <iostream>

//#define USE_TRACEPOINT

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
#else
#define REGISTER_TRACEPOINT(name, desc)
#define BEGIN_TRACEPOINT(name) 
#define END_TRACEPOINT(name) 
#endif           

void register_trace(const char* name, const char* description);
void store_trace(const char* name, unsigned long long val);
void dump_trace(std::ostream& out, unsigned long long tpersec = 0);

#endif
