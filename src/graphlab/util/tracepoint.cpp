#include <limits>
#include <string>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/tracepoint.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <boost/unordered_map.hpp>

static boost::unordered_map<size_t, graphlab::trace_count> trace_store;


namespace graphlab {
  
void trace_count::print(std::ostream& out, unsigned long long tpersec) const {
  if (tpersec == 0) {
    out << name << ": " << description << "\n";
    out << "Events:\t" << count.value << "\n";
    out << "Total:\t" << total.value << "ticks \n";
    if (count.value > 0) {
      out << "Mean:\t" << (double)total.value / count.value << "ticks \n";
      out << "Min:\t" << minimum << "ticks \n";
      out << "Max:\t" << maximum << "ticks \n";
    }
  }
  else {
    double tperms = (double)tpersec / 1000;
    out << name << ": " << description << "\n";
    out << "Events:\t" << count.value << "\n";
    out << "Total:\t" << (double)total.value / tperms << " ms \n";
    if (count.value > 0) {
      out << "Mean:\t" << (double)total.value / count.value / tperms << " ms \n";
      out << "Min:\t" << (double)minimum / tperms << " ms \n";
      out << "Max:\t" << (double)maximum / tperms << " ms \n";
    }
  }
}

} // namespace graphlab

void register_trace(const char* name, const char* description) {
  size_t hashval = boost::hash_value(name);
  trace_store[hashval] = graphlab::trace_count();
  trace_store[hashval].name = name;
  trace_store[hashval].description = description;
}

void store_trace(const char* name, unsigned long long val) {
  size_t hashval = boost::hash_value(name);
  boost::unordered_map<size_t, graphlab::trace_count>::iterator iter = trace_store.find(hashval);
  ASSERT_TRUE(iter != trace_store.end());
  graphlab::trace_count& tc = iter->second;
  tc.incorporate(val); 
}

void store_trace(const char* name, const graphlab::trace_count &val) {
  size_t hashval = boost::hash_value(name);
  boost::unordered_map<size_t, graphlab::trace_count>::iterator iter = trace_store.find(hashval);
  ASSERT_TRUE(iter != trace_store.end());
  graphlab::trace_count& tc = iter->second;
  tc.incorporate(val); 
}

void dump_trace(std::ostream& out, unsigned long long tpersec) {
  boost::unordered_map<size_t, graphlab::trace_count>::const_iterator iter = trace_store.begin();
  while (iter != trace_store.end()) {
    iter->second.print(out, tpersec);
    ++iter;
  }
}

#ifdef USE_TRACEPOINT

unsigned long long estimate_ticks_per_second() {
  unsigned long long tstart = rdtsc();
  graphlab::my_sleep(1);
  unsigned long long tend = rdtsc();
  return tend - tstart;
}

struct __print_trace_on_destruction__ {
  ~__print_trace_on_destruction__() {
    std::cout << "estimating ticks per second ..." << std::endl;
    size_t tpersec = estimate_ticks_per_second();
    std::cout << tpersec << " ticks per second." << std::endl;
    dump_trace(std::cout, tpersec);
  }
};



static __print_trace_on_destruction__ print_trace_on_end;
#endif
