#include <limits>
#include <string>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/tracepoint.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <boost/unordered_map.hpp>


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
  
  void incorporate(unsigned long long val) {
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
  
  void print(std::ostream& out, unsigned long long tpersec = 0) const {
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
};

static boost::unordered_map<size_t, trace_count> trace_store;

void register_trace(const char* name, const char* description) {
  size_t hashval = boost::hash_value(name);
  trace_store[hashval] = trace_count();
  trace_store[hashval].name = name;
  trace_store[hashval].description = description;
}

void store_trace(const char* name, unsigned long long val) {
  size_t hashval = boost::hash_value(name);
  boost::unordered_map<size_t, trace_count>::iterator iter = trace_store.find(hashval);
  ASSERT_TRUE(iter != trace_store.end());
  trace_count& tc = iter->second;
  tc.incorporate(val); 
}

void dump_trace(std::ostream& out, unsigned long long tpersec) {
  boost::unordered_map<size_t, trace_count>::const_iterator iter = trace_store.begin();
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
