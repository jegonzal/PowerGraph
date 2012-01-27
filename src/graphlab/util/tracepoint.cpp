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
  
  void print(std::ostream& out) const {
    out << name << ": " << description << "\n";
    out << "Events:\t" << count.value << "\n";
    out << "Total:\t" << total.value << "\n";
    out << "Mean:\t" << (double)total.value / count.value << "\n";
    out << "Min:\t" << minimum << "\n";
    out << "Max:\t" << maximum << "\n";
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

void dump_trace(std::ostream& out) {
  boost::unordered_map<size_t, trace_count>::const_iterator iter = trace_store.begin();
  while (iter != trace_store.end()) {
    iter->second.print(out);
    ++iter;
  }
}

struct __print_trace_on_destruction__ {
  ~__print_trace_on_destruction__() {
    dump_trace(std::cout);
  }
};

static __print_trace_on_destruction__ print_trace_on_end;
