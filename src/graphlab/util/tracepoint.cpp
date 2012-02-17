#include <limits>
#include <string>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/tracepoint.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <boost/unordered_map.hpp>


namespace graphlab {

static unsigned long long rtdsc_ticks_per_sec = 0; 
mutex rtdsc_ticks_per_sec_mutex;
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

unsigned long long estimate_ticks_per_second() {
  rtdsc_ticks_per_sec_mutex.lock();
  if (rtdsc_ticks_per_sec == 0) {
    unsigned long long tstart = rdtsc();
    graphlab::my_sleep(1);
    unsigned long long tend = rdtsc();
    rtdsc_ticks_per_sec = tend - tstart;
  }
  rtdsc_ticks_per_sec_mutex.unlock();
  return rtdsc_ticks_per_sec;
}


trace_count::~trace_count() {
#ifdef USE_TRACEPOINT
  print(std::cout, estimate_ticks_per_second());
#endif
}

} // namespace graphlab

