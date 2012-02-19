#include <limits>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/logger/assertions.hpp>

#define EVENT_BAR_WIDTH 40

namespace graphlab {

  
void dist_event_log::initialize(distributed_control& dc,
                           std::ostream &ostrm,
                           size_t flush_interval_ms,
                           event_print_type event_print) {
  rmi = new dc_dist_object<dist_event_log>(dc, this);
  rmi->barrier();
  m.lock();
  out = &ostrm;
  flush_interval = flush_interval_ms;
  if (rmi->procid() > 0) {
    flush_interval /= 10;
  }
  print_method = event_print;
  prevtime = 0;
  ti.start();
  cond.signal();
  std::cout << "init" << std::endl;
  m.unlock();
}

dist_event_log::~dist_event_log() {
  if (!finished) destroy();
}

void dist_event_log::destroy() {
  uint32_t pos;
  if (rmi->procid() == 0 && hascounter.first_bit(pos)) {
    do {
      (*out) << descriptions[pos]  << ":\t" << totalcounter[pos].value << " Events\n";
    } while(hascounter.next_bit(pos));
  }
  finished = true;
  m.lock();
  cond.signal();
  m.unlock();
  printing_thread.join();
  delete rmi;
}

void dist_event_log::close() {
  out = NULL;
  m.lock();
  flush_interval = 0;
  m.unlock();
}

void dist_event_log::add_event_type(unsigned char eventid,
                               std::string description) {
  descriptions[eventid] = description;
  max_desc_length = std::max(max_desc_length, description.length());
  ASSERT_MSG(max_desc_length <= 30, "Event Description length must be <= 30 characters");
  counters[eventid].value = 0;
  hascounter.set_bit(eventid);
  if (rmi->procid() == 0) {
    globalcounters[eventid].resize(rmi->numprocs());
  }
}

void dist_event_log::accumulate_event_aggregator(procid_t proc,
                                                 unsigned char eventid,
                                                 size_t count) {
  hasevents = true;
  totalcounter[eventid].inc(count);
  globalcounters[eventid][proc].inc(count);
}
  
struct counter_statistics{
  size_t minimum;
  size_t maximum;
  size_t average;
  size_t total;
};
  
static void compute_statistics(std::vector<atomic<size_t> > &vec,
                               counter_statistics& c) {
  c.minimum = std::numeric_limits<size_t>::max();
  c.maximum = std::numeric_limits<size_t>::min();
  c.average = 0;
  c.total = 0;
  for (size_t i = 0; i < vec.size(); ++i) {
    size_t ctr = vec[i].exchange(0);
    c.minimum = std::min(c.minimum, ctr);
    c.maximum = std::max(c.minimum, ctr);
    c.total += ctr;
  }
  c.average = c.total / vec.size();
}



static void print_bar(std::ostream& out,
                             size_t val,
                             size_t len) {
  // compute the bar lengths
  if (len == 0) return;
  val = val * EVENT_BAR_WIDTH / len;
  if (val > EVENT_BAR_WIDTH) val = EVENT_BAR_WIDTH;
  size_t i = 0;
  for (i = 0; i < val; ++i) out.put('#');
  for (; i < EVENT_BAR_WIDTH; ++i) out.put(' ');
}

static void print_triple_bar(std::ostream& out,
                             size_t low,
                             size_t mid,
                             size_t high,
                             size_t len) {
  // compute the bar lengths
  if (len == 0) return;
  low = low * EVENT_BAR_WIDTH / len;
  if (low > EVENT_BAR_WIDTH) low = EVENT_BAR_WIDTH;

  mid = mid * EVENT_BAR_WIDTH / len;
  if (mid > EVENT_BAR_WIDTH) mid = EVENT_BAR_WIDTH;
  
  high = high * EVENT_BAR_WIDTH / len;
  if (high > EVENT_BAR_WIDTH) high = EVENT_BAR_WIDTH;
  size_t i = 0;
  for (i = 0; i < low; ++i) out.put('-');
  for (; i < mid; ++i) out.put('*');
  for (; i < high; ++i) out.put('#');
  for (; i < EVENT_BAR_WIDTH; ++i) out.put(' ');
}
  
void dist_event_log::print_log() {
  uint32_t pos;
  if (!hascounter.first_bit(pos)) return;
  double curtime = ti.current_time_millis();
  double timegap = curtime - prevtime;
  prevtime = curtime;

  if (hasevents == false && noeventctr == 1) return;

  counter_statistics stats[256];
  // accumulate the statistics for printing
  do {
    compute_statistics(globalcounters[pos], stats[pos]);
  } while(hascounter.next_bit(pos));
  
  bool found_events = false;
  (*out) << "Time: " << "+"<<timegap << "\t" << curtime << "\n";
  
  // reset the counter
  hascounter.first_bit(pos);
  if (print_method == NUMBER) {
    do {
      found_events = found_events || stats[pos].total > 0;
      (*out) << pos  << ":\t" << stats[pos].minimum << "\t"
             << stats[pos].average << "\t" << stats[pos].maximum << "\t"
             << stats[pos].total << "\t" << 1000 * stats[pos].total / timegap << " /s\n";
    } while(hascounter.next_bit(pos));
  }
  else if (print_method == DESCRIPTION) {
    do {
      found_events = found_events || stats[pos].total > 0;
      (*out) << descriptions[pos]  << ":\t" << stats[pos].minimum << "\t"
             << stats[pos].average << "\t" << stats[pos].maximum << "\t"
             << stats[pos].total << "\t" << 1000 * stats[pos].total / timegap << " /s\n";
    } while(hascounter.next_bit(pos));
  }
  else if (print_method == RATE_BAR) {
    char spacebuf[60];
    memset(spacebuf, ' ', EVENT_BAR_WIDTH);
    do {
      found_events = found_events || stats[pos].total > 0;
      maxcounter[pos] = std::max(maxcounter[pos], stats[pos].maximum);
      maxproc_counter[pos] = std::max(maxcounter[pos], stats[pos].maximum);
      
      // print the description prefix
      spacebuf[max_desc_length - descriptions[pos].length() + 1] = 0;
      (*out) << descriptions[pos]  << spacebuf << "|";

      print_triple_bar(*out,
                       stats[pos].minimum,
                       stats[pos].average,
                       stats[pos].maximum,
                       maxproc_counter[pos]);
      (*out) << "| " << stats[pos].minimum << " " <<
                        stats[pos].average << " " <<
                        stats[pos].maximum << "\n";

      // print description prefix again
      (*out) << descriptions[pos]  << spacebuf << "|";

      print_bar(*out, stats[pos].total, maxcounter[pos]);
      // reset the space buffer
      spacebuf[max_desc_length - descriptions[pos].length() + 1] =' ';

       (*out) << "| " << stats[pos].total << " : " << maxcounter[pos] << " /s\n\n";
    } while(hascounter.next_bit(pos));
  }
  if (found_events == false) {
      ++noeventctr;
  }
  else {
      noeventctr = 0;
  }
  hasevents = false;
  (*out) << std::endl;
}

void dist_event_log::flush() {
  if (rmi->procid() == 0) {
    // move counters to the aggregator
    uint32_t pos;
    if (!hascounter.first_bit(pos)) return;
    do {
      size_t ctrval = counters[pos].exchange(0);
      if (ctrval != 0) {
        accumulate_event_aggregator(0, pos, ctrval);
      }
    } while(hascounter.next_bit(pos));
    print_log();
  }
  else {
    uint32_t pos;
    if (!hascounter.first_bit(pos)) return;
    do {
      size_t ctrval = counters[pos].exchange(0);
      if (ctrval != 0) {
        rmi->remote_call(0,
                         &dist_event_log::accumulate_event_aggregator,
                         rmi->procid(),
                         (unsigned char)pos, ctrval);
      }
    } while(hascounter.next_bit(pos));
  }
}

void dist_event_log::thread_loop() {
  m.lock();
  while(!finished) {
    if (flush_interval == 0) {
      cond.wait(m);
    }
    else {
      m.unlock();
      my_sleep_ms(flush_interval);
      m.lock();
      if (flush_interval > 0) flush();
    }
  }
  m.unlock();
}

} // namespace