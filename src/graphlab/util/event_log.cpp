#include <graphlab/util/event_log.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/logger/assertions.hpp>

#define PROGRESS_BAR_WIDTH 40
#define BAR_CHARACTER '#'

namespace graphlab {

void event_log::initialize(std::ostream &ostrm,
                           size_t flush_interval_ms,
                           event_print_type event_print) {
  m.lock();
  out = &ostrm;
  flush_interval = flush_interval_ms;
  print_method = event_print;
  prevtime = 0;
  ti.start();
  cond.signal();
  std::cout << "init" << std::endl;
  m.unlock();
}

event_log::~event_log() {
  finished = true;
  m.lock();
  cond.signal();
  m.unlock();
  printing_thread.join();
}

void event_log::close() {
  out = NULL;
  m.lock();
  flush_interval = 0;
  m.unlock();
}

void event_log::add_event_type(unsigned char eventid,
                               std::string description) {
  descriptions[eventid] = description;
  max_desc_length = std::max(max_desc_length, description.length());
  ASSERT_MSG(max_desc_length <= 30, "Event Description length must be <= 30 characters");
  counters[eventid].value = 0;
  hascounter.set_bit(eventid);
}


void event_log::flush() {
  uint32_t pos;
  if (!hascounter.first_bit(pos)) return;
  double curtime = ti.current_time_millis();
  double timegap = curtime - prevtime;
  prevtime = curtime;

  if (hasevents == false && noeventctr == 1) return;
  
  bool found_events = false;
  (*out) << "Time: " << "+"<<timegap << "\t" << curtime << "\n";
  if (print_method == NUMBER) {
    do {
      size_t ctrval = counters[pos].exchange(0);
      found_events = found_events | ctrval > 0;
      (*out) << pos  << ":\t" << ctrval << "\t" << 1000 * ctrval / timegap << " /s\n";
    } while(hascounter.next_bit(pos));
  }
  else if (print_method == DESCRIPTION) {
    do {
      size_t ctrval = counters[pos].exchange(0);
      found_events = found_events | ctrval > 0;
      (*out) << descriptions[pos]  << ":\t" << ctrval << "\t" << 1000 * ctrval / timegap << " /s\n";
    } while(hascounter.next_bit(pos));
  }
  else if (print_method == RATE_BAR) {
    char spacebuf[60];
    char pbuf[61];
    memset(spacebuf, ' ', PROGRESS_BAR_WIDTH);
    memset(pbuf, BAR_CHARACTER, 60);
    do {
      size_t ctrval = counters[pos].exchange(0);
      found_events = found_events | ctrval > 0;
      maxcounter[pos] = std::max(maxcounter[pos], ctrval);
      size_t barlen = 0;
      size_t mc = maxcounter[pos]; 
      if (mc > 0) barlen = ctrval * PROGRESS_BAR_WIDTH / mc;
      if (barlen > PROGRESS_BAR_WIDTH) barlen = PROGRESS_BAR_WIDTH;
      
      pbuf[barlen] = '\0';
      spacebuf[max_desc_length - descriptions[pos].length() + 1] = 0;
      (*out) << descriptions[pos]  << spacebuf << "|" << pbuf;
      spacebuf[max_desc_length - descriptions[pos].length() + 1] =' ';
      pbuf[barlen] = BAR_CHARACTER;
      // now print the remaining spaces
      spacebuf[PROGRESS_BAR_WIDTH - barlen] = '\0';
      (*out) << spacebuf << "| " << ctrval << " : " << mc << " /s\n";
      spacebuf[PROGRESS_BAR_WIDTH - barlen] = ' ';
      
    } while(hascounter.next_bit(pos));
  }
  if (found_events == false) {
      ++noeventctr;
  }
  else {
      noeventctr = 0;
  }
  hasevents = false;
  out->flush();
}

void event_log::thread_loop() {
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