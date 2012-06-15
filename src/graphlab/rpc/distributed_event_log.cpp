/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */




#include <limits>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/logger/assertions.hpp>


namespace graphlab {

  struct distributed_event_log_state {
    distributed_event_log_state():
      eventlog_file_open(false),event_timer_started(false) { }
      
    std::ofstream eventlog_file;
    mutex eventlog_file_mutex;
    bool eventlog_file_open;

    timer event_timer;
    bool event_timer_started;
    mutex event_timer_mutex;
  };

  static distributed_event_log_state& get_state() {
    static distributed_event_log_state state;
    return state;
  }
  
  void dist_event_log::initialize(distributed_control& dc,
                                  std::ostream &ostrm,
                                  size_t flush_interval_ms,
                                  event_print_type event_print) {
    rmi = new dc_dist_object<dist_event_log>(dc, this);
    rmi->barrier();

    m.lock();
    out = &ostrm;
    flush_interval = flush_interval_ms;
    // other processors synchronize here more frequently.
    if (rmi->procid() > 0) {
      flush_interval /= 10;
    }
    print_method = event_print;
  
    get_state().event_timer_mutex.lock();
    if (get_state().event_timer_started == false) {
      get_state().event_timer_started = true;
      get_state().event_timer.start();
    }
    get_state().event_timer_mutex.unlock();
    prevtime = get_state().event_timer.current_time_millis();

    cond.signal();
    m.unlock();
  
  
    if (event_print == LOG_FILE) {
      if (dc.procid() == 0) {
        get_state().eventlog_file_mutex.lock();
        if (!get_state().eventlog_file_open) {
          get_state().eventlog_file_open = true;
          get_state().eventlog_file.open("eventlog.txt");
        }
        out = &get_state().eventlog_file;
        get_state().eventlog_file_mutex.unlock();
      }
    }
  }

  dist_event_log::~dist_event_log() {
    if (!finished) destroy();
  }


  void dist_event_log::flush_and_reset_counters() {
    flush();
    rmi->full_barrier();
    uint32_t pos;
    if (print_method != LOG_FILE) {
      if (rmi->procid() == 0 && hascounter.first_bit(pos)) {
        do {
          size_t r = totalcounter[pos].value;
          totalcounter[pos].dec(r);
          (*out) << descriptions[pos]  << ":\t" << r << " Events\n";
        } while(hascounter.next_bit(pos));
      }
    }
    else {
      if (rmi->procid() == 0 && hascounter.first_bit(pos)) {
        do {
          size_t r = totalcounter[pos].value;
          totalcounter[pos].dec(r);
          std::cout << descriptions[pos]  << ":\t" << r << " Events\n";
        } while(hascounter.next_bit(pos));
      }
    }
  }


  void dist_event_log::destroy() {
    m.lock();
    finished = true;
    cond.signal();
    m.unlock();
    printing_thread.join();
    flush_and_reset_counters();
    rmi->barrier();
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

  void dist_event_log::add_immediate_event_type(unsigned char eventid,
                                                std::string description) {
    descriptions[eventid] = description;
    max_desc_length = std::max(max_desc_length, description.length());
    ASSERT_MSG(max_desc_length <= 30, "Event Description length must be <= 30 characters");
    counters[eventid].value = 0;
    if (rmi->procid() == 0) {
      globalcounters[eventid].resize(rmi->numprocs());
    }
  }


  void dist_event_log::accumulate_event_aggregator(size_t proc,
                                                   unsigned char eventid,
                                                   size_t count) {
    hasevents = true;
    totalcounter[eventid].inc(count);
    globalcounters[eventid][proc].inc(count);
  }

  void dist_event_log::immediate_event_aggregator(const std::vector<std::pair<unsigned char, size_t> >& im) {
    hasevents = true;
    m.lock();
    std::copy(im.begin(), im.end(), 
              std::inserter(immediate_events, immediate_events.end()));
    m.unlock();
  }

  void dist_event_log::immediate_event(unsigned char eventid) {
    m.lock();
    immediate_events.push_back(std::make_pair(eventid, get_state().event_timer.current_time_millis()));
    m.unlock();
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
      c.maximum = std::max(c.maximum, ctr);
      c.total += ctr;
    }
    c.average = c.total / vec.size();
  }

  void dist_event_log::print_log() {
    uint32_t pos;
    if (!hascounter.first_bit(pos)) return;
    double curtime = get_state().event_timer.current_time_millis();
    double timegap = curtime - prevtime;
    prevtime = curtime;

    if (hasevents == false && noeventctr == 1) return;

    counter_statistics stats[256];
    // accumulate the statistics for printing
    do {
      compute_statistics(globalcounters[pos], stats[pos]);
    } while(hascounter.next_bit(pos));
  
    bool found_events = false;
  
    // reset the counter
    hascounter.first_bit(pos);
    if (print_method == DESCRIPTION) {
      do {
        found_events = found_events || stats[pos].total > 0;
        (*out) << descriptions[pos]  << ":\t" << curtime << "\t" << stats[pos].minimum << "\t"
               << stats[pos].average << "\t" << stats[pos].maximum << "\t"
               << stats[pos].total << "\t" << 1000 * stats[pos].total / timegap << " /s\n";
      } while(hascounter.next_bit(pos));
      if (!immediate_events.empty()) { 
        std::vector<std::pair<unsigned char, size_t> > cur;
        cur.swap(immediate_events);
        for (size_t i = 0;i < cur.size(); ++i) {
          (*out) << descriptions[cur[i].first] << ":\t" << cur[i].second << "\t" << -1 << "\t"
                 << -1 << "\t" << -1 << "\t" << -1 << "\t" << 0 << " /s\n";
        }
      }
      out->flush();
    }
    else if (print_method == LOG_FILE) {
      get_state().eventlog_file_mutex.lock();
      do {
        found_events = found_events || stats[pos].total > 0;
        (*out) << descriptions[pos]  << ":\t" << curtime << "\t" << stats[pos].minimum << "\t"
               << stats[pos].average << "\t" << stats[pos].maximum << "\t"
               << stats[pos].total << "\t" << 1000 * stats[pos].total / timegap << "\n";
      } while(hascounter.next_bit(pos));
      if (!immediate_events.empty()) { 
        std::vector<std::pair<unsigned char, size_t> > cur;
        cur.swap(immediate_events);
        for (size_t i = 0;i < cur.size(); ++i) {
          (*out) << descriptions[cur[i].first] << ":\t" << cur[i].second << "\t" << -1 << "\t"
                 << -1 << "\t" << -1 << "\t" << -1 << "\t" << 0 << " /s\n";
        }
      }
      out->flush();
      get_state().eventlog_file_mutex.unlock();
    }
    if (found_events == false) {
      ++noeventctr;
    }
    else {
      noeventctr = 0;
    }
    hasevents = false;
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
          rmi->control_call(0,
                            &dist_event_log::accumulate_event_aggregator,
                            rmi->procid(),
                            (unsigned char)pos, ctrval);
        }
      } while(hascounter.next_bit(pos));
      if (!immediate_events.empty()) {
        std::vector<std::pair<unsigned char, size_t> > cur;
        cur.swap(immediate_events);
        rmi->control_call(0, &dist_event_log::immediate_event_aggregator, cur);
      }
    }
  }

  void dist_event_log::thread_loop() {
    m.lock();
    while(!finished) {
      if (flush_interval == 0) {
        cond.wait(m);
      }
      else {
        cond.timedwait_ns(m, flush_interval * 1000000);
        flush();
      }
    }
    m.unlock();
  }

} // namespace
