/*  
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



#include <pthread.h>
#include <string>
#include <limits>
#include <cfloat>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/ui/metrics_server.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {


// predeclaration the metric server handlers
static std::pair<std::string, std::string> 
metric_names_json(std::map<std::string, std::string>& vars);

static std::pair<std::string, std::string> 
metric_aggregate_json(std::map<std::string, std::string>& vars);


uint32_t distributed_event_logger::allocate_log_entry(log_group* group) {
  log_entry_lock.lock();
  uint32_t id = 0;
  if (has_log_entry.first_zero_bit(id) == false) {
    logger(LOG_FATAL, "More than 256 Log entries created. "
        "New log entries cannot be created");
    // does not return
  }
  logs[id] = group;
  has_log_entry.set_bit(id);
  log_entry_lock.unlock();
  return id;
}

event_log_thread_local_type* distributed_event_logger::get_thread_counter_ref() {
  void* v = pthread_getspecific(key);
  if (v == NULL) {
    // allocate a new thread local entry
    event_log_thread_local_type* entry = new event_log_thread_local_type;
    // set all values to 0
    for (size_t i = 0; i < MAX_LOG_SIZE; ++i) entry->values[i] = 0;
    // cast and write it to v. We need it later. 
    // and set the thread local store
    v = (void*)(entry);
    pthread_setspecific(key, v);

    // register the key entry against the logger
    thread_local_count_lock.lock();
    // find an unused entry
    uint32_t b = 0;
    if (thread_local_count_slots.first_zero_bit(b) == false) {
      logger(LOG_FATAL, "More than 1024 active threads. "
          "Log counters cannot be created");
      // does not return
    }
    entry->thlocal_slot = b;
    thread_local_count[b] = entry;
    thread_local_count_slots.set_bit(b);
    thread_local_count_lock.unlock();
  }

  event_log_thread_local_type* entry = (event_log_thread_local_type*)(v);
  return entry;
}

/**
 * Receives the log information from each machine
 */    
void distributed_event_logger::rpc_collect_log(size_t srcproc, double srctime,
    std::vector<double> srccounts) {
  foreach(uint32_t log, has_log_entry) {
    logs[log]->lock.lock();
    // insert the new counts
    logs[log]->machine[srcproc].push_back(log_entry(srccounts[log], srctime));
    logs[log]->lock.unlock();
  }
}

void distributed_event_logger::collect_instantaneous_log() {
 foreach(uint32_t log, has_log_entry) {
    if (logs[log]->logtype == log_type::INSTANTANEOUS) {
      logs[log]->lock.lock();
      // for each log entry which is a callback entry
      // call the callback to get the counts
      if (logs[log]->is_callback_entry) {
        logs[log]->sum_of_instantaneous_entries += logs[log]->callback();
        ++logs[log]->count_of_instantaneous_entries;
      }
      else {
        // sum it across all the threads
        foreach(uint32_t thr, thread_local_count_slots) {
          logs[log]->sum_of_instantaneous_entries += thread_local_count[thr]->values[log];
        }
        ++logs[log]->count_of_instantaneous_entries;
      }
      logs[log]->lock.unlock();
    }
  }
}

/** 
 *  Collects the machine level
 *  log entry. and sends it to machine 0
 */
void distributed_event_logger::local_collect_log() {
  // put together an aggregate of all counters 
  std::vector<double> combined_counts(MAX_LOG_SIZE, 0);

  // for each thread and for each log entry which is 
  // not a callback entry. Accumulate the number of counts
  //
  foreach(uint32_t log, has_log_entry) {
    logs[log]->lock.lock();
    // cimulative entry. just add across all threads
    if (logs[log]->logtype == log_type::CUMULATIVE) {
      foreach(uint32_t thr, thread_local_count_slots) {
        double* current_thread_counts = thread_local_count[thr]->values;
        combined_counts[log] += current_thread_counts[log];
      }
    }
    else {
      // take the average 
      combined_counts[log] = logs[log]->sum_of_instantaneous_entries / 
                                logs[log]->count_of_instantaneous_entries;
      logs[log]->sum_of_instantaneous_entries = 0;
      logs[log]->count_of_instantaneous_entries = 0;
    }
    logs[log]->lock.unlock();
  }

  // send to machine 0
  if (rmi->procid() != 0) {
    rmi->remote_call(0, &distributed_event_logger::rpc_collect_log,
        (size_t)rmi->procid(), (double)ti.current_time(), combined_counts);
  }
  else {
    rpc_collect_log((size_t)0, (double)ti.current_time(), combined_counts);
  }
}

// Called only by machine 0 to get the aggregate log
void distributed_event_logger::build_aggregate_log() {
  ASSERT_EQ(rmi->procid(), 0);
  double current_time = ti.current_time();
  foreach(uint32_t log, has_log_entry) {
    logs[log]->lock.lock();
    // what is the previous time the aggregate was computed?
    double prevtime = 0;
    if (logs[log]->aggregate.size() > 0) {
      prevtime = logs[log]->aggregate.rbegin()->time;
    }
    // if it is a CUMULATIVE log, take the latest entry from each machine
    // if it is an INSTANTANEOUS log, take the average of the last times.

    double sum = 0;
    if (logs[log]->logtype == log_type::CUMULATIVE) {
      for (procid_t p = 0; p < logs[log]->machine.size(); ++p) {
        if (logs[log]->machine[p].size() > 0) {
          sum += logs[log]->machine[p].rbegin()->value;
        }
      }
    }
    else {
      // sum all the machine logs from prevtime to current time
      for (procid_t p = 0; p < logs[log]->machine.size(); ++p) {
        size_t num_used_entries = 0;
        double partial_sum = 0;
        std::vector<log_entry>::const_reverse_iterator iter = 
                                        logs[log]->machine[p].rbegin();
        while (iter != logs[log]->machine[p].rend() &&
            iter->time >= prevtime) {
          if (iter->time < current_time) {
            partial_sum += iter->value;
            ++num_used_entries;
          }
          ++iter;
        }
        sum += partial_sum / num_used_entries;
      }
    }
    // put into the aggregate count
    logs[log]->aggregate.push_back(log_entry(sum, current_time));
    logs[log]->lock.unlock();
  }
}

void distributed_event_logger::periodic_timer() {
  periodic_timer_lock.lock();
  timer ti; ti.start();
  double prevtime = ti.current_time();
  while (!periodic_timer_stop){ 
    collect_instantaneous_log();
    if (ti.current_time() >= prevtime + RECORD_FREQUENCY) {
      prevtime = ti.current_time();
      local_collect_log();
      if (rmi->procid() == 0)  build_aggregate_log();
    }
    periodic_timer_cond.timedwait_ms(periodic_timer_lock, 1000 * TICK_FREQUENCY);
  }
  periodic_timer_lock.unlock();
}

distributed_event_logger::distributed_event_logger():rmi(NULL) {
  pthread_key_create(&key, NULL);
  // clear the bit fields
  has_log_entry.clear();
  thread_local_count_slots.clear();
  periodic_timer_stop = false;
}

void distributed_event_logger::destroy_event_logger() {
  // kill the tick thread
  bool thread_was_started = false;
  periodic_timer_lock.lock();
  // if periodic_timer_stop is false, then
  // thread was started. signal it and wait for it later to 
  // join
  if (periodic_timer_stop == false) {
    periodic_timer_stop = true;
    thread_was_started = true;
    periodic_timer_cond.signal();
  }
  periodic_timer_lock.unlock();
  if (thread_was_started) tick_thread.join();
  // make sure everyone has joined before I start freeing stuff
  rmi->full_barrier();
  delete rmi;
  pthread_key_delete(key);
  // here also free all the allocated memory!
  foreach(uint32_t thr, thread_local_count_slots) {
    if (thread_local_count[thr] != NULL) delete thread_local_count[thr];
  }
  foreach(uint32_t log, has_log_entry) {
    if (logs[log] != NULL) delete logs[log];
  }


}

void distributed_event_logger::set_dc(distributed_control& dc) {
  if (rmi == NULL) {
    rmi = new dc_dist_object<distributed_event_logger>(dc, this);
    // register a deletion callback since the distributed_event_logger
    // will be destroyed only after main

    dc.register_deletion_callback(boost::bind(
                              &distributed_event_logger::destroy_event_logger, 
                              this));

    dc.barrier();
    // everyone starts the timer at the same time
    // at the one distributed synchronization point we have
    ti.start();
    // procid 0 waits 0.2s to skew the local timer a little
    // so everyone else's log has time to show up
    if (rmi->procid() == 0) {
      timer::sleep_ms(200);
    }
    periodic_timer_stop = false;
    // spawn a thread for the tick
    tick_thread.launch(boost::bind(&distributed_event_logger::periodic_timer,
          this));

    // register the metric server callbacks
    add_metric_server_callback("names.json", metric_names_json);
    add_metric_server_callback("metrics_aggregate.json", metric_aggregate_json);
  }
}
    
size_t distributed_event_logger::create_log_entry(std::string name, 
                          log_type::log_type_enum logtype) {
  // look for an entry with the same name
  bool has_existing = false;
  size_t existingid = 0;
  log_entry_lock.lock();
  foreach(uint32_t log, has_log_entry) {
    if (logs[log]->name == name) {
      ASSERT_MSG(logs[log]->is_callback_entry == false,
                 "Cannot convert callback log to counter log");
      has_existing = true;
      existingid = log;
      break;
    }
  }
  log_entry_lock.unlock();
  if (has_existing) return existingid;

  log_group* group = new log_group;
  group->logtype = logtype;
  group->name = name;
  group->callback = NULL;
  group->is_callback_entry = false;
  group->sum_of_instantaneous_entries = 0.0;
  group->count_of_instantaneous_entries = 0;
  // only allocate the machine vector on the root machine.
  // no one else needs it 
  if (rmi->procid() == 0) {
    group->machine.resize(rmi->numprocs());
  } 
  // ok. get an ID
  uint32_t id = allocate_log_entry(group);
  // enforce that all machines are running this at the same time 
  rmi->barrier();
  return id;
}

size_t distributed_event_logger::create_callback_entry(std::string name, 
              boost::function<double(void)> callback) {
  bool has_existing = false;
  size_t existingid = 0;
  log_entry_lock.lock();
  foreach(uint32_t log, has_log_entry) {
    if (logs[log]->name == name) {
      has_existing = true;
      existingid = log;
      break;
    }
  }
  log_entry_lock.unlock();
  if (has_existing) {
    // ok... we have an existing entry. We may
    // overwrite the callback if the callback is NULL
    ASSERT_MSG(logs[existingid]->is_callback_entry == true,
                 "Cannot convert counter log to callback log");

    logs[existingid]->lock.lock();
    ASSERT_MSG(logs[existingid]->callback == NULL, 
        "Cannot create another callback log entry with"
        "the same name %s", name.c_str());
    logs[existingid]->callback = callback;
    logs[existingid]->lock.unlock();
    return existingid;
  }

  log_group* group = new log_group;
  group->logtype = log_type::INSTANTANEOUS;
  group->name = name;
  group->callback = callback;
  group->is_callback_entry = true;
  group->sum_of_instantaneous_entries = 0.0;
  group->count_of_instantaneous_entries = 0;

  // only allocate the machine vector on the root machine.
  // no one else needs it 
  if (rmi->procid() == 0) {
    group->machine.resize(rmi->numprocs());
  } 
  // ok. get an ID
  uint32_t id = allocate_log_entry(group);
  // enforce that all machines are running this at the same time 
  rmi->barrier();
  return id;
}

void distributed_event_logger::thr_inc_log_entry(size_t entry, size_t value) {
  event_log_thread_local_type* ev = get_thread_counter_ref();
  ASSERT_LT(entry, MAX_LOG_SIZE);
  ASSERT_EQ(logs[entry]->is_callback_entry, false);
  ev->values[entry] += value;
}

void distributed_event_logger::thr_dec_log_entry(size_t entry, size_t value) {
  event_log_thread_local_type* ev = get_thread_counter_ref();
  ASSERT_LT(entry, MAX_LOG_SIZE);
  // does not work for cumulative logs
  ASSERT_NE((int)logs[entry]->logtype, (int) log_type::CUMULATIVE);
  ASSERT_EQ(logs[entry]->is_callback_entry, false);
  ev->values[entry] -= value;
}


void distributed_event_logger::free_callback_entry(size_t entry) {
  ASSERT_LT(entry, MAX_LOG_SIZE);
  // does not work for cumulative logs
  logs[entry]->lock.lock();
  ASSERT_EQ(logs[entry]->is_callback_entry, true);
  logs[entry]->callback = NULL;
  logs[entry]->lock.unlock();
}

distributed_event_logger& get_event_log() {
  static distributed_event_logger dist_event_log;
  return dist_event_log;
}





/*
   Used to process the names.json request
*/
std::pair<std::string, std::string> 
static metric_names_json(std::map<std::string, std::string>& vars) {
  std::stringstream strm;
  char *pname = getenv("_");
  std::string progname;
  if (pname != NULL) progname = pname;


  distributed_event_logger& evlog = get_event_log();
  log_group** logs = evlog.get_logs_ptr();
  fixed_dense_bitset<MAX_LOG_SIZE>& has_log_entry = evlog.get_logs_bitset();

  strm << "{\n"
       << "  \"program_name\": \""<< progname << "\",\n"
       << "  \"time\": " << evlog.get_current_time() << ",\n"
       << "  \"metrics\": [\n";
  // output the metrics
  size_t nlogs = has_log_entry.popcount();

  size_t logcount = 0;
  foreach(uint32_t log, has_log_entry) {
    double currate = 0;
    size_t len = logs[log]->aggregate.size();
    if (len >= 1) { 
      double logtime = logs[log]->aggregate.rbegin()->time;
      double logval = logs[log]->aggregate.rbegin()->value;
      double prevtime = 0;
      double prevval = 0;
      if (logs[log]->aggregate.size() >= 2) {
        prevtime = logs[log]->aggregate[len - 2].time;
        prevval = logs[log]->aggregate[len - 2].value;
      }
      currate = (logval - prevval) / (logtime - prevtime);
    }

    strm << "    {\n"
         << "      \"id\":" << log << ",\n"
         << "      \"name\": \"" << logs[log]->name << "\",\n"
         << "      \"cumulative\": " << (int)(logs[log]->logtype) << ",\n"
         << "      \"rate\": " << currate << ",\n"
         << "      \"value\": " << ( logs[log]->aggregate.size() > 0 ?
                                              logs[log]->aggregate.rbegin()->value 
                                              : 0 ) << "\n"
         << "    }\n";
    ++logcount;
    if (logcount < nlogs) {
      strm << ",";
    }
  }
  strm << "  ]\n"
       << "}\n";

  return std::make_pair(std::string("text/plain"), strm.str());
}

std::pair<std::string, std::string> 
static metric_aggregate_json(std::map<std::string, std::string>& vars) {
  double tstart = 0;
  double tend = DBL_MAX;
  bool rate = false;
  std::string name;
  // see what variables there are

  if (vars.count("name")) name = vars["name"];
  if (vars.count("tstart")) tstart = atof(vars["tstart"].c_str());
  if (vars.count("tend")) tend = atof(vars["tend"].c_str());
  if (vars.count("rate")) rate = (atoi(vars["rate"].c_str()) != 0);


  // name is not optional
  name = trim(name);
  if (name.length() == 0) {
    return std::make_pair(std::string("text/plain"), std::string());
  }

  distributed_event_logger& evlog = get_event_log();
  log_group** logs = evlog.get_logs_ptr();
  fixed_dense_bitset<MAX_LOG_SIZE>& has_log_entry = evlog.get_logs_bitset();

  std::stringstream strm;

  foreach(uint32_t log, has_log_entry) {
    if (logs[log]->name == name) {
      strm << "    {\n"
           << "      \"id\":" << log << ",\n"
           << "      \"name\": \"" << logs[log]->name << "\",\n"
           << "      \"cumulative\": " << (int)(logs[log]->logtype) << ",\n"
           << "      \"record\": [";

      std::vector<log_entry> output_entries;
      // annoyingly, json does not let me put a trailing comma in the array.
      // thus I need to first write it to a vector, before dumping it to json

      for (size_t i = 0; i < logs[log]->aggregate.size(); ++i) {
        double logtime = logs[log]->aggregate[i].time;
        double logval = logs[log]->aggregate[i].value;

        if (logtime > tstart && logtime <= tend) {
          // only cumulative logs can have rate
          if (rate == 0 || logs[log]->logtype == log_type::INSTANTANEOUS) {
            output_entries.push_back(log_entry(logval, logtime));
          }
          else {
            double prevval = 0;
            double prevtime = 0;
            if (i > 0) {
              prevtime = logs[log]->aggregate[i - 1].time;
              prevval = logs[log]->aggregate[i - 1].value;
            }
            double currate = 0;
            // avoid divide by zero annoyances
            if (logtime > prevtime) {
              currate = (logval - prevval) / (logtime - prevtime);
            }
            output_entries.push_back(log_entry(currate, logtime));
          }
        }
      }

      for (size_t i = 0 ;i < output_entries.size(); ++i) {
        strm << " [" 
             << output_entries[i].time << ", " 
             << output_entries[i].value 
             << "] ";
        // add a comma if this is not the last entry
        if (i < output_entries.size() - 1) strm << ", ";
      }
      strm << "]\n"
        << "    }\n";
    }
  }
  return std::make_pair(std::string("text/plain"), strm.str());
}

} // namespace graphlab
