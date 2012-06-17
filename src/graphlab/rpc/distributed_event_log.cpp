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



#include <string>
#include <limits>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/logger/assertions.hpp>
#include <pthread.h>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/macros_def.hpp>
#define EVENT_BAR_WIDTH 40

namespace graphlab {
namespace log_impl {

/// A single entry in time
struct log_entry {
  // The time the log was taken
  double time;
  // The value at the time. If this is a CUMULATIVE entry, this
  // will contain the total number of events since the previous log entry
  size_t value; 
};



enum log_type {
  INSTANTANEOUS = 0, ///< Sum of log values over time are not meaningful 
  CUMULATIVE = 1    ///< Sum of log values over time are meaningful 
};

/// Logging information for a particular log entry (say \#updates)
struct log_group{
  mutex lock;

  /// name of the group
  std::string name;

  /// The type of log. Instantaneous or Cumulative 
  log_type logtype;

  boost::function<size_t(void)> callback;

  /// machine[i] holds a vector of entries from machine i
  std::vector<std::vector<log_entry> > machine;
  /// aggregate holds vector of totals
  std::vector<log_entry> aggregate;
};


const int MAX_LOG_SIZE = 256;
const int MAX_LOG_THREADS = 1024;

/**
 * This is the type that is held in the thread local store
 */
struct event_log_thread_local_type {
  /** The values written to by each thread. 
   * An array with max length MAX_LOG_SIZE 
   */
  size_t* values;
  /** The slot index in the thread_local_count datastructure
   * which holds the storage for the values.
   */
  size_t thlocal_slot;
};


/// The master event log implementation
/// Only one instance of this can be created
class master_event_logger {
  private:
    // a key to allow multiple threads, each to have their
    // own counter. Decreases performance penalty of the
    // the event logger.
    pthread_key_t key;

    dc_dist_object<master_event_logger>* rmi;
    
    // The array of logs. We can only have a maximum of MAX_LOG_SIZE logs
    // This is only created on machine 0
    log_group* logs[MAX_LOG_SIZE];
    // this bit field is used to identify which log entries are active
    fixed_dense_bitset<MAX_LOG_THREADS> has_log_entry;
    mutex log_entry_lock;

    // A collection of slots, one for each thread, to hold 
    // the current thread's active log counter.
    // Threads will write directly into here
    // and a master timer will sum it all up periodically
    size_t* thread_local_count[MAX_LOG_THREADS];
    // a bitset which lets me identify which slots in thread_local_counts
    // are used.
    fixed_dense_bitset<MAX_LOG_THREADS> thread_local_count_slots; 
    mutex thread_local_count_lock;

    // timer managing the frequency at which logs are transmitted to the root
    timer ti; 


    uint32_t allocate_log_entry(log_group* group) {
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
    /**
      * Returns a pointer to the current thread log counter
      * creating one if one does not already exist.
      */
    size_t* get_thread_counter_ref() {
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
        thread_local_count_slots.set_bit(b);
        thread_local_count_lock.unlock();
      }

      event_log_thread_local_type* entry = (event_log_thread_local_type*)(v);
      return entry->values;
    }

    /**
     * Receives the log information from each machine
     */    
    void rpc_collect_log(procid_t srcproc, double srctime,
                         std::vector<size_t> srccounts) {
      
    }

    /** 
     *  Collects the machine level
     *  log entry. and sends it to machine 0
     */
    void local_collect_log() {
      // put together an aggregate of all counters 
      std::vector<size_t> combined_counts(MAX_LOG_SIZE, 0);
      thread_local_count_lock.lock();
      log_entry_lock.lock();

      // for each thread and for each log entry which is 
      // not a callback entry. Accumulate the number of counts
      foreach(uint32_t thr, thread_local_count_slots) {
        size_t *current_thread_counts = thread_local_counts[thr];
        foreach(uint32_t log, has_log_entry) {
          if (logs[log]->callback != NULL) {
            combined_counts[log] += current_thread_counts[log];
          }
        }
      }
      log_entry_lock.unlock();
      thread_local_count_lock.unlock();

      // for each log entry which is a callback entry
      // call the callback to get the counts
      for (size_t i = 0; i < num_log_entries; ++i) {
        if (logs[i]->callback != NULL) {
          combined_counts[i] = logs[i]->callback();
        }
      }
      // send to machine 0
      rmi.remote_call(0, &master_event_logger::rpc_collect_log,
                      rmi.procid(), ti.current_time(), combined_counts);
    }
  public:
    master_event_logger():rmi(NULL) {
      pthread_key_create(&key, NULL);
      next_log_index.value = 0;
      // clear the bit fields
      has_log_entry.clear();
      thread_local_count_slots.clear();
    }

    void set_dc(distributed_control& dc) {
      if (rmi != NULL) {
        rmi = new dc_dist_object<master_event_logger>(dc, this);
        dc.barrer();
        // everyone starts the timer at the same time
        // at the one distributed synchronization point we have
        ti.start();
      }
    }
    /**
     * Creates a new log entry with a given name and log type.
     * Returns the ID of the log. Must be called by 
     * all machines simultaneously with the same settings.
     */
    size_t create_log_entry(std::string name, log_type logtype) {
      log_group* group = new log_group;
      group->logtype = logtype;
      group->name = name;
      group->callback = NULL;
      // only allocate the machine vector on the root machine.
      // no one else needs it 
      if (rmi.procid() == 0) {
       group->machine.resize(rmi->numprocs());
      } 
      // ok. get an ID
      uint32_t id = allocate_log_entry(group);
      // enforce that all machines are running this at the same time 
      rmi.barrier();
      return id;
    }

    /**
     * Creates a new callback log entry with a given name and log type.
     * Returns the ID of the log. Must be called by 
     * all machines simultaneously with the same settings.
     * Callback will be triggered periodically
     */
    size_t create_callback_entry(std::string name, log_type logtype,
                                 boost::function<size_t(void)> callback) {
      log_group* group = new log_group;
      group->logtype = logtype;
      group->name = name;
      group->callback = callback;
      // only allocate the machine vector on the root machine.
      // no one else needs it 
      if (rmi.procid() == 0) {
       group->machine.resize(rmi->numprocs());
      } 
      // ok. get an ID
      uint32_t id = find_free_log_entry(group);
      // enforce that all machines are running this at the same time 
      rmi.barrier();
      return id;

    }

    /**
     * Deletes a log entry created by create_log_entry()
     * or create_callback_entry();
     * Must be called by all machines simultaneously
     */
    void free_log_entry(size_t entry) {
      rmi.barrier();
      // clear the bit
      log_entry_lock.lock();
      ASSERT_TRUE(has_log_entry.get(entry));
      has_log_entry.clear_bit(entry); 
      delete logs[entry];
      log_entry_lock.unlock();
    }
    /**
     * Increments the value of a log entry
     */
    void thr_inc_log_entry(size_t entry, size_t value) {
      size_t* array = get_thread_counter_ref();
      ASSERT_LT(entry, MAX_LOG_SIZE);
      array[entry] += value;
    }

    /**
     * Sets the value of a log entry
     */
    void set_log_entry(size_t entry, size_t value) {
      size_t* array = get_thread_counter_ref();
      ASSERT_LT(entry, MAX_LOG_SIZE);
      array[entry] = value;
    }

};

} // namespace log_impl

} // namespace graphlab
