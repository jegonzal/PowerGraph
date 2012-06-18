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


#ifndef GRAPHLAB_DISTRIBUTED_EVENT_LOG_HPP
#define GRAPHLAB_DISTRIBUTED_EVENT_LOG_HPP
#include <iostream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/dense_bitset.hpp>

namespace graphlab {

// forward declaration because we need this in the
// class but we want dc_dist_object to be able
// to use this class too.
template <typename T>
class dc_dist_object;
class distributed_control;



const int MAX_LOG_SIZE = 256;
const int MAX_LOG_THREADS = 1024;
const double TICK_FREQUENCY = 0.5;
const double RECORD_FREQUENCY = 5.0;



/// A single entry in time
struct log_entry: public IS_POD_TYPE {
  // The value at the time. If this is a CUMULATIVE entry, this
  // will contain the total number of events since the start
  double value;

  // The time the log was taken
  double time;
  log_entry(double value, double time): value(value),time(time) { }
};


namespace log_type {
enum log_type_enum {
  INSTANTANEOUS = 0, ///< Sum of log values over time are not meaningful 
  CUMULATIVE = 1    ///< Sum of log values over time are meaningful 
};
}

/// Logging information for a particular log entry (say \#updates)
struct log_group{
  mutex lock;

  /// name of the group
  std::string name;

  /// Set to true if this is a callback entry
  bool is_callback_entry;

  /// The type of log. Instantaneous or Cumulative 
  log_type::log_type_enum logtype;

  boost::function<double(void)> callback;

  double sum_of_instantaneous_entries;
  size_t count_of_instantaneous_entries;

  /// machine[i] holds a vector of entries from machine i
  std::vector<std::vector<log_entry> > machine;
  /// aggregate holds vector of totals
  std::vector<log_entry> aggregate;
};


/**
 * This is the type that is held in the thread local store
 */
struct event_log_thread_local_type {
  /** The values written to by each thread. 
   * An array with max length MAX_LOG_SIZE 
   */
  double values[MAX_LOG_SIZE];
  size_t thlocal_slot;

  // These are used for time averaging instantaneous values
};


class distributed_event_logger {
  private:
    // a key to allow multiple threads, each to have their
    // own counter. Decreases performance penalty of the
    // the event logger.
    pthread_key_t key;

    dc_dist_object<distributed_event_logger>* rmi;
    
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
    event_log_thread_local_type* thread_local_count[MAX_LOG_THREADS];
    // a bitset which lets me identify which slots in thread_local_counts
    // are used.
    fixed_dense_bitset<MAX_LOG_THREADS> thread_local_count_slots; 
    mutex thread_local_count_lock;

    // timer managing the frequency at which logs are transmitted to the root
    timer ti; 
    thread tick_thread;

    uint32_t allocate_log_entry(log_group* group);
    /**
      * Returns a pointer to the current thread log counter
      * creating one if one does not already exist.
      */
    event_log_thread_local_type* get_thread_counter_ref();

    /**
     * Receives the log information from each machine
     */    
    void rpc_collect_log(size_t srcproc, double srctime,
                         std::vector<double> srccounts);

    void collect_instantaneous_log(); 
    /** 
     *  Collects the machine level
     *  log entry. and sends it to machine 0
     */
    void local_collect_log(); 
    
    // Called only by machine 0 to get the aggregate log
    void build_aggregate_log();

    mutex periodic_timer_lock;
    conditional periodic_timer_cond;
    bool periodic_timer_stop;

    /** a new thread spawns here and sleeps for 5 seconds at a time
     *  when it wakes up it will insert log entries
     */
    void periodic_timer();
  public:
    distributed_event_logger();

    // called by the destruction of distributed_control
    void destroy_event_logger();


    /**
     * Associates the event log with a DC object.
     * Must be called by all machines simultaneously.
     * Can be called more than once, but only the first call will have
     * an effect.
     */
    void set_dc(distributed_control& dc);
    /**
     * Creates a new log entry with a given name and log type.
     * Returns the ID of the log. Must be called by 
     * all machines simultaneously with the same settings.
     */
    size_t create_log_entry(std::string name, log_type::log_type_enum logtype);

    /**
     * Creates a new callback log entry with a given name and log type.
     * Returns the ID of the log. Must be called by 
     * all machines simultaneously with the same settings.
     * Callback will be triggered periodically.
     * Callback entries must be deleted once the callback goes
     * out of scope.
     */
    size_t create_callback_entry(std::string name, 
                                 boost::function<double(void)> callback);

    void free_callback_entry(size_t entry);

    /**
     * Increments the value of a log entry
     */
    void thr_inc_log_entry(size_t entry, size_t value);

    /**
     * Increments the value of a log entry
     */
    void thr_dec_log_entry(size_t entry, size_t value);


};


extern distributed_event_logger& get_event_log();


} // namespace graphlab
#define DECLARE_EVENT(name) size_t name;

#define INITIALIZE_EVENT_LOG(dc) graphlab::get_event_log().set_dc(dc);
#define ADD_CUMULATIVE_EVENT(name, desc) \
    name = graphlab::get_event_log().create_log_entry(desc, log_type::CUMULATIVE);

#define ADD_INSTANTANEOUS_EVENT(name, desc) \
    name = graphlab::get_event_log().create_log_entry(desc, log_type::INSTANTANEOUS);

#define ADD_CALLBACK_EVENT(name, desc, callback) \
    name = graphlab::get_event_log().create_callback_entry(desc, callback);


#define FREE_CALLBACK_EVENT(name) \
  graphlab::get_event_log().free_callback_entry(name);

#define INCREMENT_EVENT(name, count) graphlab::get_event_log().thr_inc_log_entry(name, count);
#define DECREMENT_EVENT(name, count) graphlab::get_event_log().thr_dec_log_entry(name, count);

#endif
