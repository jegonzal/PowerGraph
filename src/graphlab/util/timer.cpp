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


#include <iostream>
#include <signal.h>
#include <sys/time.h>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>

std::ostream&  operator<<(std::ostream& out, const graphlab::timer& t) {
  return out << t.current_time();
} 

namespace graphlab {
  void alarm_wakeup(int i);
  
  
  class hundredms_timer {
  public:
    hundredms_timer() {    
      stop = false;
      tout_val.it_interval.tv_sec = 0;
      tout_val.it_interval.tv_usec = 0;
      tout_val.it_value.tv_sec = 0;
      tout_val.it_value.tv_usec = 50000;
      signal(SIGALRM,alarm_wakeup); /* set the Alarm signal capture */    
      setitimer(ITIMER_REAL, &tout_val,0);
      ti.start();
    }
    size_t ctr; 
    timer ti;
    struct itimerval tout_val;
    bool stop;
    
    ~hundredms_timer() {  
      stop = true;
      signal(SIGALRM, SIG_IGN);
    }
  };
  
  hundredms_timer hmstimer;
  
  
  void alarm_wakeup(int i) {
    if (hmstimer.stop) return;
    signal(SIGALRM,alarm_wakeup);
    // compute the desired time till the next 100ms tick by using a real timer call
    double realtime = hmstimer.ti.current_time() ;
    // round down
    hmstimer.ctr = (size_t)(realtime * 10);
    setitimer(ITIMER_REAL, &(hmstimer.tout_val), 0);   
  }

  /**
   * Precision of deciseconds 
   */
  float timer::approx_time_seconds() {
    return float(hmstimer.ctr) / 10;
  }

  /**
   * Precision of deciseconds 
   */
  size_t timer::approx_time_millis() {
    return hmstimer.ctr * 100;
  }


  void timer::sleep(size_t sleeplen) {
    struct timespec timeout;
    timeout.tv_sec = sleeplen;
    timeout.tv_nsec = 0;
    while (nanosleep(&timeout, &timeout) == -1);
  }


  /**
  Sleeps for sleeplen milliseconds.
  */
  void timer::sleep_ms(size_t sleeplen) {
    struct timespec timeout;
    timeout.tv_sec = sleeplen / 1000;
    timeout.tv_nsec = (sleeplen % 1000) * 1000000;
    while (nanosleep(&timeout, &timeout) == -1);
  }
  

  
  
static unsigned long long rtdsc_ticks_per_sec = 0; 
static mutex rtdsc_ticks_per_sec_mutex;

unsigned long long estimate_ticks_per_second() {
  if (rtdsc_ticks_per_sec == 0) {
    rtdsc_ticks_per_sec_mutex.lock();
      if (rtdsc_ticks_per_sec == 0) {
      unsigned long long tstart = rdtsc();
      graphlab::timer::sleep(1);
      unsigned long long tend = rdtsc();
      rtdsc_ticks_per_sec = tend - tstart;
      }
    rtdsc_ticks_per_sec_mutex.unlock();
  }
  return rtdsc_ticks_per_sec;
}

}

