/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <signal.h>
#include <sys/time.h>

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
  float lowres_time_seconds() {
    return float(hmstimer.ctr) / 10;
  }


  void my_sleep(size_t sleeplen) {
    fd_set set;
    struct timeval timeout;
    FD_ZERO (&set);
    timeout.tv_sec = sleeplen;
    timeout.tv_usec = 0;
    select (FD_SETSIZE, &set, NULL, NULL, &timeout);
  }
  
  
  /**
   * Precision of deciseconds 
   */
  size_t lowres_time_millis() {
    return hmstimer.ctr * 100;
  }

  
}
