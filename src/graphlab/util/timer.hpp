#ifndef GRAPHLAB_TIMER_HPP
#define GRAPHLAB_TIMER_HPP

#include <sys/time.h>
#include <stdio.h>

#include <iostream>

namespace graphlab {
  /**
   *   \class timer A simple class that can be used for
   *   benchmarking/timing up to microsecond resolution.
   */
  class timer {
  private:
    timeval start_time_;   
  public:
    timer() { }
    
    //! Starts the timer. 
    void start() { gettimeofday(&start_time_, NULL); }
    
    /** 
     * Returns the number of seconds since start() was called Behavior
     *  is undefined if start() was not called before calling
     *  current_time()
     */
    double current_time() const {
      timeval current_time;
      gettimeofday(&current_time, NULL);
      double answer = 
       // (current_time.tv_sec + ((double)current_time.tv_usec)/1.0E6) -
       // (start_time_.tv_sec + ((double)start_time_.tv_usec)/1.0E6);
         current_time.tv_sec-start_time_.tv_sec+ ((double)(current_time.tv_usec-start_time_.tv_usec))/1.0E6;
       return answer;
    }

    double current_time_millis() const {
      return current_time() * 1000;      
    }

  }; // end of Timer
  
  float lowres_time_seconds();
  size_t lowres_time_millis();

} // end of graphlab namespace

/** 
 * Convenience function. Allows you to call "cout << ti" where ti is
 * a timer object and it will print the number of seconds elapsed
 * since ti.start() was called.
 */
std::ostream&  operator<<(std::ostream& out, const graphlab::timer& t);

/** Returns the current time accurate to 100ms*/

#endif
