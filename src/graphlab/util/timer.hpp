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


#ifndef GRAPHLAB_TIMER_HPP
#define GRAPHLAB_TIMER_HPP

#include <sys/time.h>
#include <stdio.h>

#include <iostream>

namespace graphlab {
  /**
   * \ingroup util  
   * A simple class that can be used for
   * benchmarking/timing up to microsecond resolution.
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
        (double)(current_time.tv_sec - start_time_.tv_sec) + 
        ((double)(current_time.tv_usec - start_time_.tv_usec))/1.0E6;
       return answer;
    }

    /**
     Like current_time but returns in milliseconds
    */
    double current_time_millis() const {
      return current_time() * 1000;      
    }

    /**
     * Get the number of seconds (as a floating point value)
     * since the Unix Epoch.
     */
    static double sec_of_day() {
      timeval current_time;
      gettimeofday(&current_time, NULL);
      double answer = 
        (double)current_time.tv_sec + ((double)current_time.tv_usec)/1.0E6;
//      std::cout << current_time.tv_sec << std::endl;
//      std::cout << current_time.tv_usec << std::endl;
      return answer;
    }

    /**
     * Returns only the micro-second component of the 
     * time since the Unix Epoch.
     */
    static size_t usec_of_day() {
      timeval current_time;
      gettimeofday(&current_time, NULL);
      size_t answer = 
        (size_t)current_time.tv_sec * 1000000 + (size_t)current_time.tv_usec;
      return answer;
    }


  }; // end of Timer
  
  /**
   Returns the time since program start.
   This value is only updated once every 100ms.
  */
  float lowres_time_seconds();
  
  /**
   Returns the number of milliseconds time since program start.
   This value is only updated once every 100ms.
  */
  size_t lowres_time_millis();
  
  /**
  Sleeps for sleeplen seconds
  */
  void my_sleep(size_t sleeplen);

  /**
  Sleeps for sleeplen milliseconds.
  */
  void my_sleep_ms(size_t sleeplen);

} // end of graphlab namespace

/** 
 * Convenience function. Allows you to call "cout << ti" where ti is
 * a timer object and it will print the number of seconds elapsed
 * since ti.start() was called.
 */
std::ostream&  operator<<(std::ostream& out, const graphlab::timer& t);


#if defined(__i386__)
static inline unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}
#elif defined(__x86_64__)
static inline unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#else
static inline unsigned long long rdtsc(void) {
  return 0;
}
#endif

#endif

