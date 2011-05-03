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
        (current_time.tv_sec - start_time_.tv_sec) + 
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
        current_time.tv_sec + ((double)current_time.tv_usec)/1.0E6;
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
        current_time.tv_sec * 1E6 + current_time.tv_usec;
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


} // end of graphlab namespace

/** 
 * Convenience function. Allows you to call "cout << ti" where ti is
 * a timer object and it will print the number of seconds elapsed
 * since ti.start() was called.
 */
std::ostream&  operator<<(std::ostream& out, const graphlab::timer& t);

/** Returns the current time accurate to 100ms*/

#endif
