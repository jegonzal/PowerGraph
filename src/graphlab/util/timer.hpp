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
        (current_time.tv_sec - start_time_.tv_sec) + 
        ((double)(current_time.tv_usec - start_time_.tv_usec))/1.0E6;
       return answer;
    }

    double current_time_millis() const {
      return current_time() * 1000;      
    }

    /**
     * Get the time of day in seconds.  This is useful for
     * initializing random number generators with high precision.
     */
    static double sec_of_day() {
      timeval current_time;
      gettimeofday(&current_time, NULL);
      double answer = 
        current_time.tv_sec + ((double)current_time.tv_usec)/1.0E6;
      std::cout << current_time.tv_sec << std::endl;
      std::cout << current_time.tv_usec << std::endl;
      return answer;
    }

    static size_t usec_of_day() {
      timeval current_time;
      gettimeofday(&current_time, NULL);
      size_t answer = 
        current_time.tv_sec * 1E6 + current_time.tv_usec;
      return answer;
    }


  }; // end of Timer

  //TOOD: Move into timer class namespace (as static functions)
  float lowres_time_seconds();
  //TOOD: Move into timer class namespace (as static functions)
  size_t lowres_time_millis();

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
