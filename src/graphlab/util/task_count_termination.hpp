#ifndef GRAPHLAB_TASK_COUNT_TERMINATION_HPP
#define GRAPHLAB_TASK_COUNT_TERMINATION_HPP

#include <cassert>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


namespace graphlab {
  /**
   * Shared termination checker based on task counting.  When a
   * processor decides to go to sleep, it should call finish(). If
   * finish() return true, it can end. Otherwise it must proceed.
   * Whenever a new job is created, it must call new_job() after
   * inserting into task queue. When the job is fully completed
   * (i.e. all new jobs created, it must call completed_job().
   */
  class task_count_termination {
    atomic<size_t> newtaskcount;
    atomic<size_t> finishedtaskcount;
    bool force_termination; //signal computation is aborted

  public:
    task_count_termination() : newtaskcount(0), 
                               finishedtaskcount(0), 
                               force_termination(false) { }
    
    ~task_count_termination(){ }
    
    bool finish() {
      return newtaskcount.value == finishedtaskcount.value 
        || force_termination;
    }
    
    void abort(){
      force_termination = true;
    }

    bool is_aborted(){ return force_termination; }

    void restart(){
       force_termination = false;
    }

    void new_job() {
      newtaskcount.inc();
    }
    void completed_job() {
      finishedtaskcount.inc();
      assert(finishedtaskcount.value <= newtaskcount.value);
    }
    
    void print() {
      std::cout << finishedtaskcount.value << " of "
                << newtaskcount.value << std::endl;
    }
  };

}
#endif
