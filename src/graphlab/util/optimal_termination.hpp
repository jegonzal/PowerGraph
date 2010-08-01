
// TODO: Remove for not being correct
#ifndef GRAPHLAB_OPTIMAL_TERMINATION_HPP
#define GRAPHLAB_OPTIMAL_TERMINATION_HPP

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>

namespace graphlab {
  /**
     Algorithm based on:
     Ho Fung Leung, Hing Fung Ting.
     Optimal Global termination Detection in Shared Memory Asynchronous 
     Multiprocessor Systems
  */
  class optimal_termination {
  public:
    optimal_termination(size_t _ncpus):worker(this),t(&worker) {
      ncpus = _ncpus;
      alpha = (volatile size_t*)malloc(sizeof(size_t) * ncpus);
      beta = (volatile size_t*)malloc(sizeof(size_t) * ncpus);
      for (size_t i = 0;i < ncpus; ++i) {
        alpha[i] = 1;
        beta[i] = 1;
      }
      gamma = 0;
      dead = false;
      // start the thread
      lastwakeup = 0;
      t.start();
    }
  
  
  
    ~optimal_termination(){ 
      free((void*)(alpha));
      free((void*)(beta));
    }
  
    void termination_stage1(size_t cpuid) {
      alpha[cpuid] = 0;
    }

    void cancel_stage1(size_t cpuid) {
      alpha[cpuid] = 1;
    }

    void termination_stage2(size_t cpuid) {
      beta[cpuid] = 0;
    }
    bool poll_termination(size_t cpuid) {
      return dead;
    }
  
    void cancel_stage2(size_t cpuid) {
      beta[cpuid] = 1;
      alpha[cpuid] = 1;
    }
  
    void new_job() {
      size_t a = lastwakeup++;
      lastwakeup = lastwakeup % ncpus;
      a = a % ncpus;
      new_job(a);
    }

    void new_job(size_t cpuhint) {
      while(alpha[cpuhint] == 0);
      gamma = 1;
    }

  
    class optimal_termination_worker: public runnable {      
      optimal_termination* st;
    public:
      optimal_termination_worker(optimal_termination* _st) : st(_st) {  }
      void run() {
        volatile size_t r = 1;
        size_t ncpus = st->ncpus;
        while (r == 1) {
          r = 0;
          for (size_t i = 0;i < ncpus; ++i) {
            r = r || st->beta[i];
          }
          r = r || st->gamma;
          st->gamma = 0;
          usleep(10000);
        }
        st->dead = true;
        return;
      }      
    }; 

  private:
    volatile size_t* alpha;
    volatile size_t* beta;
    volatile size_t gamma;
    size_t ncpus;
    size_t lastwakeup;
    bool dead;
    optimal_termination_worker worker;
    thread t;
  };

}
#endif
