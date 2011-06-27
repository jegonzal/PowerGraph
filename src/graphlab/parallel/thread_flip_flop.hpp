#ifndef GRAPHLAB_PARALLEL_THREAD_FLIP_FLOP_HPP
#define GRAPHLAB_PARALLEL_THREAD_FLIP_FLOP_HPP

#include <graphlab/parallel/pthread_tools.hpp>

namespace graphlab {
  /**
  
    \ingroup util
    Manages two groups of threads such that only
    one group is active and running at any one point in time.
    Provides a multi-threaded version of co-routines.
    
    Group 0 is the first thread group to run. All threads in
    group 1 calls wait() and blocks. When threads in group 0 
    calls wait(), they will block as well, but when all threads in
    group 0 enter wait(), all threads in group 1 will be woken up.

    This then repeats, swapping group 0 and group 1.
  */
  class thread_flip_flop {
   private:
    
    mutex mut;
    conditional non_running_cond[2];
    size_t running_group;
    cancellable_barrier* bar[2];
    bool init; 
    bool alive;
   public:
    /**
      Initializes the flip flop object and provides the number
      of threads in each group as well as optionally, the
      first active group (defaults to 0).
    */
    inline thread_flip_flop(size_t group_0_size, 
                     size_t group_1_size, 
                     size_t first_running_group = 0) {
      running_group = first_running_group;
      alive = true;
      bar[0] = new barrier(group_0_size);
      bar[1] = new barrier(group_0_size);
    }
                     
    inline void wait(size_t group) {
      if (!alive) return;
      ASSERT_TRUE(group == 0 || group == 1);
      
      // if I am the current running group, 
      // lets wait for everyone in the group
      // to stop running with a barrier
      if (running_group == group) {
        bar[group]->wait();
        // exit if no longer alive
        if (!alive) return;
        // now every one in the running group is here
        // I need to flip the running and the not running group
        mut.lock();
        // only one person has to do the flip and the broadcast
        if (running_group == group) {
          running_group = !group;
          non_running_cond[running_group].broadcast();
        }
        mut.unlock();
      }
      
      // If I am not the running group, we wait
      mut.lock();
      while (running_group != group && alive) {
        non_running_cond[group].wait(mut);
      }
      mut.unlock();
    }
    
    inline void stop_blocking() {
      alive = false;
      non_running_cond[0].broadcast();
      non_running_cond[1].broadcast();
      bar[0]->cancel();
      bar[1]->cancel();
    }
    
    inline ~thread_flip_flop() {
      delete bar[0];
      delete bar[1];
    }
  };
}

#endif
