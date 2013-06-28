#include <boost/bind.hpp>
#include <graphlab/parallel/fiber_group.hpp>
#include <graphlab/logger/assertions.hpp>
namespace graphlab {

void fiber_group::invoke(const boost::function<void (void)>& spawn_function, 
                         fiber_group* group) {
  spawn_function();
  group->decrement_running_counter();
}


void fiber_group::launch(const boost::function<void (void)> &spawn_function) {
  launch(spawn_function, affinity);
}


void fiber_group::launch(const boost::function<void (void)> &spawn_function,
                         affinity_type worker_affinity) {
  increment_running_counter();
  fiber_control::get_instance().launch(boost::bind(invoke, spawn_function, this), 
                                       stacksize,
                                       worker_affinity);  
}

void fiber_group::join() {
  join_lock.lock();
  // no one else is waiting
  ASSERT_EQ(join_waiting, false);
  // otherwise, we need to wait
  join_waiting = true;
  while(threads_running.value != 0) {
    join_cond.wait(join_lock);
  }
  join_waiting = false;
  join_lock.unlock();
}

} // namespace graphlab

