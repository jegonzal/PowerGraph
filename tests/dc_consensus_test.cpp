#include <iostream>
#include <string>
#include <map>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/async_consensus.hpp>
using namespace graphlab;




class simple_engine_test {
 public:
  dc_dist_object<simple_engine_test> rmi;
  blocking_queue<size_t> queue;
  async_consensus cons;
  atomic<size_t> numactive;;

  simple_engine_test(distributed_control &dc):rmi(dc, this), cons(dc, 4) {
    numactive.value = 4; 
    dc.barrier();
  }

  void add_task_local(size_t i) {
    queue.enqueue(i);
    if (numactive.value < 4) cons.cancel();
  }  
  
  void task(size_t i) {
    if (i < 5) std::cout << "Task " << i << std::endl;
    if (i > 0) {
      if (rmi.numprocs() == 1) {
        add_task_local(i - 1);
      }
      else {
        rmi.remote_call((rmi.procid() + 1) % rmi.numprocs(),
                    &simple_engine_test::add_task_local,
                    i - 1);
      }
    }
  }
  
  bool try_terminate(std::pair<size_t, bool> &job) {
    job.second = false;
    
    numactive.dec();
    cons.begin_done_critical_section();
    job = queue.try_dequeue();
    if (job.second == false) {
      bool ret = cons.end_done_critical_section(true);
      numactive.inc();
      return ret;
    }
    else {
      cons.end_done_critical_section(false);
      numactive.inc();
      return false;
    }
  }
  
  void thread() {
    while(1) {
       std::pair<size_t, bool> job = queue.try_dequeue();
       if (job.second == false) {
          bool ret = try_terminate(job);
          if (ret == true) break;
          if (ret == false && job.second == false) continue;
       }
       task(job.first);
    }
  }
  
  void start_thread() {
    thread_group thrgrp; 
    for (size_t i = 0;i < 4; ++i) {
      launch_in_new_thread(thrgrp, 
                         boost::bind(
                            &simple_engine_test::thread,
                            this), -1);
    }
    
    thrgrp.join();
    ASSERT_EQ(queue.size(), 0);
  }
};


int main(int argc, char ** argv) {
  /** Initialization */
  mpi_tools::init(argc, argv);
  global_logger().set_log_level(LOG_DEBUG);

  dc_init_param param;
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  distributed_control dc(param);
  simple_engine_test test(dc);
  test.add_task_local(1000);
  test.start_thread();
  mpi_tools::finalize();
}
