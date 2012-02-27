#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/util/timer.hpp>
using namespace graphlab;

#define SEND_LIMIT (64 * 1024 * 1024)
#define SEND_LIMIT_PRINT "64MB"
struct teststruct {

  dc_dist_object<teststruct> rmi;
  teststruct(distributed_control &dc):rmi(dc, this) {
    rmi.barrier();
  }
  
  /**
   *  Receiver
   */
  
  atomic<size_t> ctr;
  void receive_ints(size_t i0, size_t i1, size_t i2, size_t i3) {
    ctr.inc();
  }

  
  void receive_vector(const std::vector<size_t> &s) {
    ctr.inc();
  }

  
  /**
   * Short Sends With Remote Call
   */

  void perform_short_sends_0(size_t number) {
    for (size_t i = 0;i < number; ++i) {
      rmi.remote_call(1, &teststruct::receive_ints, 100,100,1000,5000000);
    }
  }

  void perform_short_pod_sends_0(size_t number) {
    for (size_t i = 0;i < number; ++i) {
      rmi.pod_call(1, &teststruct::receive_ints, 100,100,1000,5000000);
    }
  }


  void perform_long_sends_0(size_t length, size_t number) {
    std::vector<size_t> v(length, 5000000);
    for (size_t i = 0;i < number; ++i) {
      rmi.remote_call(1, &teststruct::receive_vector, v);
    }
  }

  void print_res(double t1, double t2, double t3) {
    std::cout << "Calls Sent in: " << t1 << "s\n";
    std::cout << SEND_LIMIT / t1 / 1024 / 1024 << " MB/s\n";
    std::cout << "Send Completed in: " << t2 << "s\n";
    std::cout << SEND_LIMIT / t2 / 1024 / 1024 << " MB/s\n";
    std::cout << "Receive Completed in: " << t2 << "s\n";
    std::cout << SEND_LIMIT / t3 / 1024 / 1024 << " MB/s\n\n";

  }
  void run_short_sends_0() {
    timer ti;
    std::cout << "Single Threaded " << SEND_LIMIT_PRINT << " sends, 4 integer blocks\n";
    ti.start();
    size_t numsends = SEND_LIMIT / (sizeof(size_t) * 4);
    perform_short_sends_0(numsends);
    double t1 = ti.current_time();
    rmi.dc().flush();
    double t2 = ti.current_time();
    rmi.full_barrier();
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }
  
  
  void run_threaded_short_sends_0(size_t numthreads) {
    timer ti;
    std::cout << numthreads << " threaded " << SEND_LIMIT_PRINT << " sends, 4 integer blocks\n";
    ti.start();
    thread_group thrgrp;
    size_t numsends = SEND_LIMIT / (sizeof(size_t) * 4 * numthreads);
    for (size_t i = 0; i < numthreads; ++i) {
      thrgrp.launch(boost::bind(&teststruct::perform_short_sends_0, this, numsends));
    }
    thrgrp.join();
    double t1 = ti.current_time();
    rmi.dc().flush();
    double t2 = ti.current_time();
    rmi.full_barrier();
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }


  void run_short_pod_sends_0() {
    timer ti;
    std::cout << "Single Threaded  "<< SEND_LIMIT_PRINT <<"  POD sends, 4 integers\n";
    ti.start();
    size_t numsends = SEND_LIMIT / (sizeof(size_t) * 4);
    perform_short_pod_sends_0(numsends);
    double t1 = ti.current_time();
    rmi.dc().flush();
    double t2 = ti.current_time();
    rmi.full_barrier();
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }
  
  
  void run_threaded_short_pod_sends_0(size_t numthreads) {
    timer ti;
    std::cout << numthreads << " threaded "<< SEND_LIMIT_PRINT <<" POD sends, 4 integers\n";
    size_t numsends = SEND_LIMIT / (sizeof(size_t) * 4 * numthreads);
    ti.start();
    thread_group thrgrp;
    for (size_t i = 0; i < numthreads; ++i) {
      thrgrp.launch(boost::bind(&teststruct::perform_short_pod_sends_0, this, numsends));
    }
    thrgrp.join();
    double t1 = ti.current_time();
    rmi.dc().flush();
    double t2 = ti.current_time();
    rmi.full_barrier();
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }




  void run_long_sends_0(size_t length) {
    timer ti;
    std::cout << "Single Threaded " << SEND_LIMIT_PRINT <<" sends, " << sizeof(size_t) * length << " bytes\n";
    ti.start();
    size_t numsends = SEND_LIMIT / (sizeof(size_t) * length);
    perform_long_sends_0(length, numsends);
    double t1 = ti.current_time();
    rmi.dc().flush();
    double t2 = ti.current_time();
    rmi.full_barrier();
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }
  
  
  void run_threaded_long_sends_0(size_t length, size_t numthreads) {
    timer ti;
    std::cout << numthreads << " threaded " << SEND_LIMIT_PRINT <<" sends, " 
                                            << sizeof(size_t) * length << " bytes\n";
    ti.start();
    size_t numsends = SEND_LIMIT / (sizeof(size_t) * length * numthreads);
    perform_long_sends_0(length, numsends);
    thread_group thrgrp;
    for (size_t i = 0; i < numthreads; ++i) {
      thrgrp.launch(boost::bind(&teststruct::perform_long_sends_0, this, length, numsends));
    }
    thrgrp.join();
    double t1 = ti.current_time();
    rmi.dc().flush();
    double t2 = ti.current_time();
    rmi.full_barrier();
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }

};


int main(int argc, char** argv) {
  // init MPI
  mpi_tools::init(argc, argv);
  
  if (mpi_tools::size() != 2) {
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }

  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  distributed_control dc(param);

  teststruct ts(dc);
  if (dc.procid() == 0) {
    ts.run_short_sends_0();
    ts.run_threaded_short_sends_0(2);
    ts.run_threaded_short_sends_0(4);
    ts.run_threaded_short_sends_0(8);
    ts.run_threaded_short_sends_0(16);
    ts.run_short_pod_sends_0();
    ts.run_threaded_short_pod_sends_0(2);
    ts.run_threaded_short_pod_sends_0(4);
    ts.run_threaded_short_pod_sends_0(8);
    ts.run_threaded_short_pod_sends_0(16);
    ts.run_long_sends_0(1024);
    ts.run_threaded_long_sends_0(1024, 2);
    ts.run_threaded_long_sends_0(1024, 4);
    ts.run_threaded_long_sends_0(1024, 8);
    ts.run_threaded_long_sends_0(1024, 16);
    ts.run_long_sends_0(10240);
    ts.run_threaded_long_sends_0(10240, 2);
    ts.run_threaded_long_sends_0(10240, 4);
    ts.run_threaded_long_sends_0(10240, 8);
    ts.run_threaded_long_sends_0(10240, 16);
  }
  else {
    for (size_t tests = 0;tests < 20; ++tests) ts.rmi.full_barrier();
  }
  dc.barrier();
  mpi_tools::finalize();
}