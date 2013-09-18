/*  
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
    dc.barrier();
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

  void receive_string(const std::string &s) {
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

  void perform_long_sends_0(size_t length, size_t number) {
    std::vector<size_t> v(length, 5000000);
    for (size_t i = 0;i < number; ++i) {
      rmi.remote_call(1, &teststruct::receive_vector, v);
    }
  }


  void perform_string_sends_0(size_t length, size_t number) {
    std::string s(length, 1);
    for (size_t i = 0;i < number; ++i) {
      rmi.remote_call(1, &teststruct::receive_string, s);
    }
  }



  void print_res(double t1, double t2, double t3) {
    std::cout << "Calls Sent at ";
    std::cout << SEND_LIMIT / t1 / 1024 / 1024 << " MB/s\n";
    std::cout << "Receive Completed at ";
    std::cout << SEND_LIMIT / t3 / 1024 / 1024 << " MB/s\n\n";

  }
  void run_short_sends_0() {
    if (rmi.procid() == 1) {
      rmi.full_barrier();
      return;
    }
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
    if (rmi.procid() == 1) {
      rmi.full_barrier();
      return;
    }
    timer ti;
    std::cout << numthreads << " threaded " << SEND_LIMIT_PRINT << " sends, 4 integer blocks\n";
    ti.start();
    fiber_group thrgrp;
    size_t numsends = SEND_LIMIT / (sizeof(size_t) * 4 * numthreads);
    for (size_t i = 0; i < numthreads; ++i) {
      fiber_control::affinity_type affinity;
      affinity.clear(); affinity.set_bit(i % fiber_control::get_instance().num_workers());
      thrgrp.launch(boost::bind(&teststruct::perform_short_sends_0, this, numsends), affinity);
    }
    thrgrp.join();
    double t1 = ti.current_time();
    rmi.dc().flush();
    double t2 = ti.current_time();
    rmi.full_barrier();
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }



  void run_string_sends_0(size_t length) {
    if (rmi.procid() == 1) {
      rmi.full_barrier();
      return;
    }
    timer ti;
    size_t numsends = SEND_LIMIT / (length);
    std::cout << "Single Threaded " << SEND_LIMIT_PRINT <<" sends, " << length << " bytes * "<< numsends <<  "\n";
    ti.start();
    size_t rd = rdtsc();
    perform_string_sends_0(length, numsends);
    size_t rd2 = rdtsc();
    std::cout << "Completed in: " << ti.current_time() << " seconds\n";
    std::cout << (rd2 - rd) / numsends << " cycles per call\n";
    double t1 = ti.current_time();
    rmi.dc().flush();
    std::cout << "Flush in: " << ti.current_time() << " seconds\n";
    double t2 = ti.current_time();
    rmi.full_barrier();
    std::cout << "Receive Complete in: " << ti.current_time() << " seconds\n";
    double t3 = ti.current_time();
    print_res(t1,t2,t3);
  }


  void run_threaded_string_sends_0(size_t length, size_t numthreads) {
    if (rmi.procid() == 1) {
      rmi.full_barrier();
      return;
    }
    timer ti;
    std::cout << numthreads << " threaded " << SEND_LIMIT_PRINT <<" sends, "
                                            << length << " bytes\n";
    ti.start();
    size_t numsends = SEND_LIMIT / (length * numthreads);
    size_t rd = rdtsc();
    fiber_group thrgrp;
    for (size_t i = 0; i < numthreads; ++i) {
      fiber_control::affinity_type affinity;
      affinity.clear(); affinity.set_bit(i % fiber_control::get_instance().num_workers());
      thrgrp.launch(boost::bind(&teststruct::perform_string_sends_0, this, length, numsends), affinity);
    }
    thrgrp.join();
    size_t rd2 = rdtsc();
    std::cout << (rd2 - rd) / (numthreads * numsends)  << " cycles per call\n";
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
  distributed_control dc;

  if (dc.numprocs() != 2) {
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }
  dc.barrier();
  teststruct ts(dc);
  /*
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
  */
  for (size_t i = 4; i < 24; ++i) {
    ts.run_string_sends_0(1<<i);
  }


  ts.run_threaded_string_sends_0(16, 1);
  ts.run_threaded_string_sends_0(16, 2);
  ts.run_threaded_string_sends_0(16, 4);
  ts.run_threaded_string_sends_0(16, 8);
  ts.run_threaded_string_sends_0(16, 16);
  dc.barrier();
  mpi_tools::finalize();
}
