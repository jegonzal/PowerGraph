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


#include <iostream>
#include <cstdio>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
using namespace graphlab;

template <typename T>
class distributed_vector {
 private:
  dc_dist_object<distributed_vector<T> > rmi; // The local RMI object
  std::map<size_t, T> data;   // storage
  mutex lock;   // protect the storage
 public:
  distributed_vector(distributed_control &dc):rmi(dc, this) { };
  
  ///Reads the value at key i
  T get(size_t i) {
    // find the owning machine
    procid_t owningmachine = i % rmi.dc().numprocs();
    
    if (owningmachine == rmi.dc().procid()) {
      // if I own the data. just read and return it
      T ret;
      lock.lock();
      ret = data[i];
      lock.unlock();
      return ret;
    }
    else {
      // otherwise I need to go to another machine
      return rmi.remote_request(owningmachine, 
                                &distributed_vector<T>::get, 
                                i);
    }
  }
  
  /// Sets the value at key i
  void set(size_t i, const T& val) {
    // find the owning machine
    procid_t owningmachine = i % rmi.dc().numprocs();
    
    if (owningmachine == rmi.dc().procid()) {
      // if I own the data set it
      lock.lock();
      data[i] = val;
      lock.unlock();
    }
    else {
      // forward the write to another machine
      rmi.remote_request(owningmachine, 
                         &distributed_vector<T>::set, 
                         i, 
                         val);
    }
  }
};

int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);
  distributed_control dc;

  if (dc.numprocs() != 2) {
    std::cout<< "RPC Example 7: Distributed Object\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }
  
  size_t i = 10;
  dc.all_reduce(i);
  std::cout << i << "\n";
  // create a distributed vector
  distributed_vector<std::string> vec(dc);
  dc.barrier();
  if (dc.procid() == 0) {
    vec.set(10, "set from 0");
    vec.set(11, "set from 0");
  }
  else {
    vec.set(1, "set from 1");
    vec.set(2, "set from 1");
  }
  dc.barrier();
  if (dc.procid() == 0) {
    std::cout << vec.get(1) << "\n";  
    std::cout << vec.get(2) << "\n";  
    std::cout << vec.get(10) << "\n";
    std::cout << vec.get(11) << std::endl;
  }
  dc.barrier();
  if (dc.procid() == 1) {
    std::cout << vec.get(1) << "\n";  
    std::cout << vec.get(2) << "\n";  
    std::cout << vec.get(10) << "\n";
    std::cout << vec.get(11) << std::endl;
  }
  dc.barrier();
  
  mpi_tools::finalize();
}
