#include <iostream>
#include <cstdio>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/dc_services.hpp>
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
  // init MPI
  mpi_tools::init(argc, argv);
  
  if (mpi_tools::size() != 2) {
    std::cout<< "RPC Example 7: Distributed Object\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }

  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(param);
  
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
  
  std::cout << vec.get(1) << "\n";  
  std::cout << vec.get(2) << "\n";  
  std::cout << vec.get(10) << "\n";
  std::cout << vec.get(11) << std::endl;
  dc.barrier();
  
  mpi_tools::finalize();
}
