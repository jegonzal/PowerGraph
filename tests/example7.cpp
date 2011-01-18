#include <iostream>
#include <cstdio>
#include <graphlab/serialization/podify.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
using namespace graphlab;

template <typename T>
class distributed_vector {
 private:
  dc_dist_object<distributed_vector<T> > rmi;
  std::map<size_t, T> data; 
  
 public:
  distributed_vector(distributed_control &dc):rmi(dc, this) { };
  
  T get(size_t i) {
    procid_t owningmachine = i % rmi.dc().numprocs();
    
    if (owningmachine == rmi.dc().procid()) return data[i];
    else return rmi.remote_request(owningmachine, &distributed_vector<T>::get, i);
  }
  
  
  void set(size_t i, const T& val) {
    procid_t owningmachine = i % rmi.dc().numprocs();
    
    if (owningmachine == rmi.dc().procid()) data[i] = val;
    else rmi.remote_request(owningmachine, &distributed_vector<T>::set, i, val);
  }
};

int main(int argc, char ** argv) {
  /** Initialization */
  global_logger().set_log_level(LOG_INFO);
  size_t machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");

  distributed_control dc(machines,"", machineid, 8, SCTP_COMM);
  dc_services services(dc);
  
  distributed_vector<std::string> vec(dc);
  services.barrier();
  if (dc.procid() == 0) {
    vec.set(10, "set from 0");
    vec.set(11, "set from 0");
  }
  else {
    vec.set(1, "set from 1");
    vec.set(2, "set from 1");
  }
  services.barrier();
  
  std::cout << vec.get(1) << "\n";  
  std::cout << vec.get(2) << "\n";  
  std::cout << vec.get(10) << "\n";
  std::cout << vec.get(11) << std::endl;
  services.barrier();
}
