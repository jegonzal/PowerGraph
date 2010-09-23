#include <iostream>
#include <string>
#include <vector>
#include <graphlab/rpc/dc.hpp>
using namespace graphlab;


void print(std::string val) {
  std::cout << val << std::endl;
}

std::vector<int> add_one(std::vector<int> val) {
  val.push_back(1);
  return val;
}


int main(int argc, char ** argv) {
  /** Initialization */
  size_t machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");

  distributed_control dc(machines,"", machineid);
  
  if (dc.procid() == 0) {
    dc.remote_call(1, print, "hello world!");
    /** Create a vector with a few elements */  
    std::vector<int> vec;
    vec.push_back(1); vec.push_back(2);
    /** Call the remote machine */  
    vec = dc.remote_request(1, add_one, vec);
    
    /** Print the vector */  
    for (size_t i = 0; i < vec.size(); ++i) {
      std::cout << vec[i] << ", ";
    }
    std::cout << std::endl;
  }
  getchar();
}
