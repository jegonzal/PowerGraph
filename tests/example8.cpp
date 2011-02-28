#include <iostream>
#include <string>
#include <vector>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/generics/any.hpp>
using namespace graphlab;


void print(any val) {
  val.print(std::cout);
  std::cout << std::endl;
}


int main(int argc, char ** argv) {
  /** Initialization */
  int machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");

  distributed_control dc(machines,"", machineid);
  
  if (dc.procid() == 0) {
    dc.remote_call(1, print, 15);
    dc.remote_call(1, print, 10.5);
    dc.remote_call(1, print, std::string("hello world"));    
  }
    
  int i = dc.procid() == 0 ? 10 : 100;
  dc.broadcast(i, dc.procid() == 0);
  std::cout << i << std::endl;
  assert(i == 10);
}
