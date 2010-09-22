#include <iostream>
#include <distributed/dc.hpp>
using namespace graphlab;


void print(int val) {
  std::cout << val << std::endl;
}

int add_one(int val) {
  return val + 1;
}


int main(int argc, char ** argv) {
  /** Initialization */
  size_t machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");

  distributed_control dc(machines,"", machineid,8,SCTP_COMM);
  
  if (dc.procid() == 0) {
    dc.remote_call(1, print, 10);
    std::cout << "5 plus 1 is : " << dc.remote_request(1, add_one, 5) << std::endl;
    std::cout << "11 plus 1 is : " << dc.remote_request(1, add_one, 11) << std::endl;
  }
  getchar();
}
