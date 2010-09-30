#include <iostream>
#include <string>
#include <map>
#include <graphlab/serialization/podify.hpp>
#include <graphlab/rpc/dc.hpp>
using namespace graphlab;

struct teststruct{
  int a;
  double b;
};
SERIALIZABLE_POD(teststruct);


void print(std::map<int, teststruct> &data1,  
           std::string data2) {
  std::cout << "1.a = " << data1[1].a << std::endl;
  std::cout << "10.b = " << data1[10].b << std::endl;
  std::cout << "string = " << data2 << std::endl;
}



int main(int argc, char ** argv) {
  /** Initialization */
  size_t machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");

  distributed_control dc(machines,"", machineid);
  
  if (dc.procid() == 0) {
    std::map<int, teststruct> data;
    data[1].a = 10;
    data[2].b = 15.0;
    dc.remote_call(1, print, data, "hello world!");
  }
  getchar();
}
