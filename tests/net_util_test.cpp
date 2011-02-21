#include <sys/types.h>
#include <ifaddrs.h>
#include <netinet/in.h>
#include <string.h>
#include <iostream>
#include <graphlab/util/net_util.hpp>
using namespace graphlab;

int main(int argc, char** argv) {
  std::cout << "local IP: " << get_local_ip() << std::endl;
  std::cout << "local IP str: " << get_local_ip_as_str() << std::endl;
  size_t port = get_free_tcp_port();
  std::cout << "A free port: " << port << std::endl;
  // try to bind on the port   
  int sock = socket(AF_INET, SOCK_STREAM, 0);
  sockaddr_in my_addr;
  my_addr.sin_family = AF_INET;
  my_addr.sin_port = port;
  my_addr.sin_addr.s_addr = INADDR_ANY;
  memset(&(my_addr.sin_zero), '\0', 8);
  if (bind(sock, (sockaddr*)&my_addr, sizeof(my_addr)) < 0){
    std::cout << "Failed to bind to the port" << std::endl;
  }
  else {
    std::cout << "Bind successful." << std::endl;
  }
  close(sock);
}