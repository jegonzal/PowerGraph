#include <sys/types.h>
#include <ifaddrs.h>
#include <netinet/in.h>
#include <string.h>

#include <graphlab/logger/assertions.hpp>
#include <sstream>

namespace graphlab {
  
std::string get_local_ip_as_str() {
  std::stringstream strm;
  // code adapted from
  struct ifaddrs * ifAddrStruct = NULL;
  getifaddrs(&ifAddrStruct);
  struct ifaddrs * firstifaddr = ifAddrStruct;
  ASSERT_NE(ifAddrStruct, NULL);
  while (ifAddrStruct != NULL) {
    if (ifAddrStruct->ifa_addr != NULL && 
        ifAddrStruct->ifa_addr->sa_family == AF_INET) {
      char* tmpAddrPtr = NULL;
      // check it is IP4 and not lo0.
      tmpAddrPtr = (char*)&((struct sockaddr_in *)ifAddrStruct->ifa_addr)->sin_addr;
      ASSERT_NE(tmpAddrPtr, NULL);
      if (tmpAddrPtr[0] != 127) {
        strm << int(tmpAddrPtr[0]) << "."
             << int(tmpAddrPtr[1]) << "."
             << int(tmpAddrPtr[2]) << "."
             << int(tmpAddrPtr[3]);
        return strm.str();
      }
      //break;
    }
    ifAddrStruct=ifAddrStruct->ifa_next;
  }
  freeifaddrs(firstifaddr);
  return "";
}

uint32_t get_local_ip() {
  uint32_t ip;
  // code adapted from
  struct ifaddrs * ifAddrStruct = NULL;
  getifaddrs(&ifAddrStruct);
  struct ifaddrs * firstifaddr = ifAddrStruct;
  ASSERT_NE(ifAddrStruct, NULL);
  while (ifAddrStruct != NULL) {
    if (ifAddrStruct->ifa_addr != NULL && 
        ifAddrStruct->ifa_addr->sa_family == AF_INET) {
      char* tmpAddrPtr = NULL;
      // check it is IP4 and not lo0.
      tmpAddrPtr = (char*)&((struct sockaddr_in *)ifAddrStruct->ifa_addr)->sin_addr;
      ASSERT_NE(tmpAddrPtr, NULL);
      if (tmpAddrPtr[0] != 127) {
        memcpy(&ip, tmpAddrPtr, 4);
        break;
      }
      //break;
    }
    ifAddrStruct=ifAddrStruct->ifa_next;
  }
  freeifaddrs(firstifaddr);
  return ip;
}

size_t get_free_tcp_port() {
  int sock = socket(AF_INET, SOCK_STREAM, 0);
  // uninteresting boiler plate. Set the port number and socket type
  sockaddr_in my_addr;
  my_addr.sin_family = AF_INET;
  my_addr.sin_port = 0; // port 0.
  my_addr.sin_addr.s_addr = INADDR_ANY;
  memset(&(my_addr.sin_zero), '\0', 8);
  if (bind(sock, (sockaddr*)&my_addr, sizeof(my_addr)) < 0){
    logger(LOG_FATAL, "Failed to bind to a port 0! Unable to acquire a free TCP port!");
  }
  // get the sock information
  socklen_t slen;
  sockaddr addr;
  slen = sizeof(sockaddr);
  if (getsockname(sock, &addr, &slen) < 0) {
    logger(LOG_FATAL, "Failed to get port information about bound socket");
  }
  size_t freeport = ntohs(((sockaddr_in*)(&addr))->sin_port);
  close(sock);
  return freeport;
}

} // namespace graphlab
