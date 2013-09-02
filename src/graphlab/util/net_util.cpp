/**  
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


#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <sys/types.h>
#include <ifaddrs.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/net_util.hpp>
namespace graphlab {
 

bool str_to_ip(const char* c, uint32_t& out) {
  if (c == NULL) return false;
  else return inet_pton(AF_INET, c, &out) > 0;
}

bool ip_to_str(uint32_t ip, std::string& out) {
  char ipstring[INET_ADDRSTRLEN] = {0};
  const char* ret = inet_ntop(AF_INET, &ip, ipstring, INET_ADDRSTRLEN);
  if (ret == NULL) return false;
  out = std::string(ipstring);
  return true;
}



std::string get_local_ip_as_str(bool print) {
  uint32_t ip = get_local_ip(print);
  if (ip == 0) return "127.0.0.1";
  else {
    std::string out;
    bool ip_conversion_success = ip_to_str(ip, out);
    ASSERT_TRUE(ip_conversion_success);
    return out;
  }
}

uint32_t get_local_ip(bool print) {
  // see if GRAPHLAB_SUBNET environment variable is set
  char* c_subnet_id = getenv("GRAPHLAB_SUBNET_ID");
  char* c_subnet_mask = getenv("GRAPHLAB_SUBNET_MASK");
  uint32_t subnet_id = 0;
  uint32_t subnet_mask = 0;
  std::string str_subnet_id, str_subnet_mask;
  // try to convert to a valid address when possible
  if (c_subnet_id != NULL) {
    if (!str_to_ip(c_subnet_id, subnet_id)) {
      std::cout << "Unable to convert GRAPHLAB_SUBNET_ID to a valid address. Cannot continue\n";
      exit(1); 
    }
  }
  if (c_subnet_mask != NULL) {
    if (!str_to_ip(c_subnet_mask, subnet_mask)) {
      std::cout << "Unable to convert GRAPHLAB_SUBNET_MASK to a valid address. Cannot continue\n";
      exit(1); 
    }
  }

  // error checking. 
  // By the end of this block, we should either have both subnet_id and subnet_mask filled
  // to reasonable values, or are dead.
  
  if (c_subnet_id == NULL && c_subnet_mask != NULL) {
    // If subnet mask specified but not subnet ID, we cannot continue.
    std::cout << "GRAPHLAB_SUBNET_MASK specified, but GRAPHLAB_SUBNET_ID not specified.\n";
    std::cout << "We cannot continue\n";
    exit(1);
  } 
  if (c_subnet_id != NULL && c_subnet_mask == NULL) {
    if (print) {
      std::cout << "GRAPHLAB_SUBNET_ID specified, but GRAPHLAB_SUBNET_MASK not specified.\n";
      std::cout << "We will try to guess a subnet mask\n";
    }
    // if subnet id specified, but not subnet mask. We can try to guess a mask 
    // by finding the first "on" bit in the subnet id, and matching everything
    // to the left of it.
    // easiest way to do that is to left extend the subnet_id
    subnet_mask = subnet_id;
    subnet_mask = ntohl(subnet_mask);
    subnet_mask = subnet_mask | (subnet_mask << 1);
    subnet_mask = subnet_mask | (subnet_mask << 2);
    subnet_mask = subnet_mask | (subnet_mask << 4);
    subnet_mask = subnet_mask | (subnet_mask << 8);
    subnet_mask = subnet_mask | (subnet_mask << 16);
    subnet_mask = htonl(subnet_mask);
  }
  else {
    if (print) {
      std::cout << "GRAPHLAB_SUBNET_ID/GRAPHLAB_SUBNET_MASK environment variables not defined.\n";
      std::cout << "Using default values\n";
    }
  }
  ip_to_str(subnet_id, str_subnet_id);
  ip_to_str(subnet_mask, str_subnet_mask);
  
  // make sure this is a valid subnet address.
  if (print) {
      std::cout << "Subnet ID: " << str_subnet_id << "\n";
      std::cout << "Subnet Mask: " << str_subnet_mask << "\n";
      std::cout << "Will find first IPv4 non-loopback address matching the subnet" << std::endl;
  }
  uint32_t ip(0);
  // code adapted from
  struct ifaddrs * ifAddrStruct = NULL;
  getifaddrs(&ifAddrStruct);
  struct ifaddrs * firstifaddr = ifAddrStruct;
  ASSERT_NE(ifAddrStruct, NULL);
  bool success = false;
  while (ifAddrStruct != NULL) {
    if (ifAddrStruct->ifa_addr != NULL && 
        ifAddrStruct->ifa_addr->sa_family == AF_INET) {
      char* tmpAddrPtr = NULL;
      // check it is IP4 and not lo0.
      tmpAddrPtr = (char*)&((struct sockaddr_in *)ifAddrStruct->ifa_addr)->sin_addr;
      ASSERT_NE(tmpAddrPtr, NULL);
      if (tmpAddrPtr[0] != 127) {
        memcpy(&ip, tmpAddrPtr, 4);
        // test if it matches the subnet
        if ((ip & subnet_mask) == subnet_id) {
          success = true;
          break;
        }
      }
      //break;
    }
    ifAddrStruct=ifAddrStruct->ifa_next;
  }
  freeifaddrs(firstifaddr);
  if (!success) {
    // if subnet addresses specified, and if we cannot find a valid network. Fail."
    if (c_subnet_id!= NULL) {
      std::cout << "Unable to find a network matching the requested subnet\n";
      exit(1);
    } else {
      std::cout << "Unable to find any valid IPv4 address. Defaulting to loopback\n";
    }
  }
  return ip;
}

std::pair<size_t, int> get_free_tcp_port() {
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
  std::pair<size_t, int> ret(freeport, sock);
  return ret;
}

} // namespace graphlab

