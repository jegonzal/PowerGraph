/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DC_TCP_COMM_HPP
#define DC_TCP_COMM_HPP

#include <sys/socket.h>
#include <netinet/in.h>

#include <vector>
#include <string>
#include <map>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>

namespace graphlab {
namespace dc_impl {
  
/**
 \ingroup rpc_internal
TCP implementation of the communications subsystem.
Provides a single object interface to sending/receiving data streams to
a collection of machines.
*/
class dc_tcp_comm:public dc_comm_base {
 public:
   
  dc_tcp_comm() {}
  
  size_t capabilities() const {
    return COMM_STREAM;
  }
  
  /**
   this fuction should pause until all communication has been set up
   and returns the number of systems in the network.
   After which, all other remaining public functions (numprocs(), send(), etc)
   should operate normally. Every received message should immediate trigger the 
   attached receiver
   
   machines: a vector of strings where each string is of the form [IP]:[portnumber]
   initopts: unused
   curmachineid: The ID of the current machine. machines[curmachineid] will be 
                 the listening address of this machine
   
   recvcallback: A function pointer to the receiving function. This function must be thread-safe
   tag: An additional pointer passed to the receiving function.
  */
  void init(const std::vector<std::string> &machines,
            const std::map<std::string,std::string> &initopts,
            procid_t curmachineid,
            std::vector<dc_receive*> receiver);

  /** shuts down all sockets and cleans up */
  void close();
  
  ~dc_tcp_comm() {
    close();
  }
  
  inline bool channel_active(size_t target) const {
    return (outsocks[target] != -1);
  }

  /**
    Returns the number of machines in the network.
    Only valid after call to init()
  */
  inline procid_t numprocs() const {
    return nprocs;
  }
  
  /**
   * Returns the current machine ID.
   * Only valid after call to init()
   */
  inline procid_t procid() const {
    return curid;
  }
  
  /**
   * Returns the total number of bytes sent
   */
  inline size_t network_bytes_sent() const {
    return network_bytessent.value;
  }

  /**
   * Returns the total number of bytes received
   */
  inline size_t network_bytes_received() const {
    return network_bytesreceived.value;
  }
 
  /// Flushes the TCP stream. \note Not Implemented.
  void flush(size_t target);
  
  /**
   Sends the string of length len to the target machine dest.
   Only valid after call to init();
   Establishes a connection if necessary
  */
  void send(size_t target, const char* buf, size_t len);
  
  /**
   * Sends two buffers one after another to the target machine. 
   * Receiver will receive both buf1 and buf2 consecutively.
   * For instance, buf1 could contain a packet header, while buf2 could
   * contain the packet contents.
   */
  void send2(size_t target, 
             const char* buf1, const size_t len1,
             const char* buf2, const size_t len2); 
  
  
  // receiving socket handler
  class socket_handler {
   public:
    dc_tcp_comm &owner;
    int fd;
    procid_t sourceid;
    socket_handler(dc_tcp_comm& owner, int fd, procid_t id):owner(owner), fd(fd), sourceid(id) {}
    
    void run();
  };
  
  // listening socket handler
  class accept_handler {
   public:
    dc_tcp_comm &owner;
    int listensock;
    
    accept_handler(dc_tcp_comm& owner, int listensock):owner(owner),listensock(listensock) {}
    void run();
  };
  
 private:
  /// Sets TCP_NO_DELAY on the socket passed in fd
  void set_socket_options(int fd);

  /// called when listener receives an incoming socket request
  void new_socket(int newsock, sockaddr_in* otheraddr, procid_t remotemachineid);
  
  /** opens the listening sock and spawns a thread to listen on it.
   * Uses sockhandle if non-zero
   */
  void open_listening(int sockhandle = 0);
  
  
  /// constructs a connection to the target machine
  void connect(size_t target);

  /// wrapper around the standard send. but loops till the buffer is all sent
  int sendtosock(int sockfd, const char* buf, size_t len);
  
  /** checks for the existance of an outgoing connectino to the target
   if none exists, it will create one
     */
  void check_for_out_connection(size_t target);
  
  /// all_addrs[i] will contain the IP address of machine i
  std::vector<uint32_t> all_addrs;
  std::map<uint32_t, procid_t> addr2id;
  std::vector<uint16_t> portnums;

  
  procid_t curid; 
  procid_t nprocs;
  
  /// the socket we use to listen on 
  int listensock;
  accept_handler* listenhandler;
  thread* listenthread;
  
  std::vector<dc_receive*> receiver;
  
  /// socks[i] is the socket to machine i.
  /// There is no socket to the local process ( socks[procid()] is invalid )
  /// If socks[i] == int(-1) then the sock is invalid
  std::vector<int> socks; 
  std::vector<socket_handler*> handlers;
  std::vector<thread*> handlerthreads;
  
  std::vector<int> outsocks; 
  
  atomic<size_t> network_bytessent;
  atomic<size_t> network_bytesreceived;

};

} // namespace dc_impl
} // namespace graphlab
#endif

