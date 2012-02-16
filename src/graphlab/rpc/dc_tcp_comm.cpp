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


#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <netinet/tcp.h>
#include <ifaddrs.h>
#include <poll.h>

#include <limits>
#include <vector>
#include <string>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>

//#define COMM_DEBUG
namespace graphlab {
 
  namespace dc_impl {
  
    void dc_tcp_comm::init(const std::vector<std::string> &machines,
                           const std::map<std::string,std::string> &initopts,
                           procid_t curmachineid,
                           std::vector<dc_receive*> receiver_){ 

      curid = curmachineid;
      ASSERT_LT(machines.size(), std::numeric_limits<procid_t>::max());
      nprocs = (procid_t)(machines.size());
      receiver = receiver_;
      listenthread = NULL;
      // insert machines into the address map
      all_addrs.resize(nprocs);
      portnums.resize(nprocs);
      // fill all the socks
      socks.resize(nprocs, -1);
      handlers.resize(nprocs, NULL);
      handlerthreads.resize(nprocs, NULL);
      outsocks.resize(nprocs, -1);
      // parse the machines list, and extract the relevant address information
      for (size_t i = 0;i < machines.size(); ++i) {
        // extract the port number
        size_t pos = machines[i].find(":");
        ASSERT_NE(pos, std::string::npos);
        std::string address = machines[i].substr(0, pos);
        size_t port = boost::lexical_cast<size_t>(machines[i].substr(pos+1));
    
        struct hostent* ent = gethostbyname(address.c_str());
        ASSERT_EQ(ent->h_length, 4);
        uint32_t addr = *reinterpret_cast<uint32_t*>(ent->h_addr_list[0]);
    
        all_addrs[i] = addr;
        ASSERT_LT(port, 65536);
        portnums[i] = (uint16_t)(port);
      }
      network_bytessent = 0;
      // if sock handle is set
      std::map<std::string, std::string>::const_iterator iter = 
        initopts.find("__sockhandle__");
      if (iter != initopts.end()) {
        open_listening(atoi(iter->second.c_str()));
      } else {
        open_listening();
      }
    }

    void dc_tcp_comm::close() {
      logstream(LOG_INFO) << "Closing listening socket" << std::endl;
      // close the listening socket
      if (listensock > 0) {
        ::close(listensock);
        listensock = -1;
      }
      // shutdown the listening thread
      // remember that the listening handler is self deleting
      if (listenthread != NULL) {
        listenthread->join();
        delete listenthread;
        listenthread = NULL;
      }
      listenhandler = NULL;
      logstream(LOG_INFO) << "Closing outgoing sockets" << std::endl;
      // close all outgoing sockets
      for (size_t i = 0;i < outsocks.size(); ++i) {
        if (outsocks[i] > 0) {
          ::close(outsocks[i]);
          outsocks[i] = -1;
        }
      }
      logstream(LOG_INFO) << "Closing incoming sockets" << std::endl;
      // close all incoming sockets
      for (size_t i = 0;i < socks.size(); ++i) {
        if (socks[i] > 0) {
          ::close(socks[i]);
          socks[i] = -1;
          // join the receiving threads
          // remember that the receiving handler is self deleting
          if (handlerthreads[i] != NULL) {
            handlerthreads[i]->join();
            delete handlerthreads[i];
          }
          handlerthreads[i] = NULL;
          handlers[i] = NULL;
        }
      }
    }
    void dc_tcp_comm::check_for_out_connection(size_t target) {
      // do we have an outgoing socket to that target?
      // if we don't try to establish a connection
      if (outsocks[target] == -1) {
#ifdef COMM_DEBUG
        logstream(LOG_INFO) << "No existing connection to " << target 
                            << ". Creating now." << std::endl;
#endif
        connect(target);
      }
      ASSERT_NE(outsocks[target], -1);
    }
  
    void dc_tcp_comm::send(size_t target, const char* buf, size_t len) {
      network_bytessent.inc(len);
      check_for_out_connection(target);
#ifdef COMM_DEBUG
      logstream(LOG_INFO) << len << " bytes --> " << target  << std::endl;
#endif

      int err = sendtosock(outsocks[target], buf, len);
      ASSERT_EQ(err, 0);
    }

    void dc_tcp_comm::send2(size_t target, 
                            const char* buf1, const size_t len1,
                            const char* buf2, const size_t len2) {
      network_bytessent.inc(len1 + len2);
      check_for_out_connection(target);
      struct msghdr data;
      struct iovec vec[2];
      vec[0].iov_base = (void*)buf1;
      vec[0].iov_len = len1;
      vec[1].iov_base = (void*)buf2;
      vec[1].iov_len = len2;
 
      data.msg_name = NULL;
      data.msg_namelen = 0;
      data.msg_control = NULL;
      data.msg_controllen = 0;
      data.msg_flags = 0;
      data.msg_iovlen = 2;
      data.msg_iov = vec;
  
 
#ifdef COMM_DEBUG
      logstream(LOG_INFO) << len << " bytes --> " << target  << std::endl;
#endif
      // amount of data to transmit
      size_t dataleft = len1 + len2;
      // while there is still data to be sent
      while(dataleft > 0) {
        size_t ret = sendmsg(outsocks[target], &data, 0);
        // decrement the counter
        dataleft -= ret;
        // restructure the msghdr depending on how much was sent
        struct iovec* newiovecptr = data.msg_iov;
        size_t newiovlen = data.msg_iovlen;
        for (size_t i = 0;i < (size_t)(data.msg_iovlen); ++i) {
          // amount sent was less than this entry.
          // shift the entry and retry
          if (ret < data.msg_iov[i].iov_len) {
            // shift the data
            data.msg_iov[i].iov_len -= ret;
            char* tmp = (char*) data.msg_iov[i].iov_base;
            tmp += ret;
            data.msg_iov[i].iov_base = (void*)tmp;
            break;
          }
          else {
            // amount sent exceeds this entry. we need to 
            // erase this entry (increment the iovec_ptr)
            // and go on to the next entry
            size_t l = std::min(ret, data.msg_iov[i].iov_len);
            newiovlen--;
            newiovecptr++;
            ret -= l;
            if (ret == 0) break;
          }
        }
        data.msg_iov = newiovecptr;
        data.msg_iovlen = newiovlen;
      }
    }

    int dc_tcp_comm::sendtosock(int sockfd, const char* buf, size_t len) {
      size_t numsent = 0;
      while (numsent < len) {
        ssize_t ret = ::send(sockfd, buf + numsent, len - numsent, 0);
        if (ret < 0) {
          logstream(LOG_ERROR) << "send error: " << strerror(errno) << std::endl;
          return errno;
        }
        numsent += ret;
      }
      return 0;
    }
  
    void dc_tcp_comm::set_socket_options(int fd) {
      int flag = 1;
      int result = setsockopt(fd,            /* socket affected */
                              IPPROTO_TCP,     /* set option at TCP level */
                              TCP_NODELAY,     /* name of option */
                              (char *) &flag,  
                              sizeof(int));   
      if (result < 0) {
        logstream(LOG_WARNING) 
          << "Unable to disable Nagle. Performance may be signifantly reduced"
          << std::endl;
      }
    }

    void dc_tcp_comm::flush(size_t target) {
      /*  ASSERT_NE(outsocks[target], -1);
          int one = 1;
          int zero = 0;
          setsockopt(outsocks[target], IPPROTO_TCP, TCP_CORK, &zero, sizeof(zero));
          setsockopt(outsocks[target], IPPROTO_TCP, TCP_CORK, &one, sizeof(one)); */
    }


    void dc_tcp_comm::new_socket(int newsock, sockaddr_in* otheraddr, 
                                 procid_t id) {
      // figure out the address of the incoming connection
      uint32_t addr = *reinterpret_cast<uint32_t*>(&(otheraddr->sin_addr));
      // locate the incoming address in the list
      logstream(LOG_INFO) << "Incoming connection from " 
                          << inet_ntoa(otheraddr->sin_addr) << std::endl;
      ASSERT_LT(id, all_addrs.size());
      ASSERT_EQ(all_addrs[id], addr);
      ASSERT_EQ(socks[id], -1);
      socks[id] = newsock;
      logstream(LOG_INFO) << "Proc " << procid() << " accepted connection "
                          << "from machine " << id << std::endl;
  
      handlers[id] = new socket_handler(*this, newsock, (procid_t)id);
      if (handlerthreads[id] != NULL) delete handlerthreads[id];
      handlerthreads[id] = new thread();
      handlerthreads[id]->launch(boost::bind(&socket_handler::run, 
                                             handlers[id]));
    }



    void dc_tcp_comm::open_listening(int sockhandle) {
      // open listening socket
      if (sockhandle == 0) {
        listensock = socket(AF_INET, SOCK_STREAM, 0);
        // uninteresting boiler plate. Set the port number and socket type
        sockaddr_in my_addr;
        my_addr.sin_family = AF_INET;
        my_addr.sin_port = htons(portnums[curid]);
        my_addr.sin_addr.s_addr = INADDR_ANY;
        memset(&(my_addr.sin_zero), '\0', 8);
        logstream(LOG_INFO) << "Proc " << procid() << " Bind on " 
                            << portnums[curid] << "\n";
        if (bind(listensock, (sockaddr*)&my_addr, sizeof(my_addr)) < 0)
          {
            logstream(LOG_FATAL) << "bind: " << strerror(errno) << "\n";
            ASSERT_TRUE(0);
          }
      }
      else {
        listensock = sockhandle;
      }
      logstream(LOG_INFO) << "Proc " << procid() 
                          << " listening on " << portnums[curid] << "\n";
      ASSERT_EQ(0, listen(listensock, 128));
      // spawn a thread which loops around accept
      listenhandler = new accept_handler(*this, listensock);
      listenthread = new thread();
      listenthread->launch(boost::bind(&accept_handler::run, listenhandler));
    } // end of open_listening


    void dc_tcp_comm::connect(size_t target) {
      if (outsocks[target] != -1) {
        return;
      } else {
        int newsock = socket(AF_INET, SOCK_STREAM, 0);
        set_socket_options(newsock);
        sockaddr_in serv_addr;
        serv_addr.sin_family = AF_INET;
        // set the target port
        serv_addr.sin_port = htons(portnums[target]);
        // set the target address
        serv_addr.sin_addr = *(struct in_addr*)&(all_addrs[target]);
        memset(&(serv_addr.sin_zero), '\0', 8);
        // Connect!
        logstream(LOG_INFO) << "Trying to connect from "
                            << curid << " -> " << target
                            << " on port " << portnums[target] << "\n";
        logger(LOG_INFO, "Destination IP = %s", inet_ntoa(serv_addr.sin_addr));
        // retry 10 times at 1 second intervals
        bool success = false;
        for (size_t i = 0;i < 10; ++i) {
          if (::connect(newsock, (sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
            logstream(LOG_WARNING) 
              << "connect " << curid << " to " << target << ": "
              << strerror(errno) << ". Retrying...\n";
            sleep(1);
            // posix says that 
            /* If connect() fails, the state of the socket is unspecified. 
               Conforming applications should close the file descriptor and 
               create a new socket before attempting to reconnect. */
            ::close(newsock);
            newsock = socket(AF_INET, SOCK_STREAM, 0);
            set_socket_options(newsock);

          } else {
            // send my machine id
            sendtosock(newsock, reinterpret_cast<char*>(&curid), sizeof(curid));
            success = true;
            break;
          }
        }
        if (!success) {
          logstream(LOG_FATAL) << "Failed to establish connection" << std::endl;
        }
        // remember the socket
        outsocks[target] = newsock;
        logstream(LOG_INFO) << "connection from " << curid << " to " << target
                            << " established." << std::endl;

      }
    } // end of connect



    void dc_tcp_comm::socket_handler::run() {
      // get a direct pointer to my receiver
      dc_receive* receiver = owner.receiver[sourceid];
  
      if (receiver->direct_access_support()) {
        // we have direct buffer access!
        size_t buflength;
        char *c = receiver->get_buffer(buflength);
        while(1) {      
          ssize_t msglen = recv(fd, c, buflength, 0);
          // if msglen == 0, the scoket is closed
          if (msglen <= 0) {
            owner.socks[sourceid] = -1;
            // self deleting
            delete this;
            break;
          }
          owner.network_bytesreceived.inc(msglen);
#ifdef COMM_DEBUG
          logstream(LOG_INFO) << msglen << " bytes <-- " 
                              << sourceid  << std::endl;
#endif
          c = receiver->advance_buffer(c, msglen, buflength);
        }
      } else {
        // fall back to using my own buffer
        while(1) {
          char c[10240];
      
          ssize_t msglen = recv(fd, c, 10240, 0);
          // if msglen == 0, the scoket is closed
          if (msglen <= 0) {
            owner.socks[sourceid] = -1;
            // self deleting
            delete this;
            break;
          }
#ifdef COMM_DEBUG
          logstream(LOG_INFO) << msglen << " bytes <-- " 
                              << sourceid  << std::endl;
#endif
          owner.network_bytesreceived.inc(msglen);
          receiver->incoming_data(sourceid, c, msglen);
        }
      }
    } // end of run


    void dc_tcp_comm::accept_handler::run() {
      pollfd pf;
      pf.fd = listensock;
      pf.events = POLLIN;
      pf.revents = 0;

      while(1) {
        // wait for incoming event
        poll(&pf, 1, 1000);
        // if we have a POLLIN, we have an incoming socket request
        if (pf.revents && POLLIN) {
          // accept the socket
          sockaddr_in their_addr;
          socklen_t namelen = sizeof(sockaddr_in);
          int newsock = accept(listensock, (sockaddr*)&their_addr, &namelen);
          if (newsock < 0) {
            break;
          }
          // set the socket options and inform the owner (dc_tcp_comm) that 
          // a socket has been established
          owner.set_socket_options(newsock);
          // before accepting the socket, get the machine number
          procid_t remotemachineid = (procid_t)(-1);
          ssize_t msglen = 0;
          while(msglen != sizeof(procid_t)) {
            msglen += recv(newsock, (char*)(&remotemachineid) + msglen, 
                           sizeof(procid_t) - msglen, 0);
          }
          owner.new_socket(newsock, &their_addr, remotemachineid);
        }
        if (listensock == -1) {
          // the owner has closed
          break;
        }
      }
      logstream(LOG_INFO) << "Listening thread quitting" << std::endl;
      delete this;
    } // end of run


  }; // end of namespace dc_impl
}; // end of namespace graphlab

