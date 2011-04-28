#ifndef DC_SCTP_COMM_HPP
#define DC_SCTP_COMM_HPP

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
SCTP implementation of the communications subsystem
This is experimental
*/
class dc_sctp_comm:public dc_comm_base {
 public:
   
  dc_sctp_comm();
  
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
  
  ~dc_sctp_comm();
  
  // always true. SCTP can send anywhere
  inline bool channel_active(size_t target) const {
    return true;
  }
    
  /**
    Returns the number of machines in the network.
    Only valid after call to init();
  */
  inline procid_t numprocs() const {
    return nprocs;
  }
  
  inline procid_t procid() const {
    return curid;
  }
  
  inline size_t network_bytes_sent() const {
    return network_bytessent.value;
  }

 inline size_t network_bytes_received() const {
    //TODO
    return 0;
  }

  void flush(size_t target);
  /**
   Sends the string of length len to the target machine dest.
   Only valid after call to init();
   Establishes a connection if necessary
  */
  void send(size_t target, const char* buf, size_t len);
  
  void send2(size_t target, 
             const char* buf1, const size_t len1,
             const char* buf2, const size_t len2); 


 private:
 
  void set_socket_options(int fd);

  // opens the listening sock and spawns a thread to listen on it
  void open_listening();  
  void open_sending();  
  
  // constructs a connection to the target machine
  void connect(size_t target);

  // wrapper around the standard send. but loops till the buffer is all sent
  int sendtosock(int sockfd, size_t target, const char* buf, size_t len, uint16_t stream);
  
  /// all_addrs[i] will contain the IP address of machine i
  std::vector<uint32_t> all_addrs;
  std::vector<struct sockaddr_in> all_sock_addrs;
  std::map<uint32_t, procid_t> addr2id;
  std::vector<uint32_t> portnums;
  
  procid_t curid; 
  procid_t nprocs;
  
  /// the socket we use to listen on (server socket)
  int listensock;
  thread listenthread;
  
  comm_recv_callback_type recvcallback;
  void* tag;
  
  int sendsock;
  
  atomic<size_t> network_bytessent;
  
  void server_handler_loop();
  
  std::vector<dc_receive*> receiver;
  
  std::vector<char> machines_started;
  /// waits for all machines to start up
  void wait_for_all_machines();
  void handle_control(procid_t src, 
                      const char* buf, size_t len);

};

} // namespace dc_impl
} // namespace graphlab
#endif
