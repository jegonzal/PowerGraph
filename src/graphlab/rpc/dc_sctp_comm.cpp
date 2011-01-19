#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <netinet/sctp.h>
#include <ifaddrs.h>
#include <poll.h>

#include <vector>
#include <string>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

#include <graphlab/logger/logger.hpp>
#include <graphlab/rpc/dc_sctp_comm.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>


const uint16_t STREAM_ALL = 0;

const uint16_t STREAM_CONTROL = 1;
//#define COMM_DEBUG
namespace graphlab {
 
namespace dc_impl {
  

dc_sctp_comm::dc_sctp_comm() { }

dc_sctp_comm::~dc_sctp_comm() {
  close();
}
  

void dc_sctp_comm::init(const std::vector<std::string> &machines,
                       const std::map<std::string,std::string> &initopts,
                       procid_t curmachineid,
                       comm_recv_callback_type _recvcallback,
                       void* _tag){
  curid = curmachineid;
  nprocs = machines.size(),
  tag = _tag;
  recvcallback = _recvcallback; 

  // reset the machines_started vector.
  // this is used to implement the low-level barrier
  machines_started.resize(machines.size());
  for (size_t i = 0;i < machines_started.size(); ++i) machines_started[i] = 0;
  
  // insert machines into the address map
  all_addrs.resize(nprocs);
  portnums.resize(nprocs);
  all_sock_addrs.resize(nprocs);
  sendsock = NULL;
  listensock = NULL;
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
    portnums[i] = port;
    
    all_sock_addrs[i].sin_family = AF_INET;
    all_sock_addrs[i].sin_port = htons(port);
    all_sock_addrs[i].sin_addr = *(struct in_addr*)&(all_addrs[i]);
    memset(&(all_sock_addrs[i].sin_zero), '\0', 8);

  }
  
  // set the callback to the "wait for machines call back"
  
  open_listening();
  open_sending();
  wait_for_all_machines();
  
}

void dc_sctp_comm::close() {
  logstream(LOG_INFO) << "Closing listening socket" << std::endl;
  // close the listening socket
  if (listensock > 0) {
    ::close(listensock);
  }
  logstream(LOG_INFO) << "Closing outgoing sockets" << std::endl;
  ::close(sendsock);
}

void dc_sctp_comm::send(size_t target, const char* buf, size_t len) {

  int err = sendtosock(sendsock,target, buf, len, STREAM_ALL);
  ASSERT_EQ(err, 0);
}

void dc_sctp_comm::send2(size_t target, 
                       const char* buf1, const size_t len1,
                       const char* buf2, const size_t len2) {
  send(target,buf1,len1);
  send(target,buf2,len2);
}

int dc_sctp_comm::sendtosock(int sockfd, size_t target, const char* buf, size_t len, uint16_t stream) {
  size_t numsent = 0;
  
  while (numsent < len) {
    
    int ret = sctp_sendmsg(sockfd, buf + numsent, len - numsent,
                           (struct sockaddr*)(&(all_sock_addrs[target])), sizeof(all_sock_addrs[target]),
                           curid,SCTP_ADDR_OVER,stream,0,0);
    #ifdef COMM_DEBUG
    logstream(LOG_INFO) << len << " bytes --> " << target <<"(" << stream << ")" << std::endl;
    #endif

    if (ret < 0) {
      logstream(LOG_ERROR) << "send error: " << strerror(errno) << std::endl;
      return errno;
    }
    numsent += ret;
  }
  return 0;
}
  
void dc_sctp_comm::set_socket_options(int fd) {
  // set events so recvmsg gets all the data it needs from the sndrcvinfo struct
  
   struct sctp_event_subscribe events; 
  bzero(&events, sizeof(events)); 
  events.sctp_data_io_event = 1; 
  setsockopt(fd, IPPROTO_SCTP, 
             SCTP_EVENTS, &events, sizeof(events));

  // set nodelay to minimize latency
    int flag = 1;
   int result = setsockopt(fd,            /* socket affected */
                           IPPROTO_SCTP,     /* set option at TCP level */
                           SCTP_NODELAY,     /* name of option */
                           (char *) &flag,  
                           sizeof(int));   

}

void dc_sctp_comm::flush(size_t target) {
//  ASSERT_NE(outsocks[target], -1);
//  int one = 1;
//  int zero = 0;
//  setsockopt(outsocks[target], IPPROTO_TCP, TCP_CORK, &zero, sizeof(zero));
//  setsockopt(outsocks[target], IPPROTO_TCP, TCP_CORK, &one, sizeof(one));
}


void dc_sctp_comm::open_listening() {
  // open listening socket
  listensock = socket(AF_INET, SOCK_SEQPACKET, IPPROTO_SCTP);
  set_socket_options(listensock);
  // uninteresting boiler plate. Set the port number and socket type
  sockaddr_in my_addr;
  my_addr.sin_family = AF_INET;
  my_addr.sin_port = htons(portnums[curid]);
  my_addr.sin_addr.s_addr = INADDR_ANY;
  memset(&(my_addr.sin_zero), '\0', 8);
  logstream(LOG_INFO) << "Proc " << procid() << " Bind on " << portnums[curid] << "\n";
  if (bind(listensock, (sockaddr*)&my_addr, sizeof(my_addr)) < 0)
  {
    logstream(LOG_FATAL) << "bind: " << strerror(errno) << "\n";
    ASSERT_TRUE(0);
  }
  logstream(LOG_INFO) << "Proc " << procid() << " listening on " << portnums[curid] << "\n";
  ASSERT_EQ(0, listen(listensock, 10));
  // spawn a thread which loops around accept
  listenthread = launch_in_new_thread(boost::bind(&dc_sctp_comm::server_handler_loop, 
                                                 this));
}

void dc_sctp_comm::open_sending() {
  sendsock = socket(AF_INET, SOCK_SEQPACKET, IPPROTO_SCTP);
  set_socket_options(sendsock);
}


void dc_sctp_comm::server_handler_loop() {
  sctp_sndrcvinfo info;
  int flags;
  while(1) {
    char c[10240];
    flags = 0;
    int msglen = sctp_recvmsg(listensock, c, 10240, NULL, NULL, &info, &flags);
    // srcid is set in info.sinfo_ppid
    size_t sourceid = info.sinfo_ppid;
    uint16_t stream = info.sinfo_stream;
    
    // if msglen == 0, the scoket is closed
    if (msglen == 0) {
      break;
    }
    else if (msglen < 0) {
      logstream(LOG_INFO) << "bind: " << strerror(errno) << "\n";
      sched_yield();
      continue;
    }
    #ifdef COMM_DEBUG
    logstream(LOG_INFO) << msglen << " bytes <-- " << sourceid  << "(" << stream << ","<< flags << ")" << std::endl;
    #endif
    if (stream == STREAM_ALL) recvcallback(tag, sourceid, c, msglen);
    else if (stream == STREAM_CONTROL) handle_control(sourceid, c, msglen);
    else logstream(LOG_FATAL) << "Unexpected stream number\n";
  }
}



void dc_sctp_comm::handle_control(procid_t src, 
                                  const char* buf, size_t len) {
  machines_started[src] = 1;
}


/**
 * A low-level barrier implemented directly on top of the communication layer
 * Used to synchronize and ensure all machines are alive. This is necessary
 * because SCTP is connection-less.
 * 
 * \warning This implementation is still not quite correct since I am making
 * the assumption that the packet-header is too small to fragment.
 * (which is probably is)
 */
void dc_sctp_comm::wait_for_all_machines() {
  logger(LOG_INFO, "Waiting for machines...");
  // all machines will repeatedly try to to signal machine 0
  // create an empty packet
  packet_hdr pack;
  pack.len = 0;
  pack.src = procid();
  pack.packet_type_mask = 0;
  
  if (procid() == 0) {
    machines_started[0] = 1;
    bool allstarted = false;
    
    while(allstarted == false) {
      allstarted = true;
      for (size_t i = 0;i < machines_started.size(); ++i) {
        volatile char &c = machines_started[i];
        if (c == 0) {
          allstarted = false;
          break;
        }
      }
      if (allstarted == false) usleep(100000);
    }
    // ok. all have started!
    // send a wake up to everyone
    for (size_t i = 1; i < numprocs(); ++i) {
      sendtosock(sendsock, i, reinterpret_cast<char*>(&pack) ,sizeof(packet_hdr), STREAM_CONTROL);
    }
  }
  else {
    volatile char &c0 = machines_started[0];
    while(c0 == 0) {
      sendtosock(sendsock,0, reinterpret_cast<char*>(&pack) ,sizeof(packet_hdr), STREAM_CONTROL);
      usleep(100000);
    }
  }
}




} // dc_impl
} // graphlab
