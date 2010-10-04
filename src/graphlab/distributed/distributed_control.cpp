#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/metrics/distributed_metrics.hpp>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <netinet/tcp.h>

#include <ifaddrs.h>
#include <limits>
namespace graphlab {

distributed_control *dc_singleton_ptr = NULL;

distributed_control::distributed_control(int *pargc, char*** pargv) {
  assert(dc_singleton_ptr == NULL);
  msgsent.value = 0;
  msgprocessed.value = 0;
  done = 0;
  int provided;
  MPI_Init_thread(pargc, pargv, MPI_THREAD_MULTIPLE, &provided);
  if (provided >= MPI_THREAD_MULTIPLE) {
    logger(LOG_INFO, "MPI threading available");
  }
  else {
    logger(LOG_FATAL, "MPI threading not available. Threading level %d < %d", 
                      provided, MPI_THREAD_MULTIPLE);
    ASSERT_TRUE(provided);
  }

  // create the thread specific data



  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  // get the number of processes
  int i_nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &i_nprocs);
  nprocs = (procid_t)i_nprocs;
  
  logger(LOG_INFO, "%d Processes", nprocs);
  int i_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_id);
  id = (procid_t)i_id;

  // initialize the local UDT socket
  // Now to pick up the IP address of everyone
  sync_ip_list();

  connect_udt();
  send_thread = new background_send_thread(*this);
  send_thread->start();
  dc_singleton_ptr = this;
  terminator = new distributed_terminator(*this);
  terminator->set_use_control_packets(true);
  mpi_barrier();
  
}




distributed_control::~distributed_control() {

   logger(LOG_INFO, "Shutting down distributed control");
  close_all_connections();
  // close all sockets


  mpi_barrier();
  // kill the send thread
  send_requests.stop_blocking();
  send_thread->join();
  logstream(LOG_INFO) << "Total Bytes Transmitted: " << send_thread->bytes_sent << std::endl;
  delete send_thread;

  
  // kill the message handling thread
  dispatch_requests.stop_blocking();
  for (size_t i = 0; i < dispatch_thread.size(); ++i) {
    dispatch_thread[i]->join();
    delete dispatch_thread[i];
  }
  
  delete [] all_addrs;
  // delete buffers
  for (size_t i = 0; i <buffer.size(); ++i) {
    if (buffer[i].buffer != NULL) free(buffer[i].buffer);
  }
  logger(LOG_INFO, "MPI_Finalize");
  MPI_Finalize();
}


void distributed_control::init_message_processing(size_t nummsgthreads) {
  ASSERT_GT(nummsgthreads, 0);
  create_receive_buffers(131072);
  // begin message processing threads
  dispatch_thread.resize(nummsgthreads);
  for (size_t i = 0; i < nummsgthreads; ++i) {
    dispatch_thread[i] = new message_dispatch_thread(*this);
    dispatch_thread[i]->start();
  }
  
  // begin socket receive threads 
  size_t numrecvthreads = socks.size() - 1;
  procs.resize(numrecvthreads );
  for (size_t i = 0;i < numrecvthreads ; ++i) {
    procs[i].done = &done;
    procs[i].dc = this;
    procthreads.launch(&(procs[i]));
  }
  mpi_barrier();
}

void distributed_control::set_socket_options(int fd) {
  int flag = 1;
  int result = setsockopt(fd,            /* socket affected */
                          IPPROTO_TCP,     /* set option at TCP level */
                          TCP_NODELAY,     /* name of option */
                          (char *) &flag,  
                          sizeof(int));   
  if (result < 0) {
    logger(LOG_WARNING, "Unable to disable Nable. Performance may be signifantly reduced");
  }
}

void distributed_control::report_stats() {
     distributed_metrics::instance(this)->set_value("total_bytes", (double) send_thread->bytes_sent);
     barrier();
     
    if (id == 0) distributed_metrics::instance(this)->report();
}

size_t distributed_control::read_tcp_buflen(int fd) {
  int rcvbuflen;
  socklen_t size = sizeof(int);
  int ret = getsockopt(fd, 0,SO_RCVBUF, (void*)(&rcvbuflen), &size);
  if (ret < 0) {
    logstream(LOG_FATAL) << "Getsockopt failure: " << strerror(errno) << std::endl;
    ASSERT_FALSE(true);
  }
  return rcvbuflen;
}

void distributed_control::get_local_ip(char ip[4]) {
  // code adapted from
  // http://stackoverflow.com/questions/212528/linux-c-get-the-ip-address-of-local-computer
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
        memcpy(ip, tmpAddrPtr, 4);
        break;
      }
      //break;
    }
    ifAddrStruct=ifAddrStruct->ifa_next;
  }
  freeifaddrs(firstifaddr);
}

void distributed_control::sync_ip_list() {
  // get the local IP
  char localip[4];
  get_local_ip(localip);
  struct in_addr* localipinaddr = (struct in_addr*)(localip);
  logger(LOG_INFO, "Local IP = %s", inet_ntoa(*localipinaddr));

  // build the gather op for all addresses
  all_addrs = new uint32_t[numprocs()];

  MPI_Allgather(localip, 4,  MPI_CHAR,
                all_addrs , 4, MPI_CHAR, MPI_COMM_WORLD);
  
  for (procid_t i = 0;i < numprocs(); ++i) {
    logger(LOG_INFO, "Node %d IP = %s",i, inet_ntoa(*(struct in_addr*)&(all_addrs[i])));
  }
  localport = DC_LOCAL_BASE_PORT_NUM + procid();
}


void distributed_control::connect_udt() {
  open_listening_udt();
  socks.resize(numprocs());
  socks[procid()] = 0;
  mpi_barrier();
  // this is a little tricky. both connect() and accept() blocks.
  // so I have to do this in phases
  for (procid_t i = 0;i < numprocs(); ++i) {
    if (procid() == i) {
      // my turn to connect to everyone. I only need to connect to those with
      // id above me
      for (procid_t j = i+1;j < numprocs(); ++j) {
      
        sockfd_t newsock = socket(AF_INET, SOCK_STREAM, 0);
        sockaddr_in serv_addr;
        serv_addr.sin_family = AF_INET;
        // set the target port
        serv_addr.sin_port = htons(DC_LOCAL_BASE_PORT_NUM + j);
        // set the target address
        serv_addr.sin_addr = *(struct in_addr*)&(all_addrs[j]);
        memset(&(serv_addr.sin_zero), '\0', 8);
        // Connect!
        logstream(LOG_INFO) << "Trying to connect from "
                            << i << " -> " << j
                            << " on port " << DC_LOCAL_BASE_PORT_NUM + j << "\n";
        logger(LOG_INFO, "Destination IP = %s", inet_ntoa(serv_addr.sin_addr));
        if (connect(newsock, (sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
          logstream(LOG_FATAL) << "connect " << i << " to " << j << ": "
                               << strerror(errno) << "\n";
          ASSERT_TRUE(false);
        }
        logstream(LOG_INFO) << "Sent connection " << i << " -> " << j << "\n";
        // store in the socks table
        set_socket_options(newsock);
        socks[j] = newsock;

      }
    }
    else if (procid() > i) {
      logstream(LOG_INFO) << "Proc " << procid() << " Waiting to accept\n";
      
      sockaddr_in their_addr;
      socklen_t namelen = sizeof(sockaddr_in);
      // this is a socket to i
      socks[i] = accept(listensock, (sockaddr*)&their_addr, &namelen);
      set_socket_options(socks[i]);
      logstream(LOG_INFO) << "Accepted connection " << i << " -> " << procid() << "\n";
    }
    mpi_barrier();
  }
  // allocate the locks
  sendlocks.resize(socks.size());
  recvlocks.resize(socks.size());
  logstream(LOG_INFO) << "All connections constructed\n";
}



void distributed_control::open_listening_udt() {
  // open listening socket
  listensock = socket(AF_INET, SOCK_STREAM, 0);
  sockaddr_in my_addr;
  my_addr.sin_family = AF_INET;
  my_addr.sin_port = htons(localport);
  my_addr.sin_addr.s_addr = INADDR_ANY;
  memset(&(my_addr.sin_zero), '\0', 8);
  logstream(LOG_INFO) << "Proc " << procid() << " Bind on " << localport << "\n";
  if (bind(listensock, (sockaddr*)&my_addr, sizeof(my_addr)) < 0)
  {
    logstream(LOG_FATAL) << "bind: " << strerror(errno) << "\n";
    ASSERT_TRUE(0);
  }
  logstream(LOG_INFO) << "Proc " << procid() << " listening on " << localport << "\n";
  if (numprocs() > 1) {
    ASSERT_EQ(0, listen(listensock, numprocs() - 1));
  }
}

void distributed_control::close_all_connections() {
  mpi_barrier();
  logger(LOG_INFO, "Closing sockets");
  done = 1;
  mpi_barrier();
  // wake all receivers up everyone up;
  for (procid_t j = 0 ;j < numprocs(); ++j) {
    if (j != procid()) {
      char *a = (char*)malloc(16);
      memset(a, 0, 16);
      send_to_sock(j, a, 1);
    }
  }
  procthreads.join();
  mpi_barrier();
  // socket closing order.to avoid TIME_WAITs
  // close initiators first
  for (procid_t j = procid()+1;j < numprocs(); ++j) {
    close(socks[j]);
  }
  mpi_barrier();
  // then close everything else
  for (procid_t j = 0 ;j < procid(); ++j) {
    close(socks[j]);
  }
  close(listensock);
  // wake up all threads
  procthreads.signalall(SIGURG);
  logger(LOG_INFO, "All sockets closed");
}

void distributed_control::print_stats(procid_t target) {
  // TODO:
}

void distributed_control::create_receive_buffers(size_t rcvbuflen) {
  // create the receive buffers
  buffer.resize(socks.size());

  logger(LOG_INFO, "Allocating buffer of size %d", rcvbuflen);
  for (size_t i = 0;i < buffer.size(); ++i) {
    buffer[i].buffer = (char*)(malloc(rcvbuflen));
    ASSERT_NE(buffer[i].buffer, NULL);
    buffer[i].buflen = rcvbuflen;
    buffer[i].buftail = 0;
    buffer[i].minimum_buflen = rcvbuflen;
    buffer[i].num_recvcalls = 0;
    buffer[i].weighted_bufferutilization = 0;
  }
}


void distributed_control::comms_barrier() {
  logger(LOG_INFO, "%d comms barrier", procid());
  logstream(LOG_INFO) << msgsent.value << " " << msgprocessed.value << std::endl;
  terminator->reset();
  mpi_barrier();
  while (!terminator->done(msgsent.value, msgprocessed.value)) {
    sched_yield();
  }
  terminator->reset();
}

void set_ptr_to_value_1_handler(distributed_control& dc, 
                                procid_t source,  
                                void* ptr,    //serialized any
                                size_t len,   
                                handlerarg_t sizet_ptr) {
    *(reinterpret_cast<volatile size_t*>(sizet_ptr)) = 1;
}


}

