#ifndef GRAPHLAB_DISTRIBUTED_CONTROL_HPP
#define GRAPHLAB_DISTRIBUTED_CONTROL_HPP
#include <vector>
#include <mpi.h>
#include <boost/typeof/typeof.hpp>
#include <logger/logger.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <boost/static_assert.hpp>

#include <graphlab/distributed/repack_dispatch.hpp>
#include <graphlab/distributed/serialize_dispatch.hpp>
#include <graphlab/distributed/distributed_terminator.hpp>
#include <graphlab/distributed/distributed_control_types.hpp>

#include <serialization/oarchive.hpp>
#include <serialization/iarchive.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/blocking_queue.hpp>


namespace graphlab {

#define DC_LOCAL_BASE_PORT_NUM 13000

  /*
  Handler functions take the following structure
  There can be up to 8 integer arguments
  */
typedef void (*handler_type)(distributed_control& dc, size_t source,
                            void* ptr, size_t len,  ... );

typedef void (*xdispatcher)(distributed_control& dc, size_t source,
                            void* ptr, size_t len, void* stack);

typedef void (*xsdispatcher)(distributed_control& dc, size_t source,
                            void* ptr, size_t len, iarchive &arc);



class distributed_control {
 public:
  typedef unsigned char handler_t;
  

  /**
    Begins network initialization. pargc, pargv may be NULL
  */
  distributed_control(int *pargc, char*** pargv);
  
  ~distributed_control();

  /// Begins background message processing threads
  void init_message_processing(size_t numthreads = 1);
  
  inline const procid_t procid() const {
    return id;
  }
  
  inline const procid_t numprocs() const {
    return nprocs;
  }


  inline void mpi_barrier() {
  logger(LOG_INFO, "%d barrier", procid());
    MPI_Barrier(MPI_COMM_WORLD);
  }

  void comms_barrier();

  inline void barrier() {
     comms_barrier();
  }
  
  void report_stats();
  
  // ----------    Standard Remote Call with integer parameters -------------

  template<typename F>
  void remote_call(procid_t target, F remote_function, void* ptr, size_t len) {
    // this also serves as a typecheck for F
    if (target == id) {remote_function(*this, id, ptr, len); return;}
    handlerarg_t arr[8];
    send_call_message(target, (void*)remote_function, ptr, len, 0, arr);
  }
  
  template<typename F>
  void remote_call(procid_t target, F remote_function, void* ptr, size_t len, 
                     handlerarg_t i0) {
    if (target == id) {remote_function(*this, id, ptr, len, i0);return;}
    handlerarg_t arr[8]; 
    arr[0] = i0;
    send_call_message(target, (void*)remote_function, ptr, len, 1, arr);    
  }

  template<typename F>
  void remote_call(procid_t target, F remote_function, void* ptr, size_t len, 
                     handlerarg_t i0, handlerarg_t i1) {
    if (target == id) {remote_function(*this,id, ptr, len, i0, i1);  return;}
    handlerarg_t arr[8]; 
    arr[0] = i0; arr[1] = i1;
    send_call_message(target, (void*)remote_function, ptr, len, 2, arr);
  }
  
  template<typename F>
  void remote_call(procid_t target, F remote_function, void* ptr, size_t len, 
                     handlerarg_t i0, handlerarg_t i1, handlerarg_t i2) {
    if (target == id){ remote_function(*this, id, ptr, len, i0, i1, i2);
                        return;}
    handlerarg_t arr[8]; 
    arr[0] = i0; arr[1] = i1; arr[2] = i2;
    send_call_message(target, (void*)remote_function, ptr, len, 3, arr);
  }
  
  template<typename F>
  void remote_call(procid_t target, F remote_function, void* ptr, size_t len, 
                     handlerarg_t i0, handlerarg_t i1, handlerarg_t i2, 
                     handlerarg_t i3) {
    if (target == id){ remote_function(*this, id, ptr, len, i0, i1, i2, i3);
                       return;}
    handlerarg_t arr[8]; 
    arr[0] = i0; arr[1] = i1; arr[2] = i2;
    arr[3] = i3;
    send_call_message(target, (void*)remote_function, ptr, len, 4, arr);
  }
  
  template<typename F>
  void remote_call(procid_t target, F remote_function, void* ptr, size_t len, 
                     handlerarg_t i0, handlerarg_t i1, handlerarg_t i2, 
                     handlerarg_t i3, handlerarg_t i4) {
    if (target == id){ remote_function(*this, id, ptr, len, i0, i1, i2, i3, 
                                                           i4);  return;}
    handlerarg_t arr[8]; 
    arr[0] = i0; arr[1] = i1; arr[2] = i2;
    arr[3] = i3; arr[4] = i4;
    send_call_message(target, (void*)remote_function, ptr, len, 5, arr);
  }

  template<typename F>
  void remote_call(procid_t target, F remote_function, void* ptr, size_t len, 
                     handlerarg_t i0, handlerarg_t i1, handlerarg_t i2, 
                     handlerarg_t i3, handlerarg_t i4, handlerarg_t i5) {
    if (target == id){ remote_function(*this, id, ptr, len, i0, i1, i2, i3, 
                                                           i4, i5);  return;}
    handlerarg_t arr[8]; 
    arr[0] = i0; arr[1] = i1; arr[2] = i2;
    arr[3] = i3; arr[4] = i4; arr[5] = i5;
    send_call_message(target, (void*)remote_function, ptr, len, 6, arr);
  }
  
  // ----------  Extended Remote Call with stack parameters -------------
  
  template<typename F>
  void remote_callx(procid_t target, F remote_function, void* ptr, size_t len) {
    //BOOST_STATIC_ASSERT((boost::function_types::is_function<F,  stdcall_cc>::value == true));
    if (target == id) { remote_function(*this, id, ptr, len);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function);
    xdispatcher dispptr = (DISPATCH0<struct_type>);
    send_callx_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0>
  void remote_callx(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0) {
    if (target == id) { remote_function(*this, id, ptr, len, i0);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0);
    xdispatcher dispptr = (DISPATCH1<struct_type>);
    send_callx_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename  T1>
  void remote_callx(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1);
    xdispatcher dispptr = (DISPATCH2<struct_type>);
    send_callx_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2>
  void remote_callx(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2);
    xdispatcher dispptr = (DISPATCH3<struct_type>);
    send_callx_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3>
  void remote_callx(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2, T3 i3) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2, i3);
                       return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2, i3)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2, i3);
    xdispatcher dispptr = (DISPATCH4<struct_type>);
    send_callx_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3, 
                      typename T4>
  void remote_callx(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2, T3 i3, T4 i4) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2, i3, i4);
                       return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4);
    xdispatcher dispptr = (DISPATCH5<struct_type>);
    send_callx_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3, 
                      typename T4, typename T5>
  void remote_callx(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2, T3 i3, T4 i4, T5 i5) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2, i3, i4, i5);
                      return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4, i5)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4, i5);
    xdispatcher dispptr = (DISPATCH6<struct_type>);
    send_callx_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }


  // ----------  Extended Remote Call with serializable parameters -------------

  
  template<typename F>
  void remote_callxs(procid_t target, F remote_function, void* ptr, size_t len) {
    //BOOST_STATIC_ASSERT((boost::function_types::is_function<F,  stdcall_cc>::value == true));
    std::stringstream str;
    oarchive arc(str);
    if (target == id) { remote_function(*this, id, ptr, len);  return;}
    typedef typename SER_GET_STRUCT_TYPE0<F>::struct_type struct_type;
    SERREPACKSTRUCT(remote_function, arc);
    xsdispatcher dispptr = (SERDISPATCH0<struct_type>);
    send_callxs_message(target, (void*)(dispptr), ptr, len, str.str());
  }

  template<typename F, typename T0>
  void remote_callxs(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0) {
    std::stringstream str;
    oarchive arc(str);
    if (target == id) { remote_function(*this, id, ptr, len, i0);  return;}
    typedef typename SER_GET_STRUCT_TYPE1<F>::struct_type struct_type;
    SERREPACKSTRUCT(remote_function, arc, i0);
    xsdispatcher dispptr = (SERDISPATCH1<struct_type>);
    send_callxs_message(target, (void*)(dispptr), ptr, len, str.str());
  }

  template<typename F, typename T0, typename  T1>
  void remote_callxs(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1) {
    std::stringstream str;
    oarchive arc(str);
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1);  return;}
    typedef typename SER_GET_STRUCT_TYPE2<F>::struct_type struct_type;
    SERREPACKSTRUCT(remote_function, arc, i0, i1);
    xsdispatcher dispptr = (SERDISPATCH2<struct_type>);
    send_callxs_message(target, (void*)(dispptr), ptr, len, str.str());
  }

  template<typename F, typename T0, typename T1, typename T2>
  void remote_callxs(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2) {
    std::stringstream str;
    oarchive arc(str);
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2);  return;}
    typedef typename SER_GET_STRUCT_TYPE3<F>::struct_type struct_type;
    SERREPACKSTRUCT(remote_function, arc, i0, i1, i2);
    xsdispatcher dispptr = (SERDISPATCH3<struct_type>);
    send_callxs_message(target, (void*)(dispptr), ptr, len, str.str());
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3>
  void remote_callxs(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2, T3 i3) {
    std::stringstream str;
    oarchive arc(str);
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2, i3);  return;}
    typedef typename SER_GET_STRUCT_TYPE4<F>::struct_type struct_type;
    SERREPACKSTRUCT(remote_function, arc, i0, i1, i2, i3);
    xsdispatcher dispptr = (SERDISPATCH4<struct_type>);
    send_callxs_message(target, (void*)(dispptr), ptr, len, str.str());
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3, 
                      typename T4>
  void remote_callxs(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2, T3 i3, T4 i4) {
    std::stringstream str;
    oarchive arc(str);
    if (target == id) {
      remote_function(*this, id, ptr, len, i0, i1, i2, i3, i4);
      return;
    }
    typedef typename SER_GET_STRUCT_TYPE5<F>::struct_type struct_type;
    SERREPACKSTRUCT(remote_function, arc, i0, i1, i2, i3, i4);
    xsdispatcher dispptr = (SERDISPATCH5<struct_type>);
    send_callxs_message(target, (void*)(dispptr), ptr, len, str.str());
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3, 
                      typename T4, typename T5>
  void remote_callxs(procid_t target, F remote_function, void* ptr, size_t len, 
                      T0 i0, T1 i1, T2 i2, T3 i3, T4 i4, T5 i5) {
    std::stringstream str;
    oarchive arc(str);
    if (target == id) {
      remote_function(*this, id, ptr, len, i0, i1, i2, i3, i4, i5);
      return;
    }
    typedef typename SER_GET_STRUCT_TYPE6<F>::struct_type struct_type;
    SERREPACKSTRUCT(remote_function, arc, i0, i1, i2, i3, i4, i5);
    xsdispatcher dispptr = (SERDISPATCH6<struct_type>);
    send_callxs_message(target, (void*)(dispptr), ptr, len, str.str());
  }


  // ------------------------ Remote Call of Control messages: Only stack parameters allowed
  template<typename F>
  void remote_call_control(procid_t target, F remote_function, void* ptr, size_t len) {
    //BOOST_STATIC_ASSERT((boost::function_types::is_function<F,  stdcall_cc>::value == true));
    if (target == id) { remote_function(*this, id, ptr, len);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function);
    xdispatcher dispptr = (DISPATCH0<struct_type>);
    send_call_control_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0>
  void remote_call_control(procid_t target, F remote_function, void* ptr, size_t len,
                      T0 i0) {
    if (target == id) { remote_function(*this, id, ptr, len, i0);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0);
    xdispatcher dispptr = (DISPATCH1<struct_type>);
    send_call_control_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename  T1>
  void remote_call_control(procid_t target, F remote_function, void* ptr, size_t len,
                      T0 i0, T1 i1) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1);
    xdispatcher dispptr = (DISPATCH2<struct_type>);
    send_call_control_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2>
  void remote_call_control(procid_t target, F remote_function, void* ptr, size_t len,
                      T0 i0, T1 i1, T2 i2) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2);  return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2);
    xdispatcher dispptr = (DISPATCH3<struct_type>);
    send_call_control_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3>
  void remote_call_control(procid_t target, F remote_function, void* ptr, size_t len,
                      T0 i0, T1 i1, T2 i2, T3 i3) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2, i3);
                       return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2, i3)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2, i3);
    xdispatcher dispptr = (DISPATCH4<struct_type>);
    send_call_control_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3,
                      typename T4>
  void remote_call_control(procid_t target, F remote_function, void* ptr, size_t len,
                      T0 i0, T1 i1, T2 i2, T3 i3, T4 i4) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2, i3, i4);
                       return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4);
    xdispatcher dispptr = (DISPATCH5<struct_type>);
    send_call_control_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  template<typename F, typename T0, typename T1, typename T2, typename T3,
                      typename T4, typename T5>
  void remote_call_control(procid_t target, F remote_function, void* ptr, size_t len,
                      T0 i0, T1 i1, T2 i2, T3 i3, T4 i4, T5 i5) {
    if (target == id) { remote_function(*this, id, ptr, len, i0, i1, i2, i3, i4, i5);
                      return;}
    typedef BOOST_TYPEOF(REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4, i5)) struct_type;
    struct_type ps = REPACKSTRUCT(remote_function, i0, i1, i2, i3, i4, i5);
    xdispatcher dispptr = (DISPATCH6<struct_type>);
    send_call_control_message(target, (void*)(dispptr), ptr, len, &ps, sizeof(ps));
  }

  struct recv_buffer{
    char* buffer;
    size_t buflen;
    size_t buftail;
    size_t minimum_buflen;
    
    size_t num_recvcalls;
    size_t weighted_bufferutilization;
  };
  std::vector<recv_buffer> buffer;

  class messageproc_thread :public runnable {
  public:
    size_t* done;
    distributed_control *dc;
    
    // constructor /destructors
    messageproc_thread();
    ~messageproc_thread();

    void run();
  };

  void print_stats(procid_t target);

  
  struct send_req_data{
    send_req_data() {}
    send_req_data(procid_t target, char* buf, size_t len):
                            target(target), buf(buf), len(len) { }
    procid_t target;
    char* buf;
    size_t len;
  };

  
  class background_send_thread:public thread {
    distributed_control &dc;
   public: 
    size_t bytes_sent;
    background_send_thread(distributed_control &dc):dc(dc),bytes_sent(0) { }
    void run() {
      logger(LOG_INFO, "send thread started");
      while(1) {
        std::pair<send_req_data, bool> ret = dc.send_requests.dequeue();
        if (ret.second == false) break;
        bytes_sent+=ret.first.len;
        dc.send_to_sock(ret.first.target, ret.first.buf, ret.first.len);
      }
    }
  };
  
  struct dispatch_req_data{
    char* buf;
    size_t len;
  };

  class message_dispatch_thread:public thread {
    distributed_control &dc;
   public: 
    message_dispatch_thread(distributed_control &dc):dc(dc) { }
    void run();
  };

 private:
  /// number of processes in the network
  procid_t nprocs;
  /// id of the current process (MPI rank)
  procid_t id;

  /// all_addrs[i] will contain the IP address of machine i
  uint32_t *all_addrs;
  /// The local port we are listening on
  size_t localport;

  /// the socket we use to listen on 
  sockfd_t listensock;
  
  /// socks[i] is the socket to machine i.
  /// There is no socket to the local process ( socks[procid()] is invalid )
  std::vector<sockfd_t> socks; 

  size_t done;
  std::vector<messageproc_thread> procs;
  thread_group procthreads;


  std::vector<spinlock> sendlocks;  // TODO: unused now since only 1 send thread
  std::vector<spinlock> recvlocks;  // 
   
  // background sending
  background_send_thread *send_thread;
  blocking_queue<send_req_data> send_requests;
  
  // background receiving
  std::vector<message_dispatch_thread*> dispatch_thread;
  blocking_queue<dispatch_req_data> dispatch_requests;

  atomic<size_t> msgsent;
  atomic<size_t> msgprocessed;
  distributed_terminator* terminator;
    
  void set_socket_options(int fd);
  size_t read_tcp_buflen(int fd);
  /// returns the local ip in ip
  void get_local_ip(char ip[4]);
  /// synchronizes the IP list (all_addrs) among the entire MPI ring
  void sync_ip_list();
  /// connects all the machines in the MPI ring
  void connect_udt();
  /// Opens the local listening UDT port
  void open_listening_udt();
  /// Closes all UDT connections
  void close_all_connections();
  void create_receive_buffers(size_t rcvbuflen);
  
  void send_to_sock(procid_t sock, char* buf, size_t len);
  // receive functions
  void receive_handler(sockfd_t sock, recv_buffer &buffer);
  void receive_call_message(char* msg, size_t len);
  void receive_callx_message(char* msg, size_t len);
  void receive_callxs_message(char* msg, size_t len);
  // send functions
  void send_call_message(procid_t target, void* remote_function, void* ptr, 
                      size_t len,size_t numargs, handlerarg_t *arr);
  void send_callx_message(procid_t target, void* remote_function, void* ptr, 
                      size_t len,void* stackbegin, size_t stacklen);
  void send_callxs_message(procid_t target, void* remote_function, void* ptr, 
                      size_t len,const std::string &s);
  // control messages are the same as callx messages but they do not count
  // towards the message counter
  void send_call_control_message(procid_t target, void* remote_function, void* ptr,
                      size_t len,void* stackbegin, size_t stacklen);

//  void pack_dgram(size_t handlerid, std::vector<uint64_t> &ints,
//                  void* b, int len);
  
};



void set_ptr_to_value_1_handler(distributed_control& dc, 
                                procid_t source,  
                                void* ptr,    //serialized any
                                size_t len,   
                                handlerarg_t sizet_ptr);

extern distributed_control *dc_singleton_ptr;

}

#endif
