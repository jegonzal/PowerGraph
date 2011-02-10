#ifndef GRAPHLAB_DC_HPP
#define GRAPHLAB_DC_HPP
#include <iostream>
#include <boost/iostreams/stream.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/resizing_array_sink.hpp>
#include <graphlab/util/blocking_queue.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>

#include <graphlab/rpc/dc_receive.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>
#include <graphlab/rpc/dc_dist_object_base.hpp>

#include <graphlab/rpc/is_rpc_call.hpp>
#include <graphlab/rpc/portable_dispatch.hpp>
#include <graphlab/rpc/portable_issue.hpp>
#include <graphlab/rpc/function_call_issue.hpp>
#include <graphlab/rpc/request_issue.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/function_ret_type.hpp>

#include <boost/preprocessor.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>

namespace graphlab {


/**
Struct form of the distributed control constructor parameters.
This will allow curious initialization functions to be written.
For instance:
distributed_control dc(dc_MPI_init(argc, argv));
*/
struct dc_init_param{
  std::vector<std::string> machines;
  std::string initstring;
  procid_t curmachineid;
  size_t numhandlerthreads;
  dc_comm_type commtype;
};

#define DEFAULT_NUMHANDLERTHREADS 8
#define DEFAULT_COMMTYPE TCP_COMM

// forward declaration for dc services
class dc_services;

/**
The primary distributed RPC object.
The basic operation goes like this:
the distributed_control will initialize the communication type using the template
argument provided. It will then spawn a receive multiplxer and a send multiplexer.
The receive multiplexer is attached the communication object using the call back
interface, while the send multiplexer is provided directly access to the comm object.

*/
class distributed_control{
  public:
        /**  Each element of the function call queue is a data/len pair */
    struct function_call_block{
      function_call_block() {}
      function_call_block(procid_t source, unsigned char packet_type_mask, 
                          char* data, size_t len): 
                          source(source), packet_type_mask(packet_type_mask), 
                          data(data), len(len) {}
      procid_t source;
      unsigned char packet_type_mask;
      char* data;
      size_t len;
    };
  private:
   /// initialize receiver threads. private form of the constructor
   void init(const std::vector<std::string> &machines,
             const std::string &initstring,
             procid_t curmachineid,
             size_t numhandlerthreads,
             dc_comm_type commtype = DEFAULT_COMMTYPE);
   
  /// a pointer to the communications subsystem
  dc_impl::dc_comm_base* comm; 
 
  /// senders and receivers to all machines
  std::vector<dc_impl::dc_receive*> receivers;
  std::vector<dc_impl::dc_send*> senders;
  
  /// A thread group of function call handlers
  thread_group fcallhandlers;
  
  /// a queue of functions to be executed
  blocking_queue<function_call_block> fcallqueue;
  
  /// A map of function name to dispatch function. Used for "portable" calls
  dc_impl::dispatch_map_type portable_dispatch_call_map;
  dc_impl::dispatch_map_type portable_dispatch_request_map;

  
  /// object registrations;
  std::vector<void*> registered_objects;
  std::vector<dc_impl::dc_dist_object_base*> registered_rmi_instance;

  /// For convenience, we provide a instance of dc_services 
  dc_services* distributed_services;

  /// ID of the local machine
  procid_t localprocid;
  /// Number of machines
  procid_t localnumprocs;
  
  atomic<size_t> global_calls_sent;
  atomic<size_t> global_calls_received;
  
  
  /// the callback given to the comms class. Called when data is inbound
  friend void dc_recv_callback(void* tag, procid_t src, const char* buf, size_t len);
  
  
  /// the callback given to the comms class. Called when data is inbound
  template <typename T> friend class dc_dist_object;
  
  
  /// disable the operator= by placing it in private 
  distributed_control& operator=(const distributed_control& dc) { return *this; }


  std::map<std::string, std::string> parse_options(std::string initstring);
  
  volatile inline size_t num_registered_objects() {
    return registered_objects.size();
  }
  
  
  // this stores the temporary results for the blocking send_to and recv_from operations
  
  
 public:
   

  
  distributed_control(dc_init_param initparam) {
    init(initparam.machines, 
         initparam.initstring, 
         initparam.curmachineid, 
         initparam.numhandlerthreads,
         initparam.commtype);
  }

  distributed_control(const std::vector<std::string> &machines,
                      const std::string &initstring, 
                      procid_t curmachineid, 
                      size_t numhandlerthreads = DEFAULT_NUMHANDLERTHREADS,
                      dc_comm_type commtype = DEFAULT_COMMTYPE) {
    init(machines, initstring, curmachineid, numhandlerthreads, commtype);
  }

  ~distributed_control();

  /// returns the id of the current processor
  inline procid_t procid() const {
    return localprocid;
  }
  
  /// returns the number of processors in total.
  inline procid_t numprocs() const {
    return localnumprocs;
  }
  
  /**
  This generates the interface functions for the standard calls, basic calls, and fast calls
  */
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
  #define GENI(Z,N,_) BOOST_PP_CAT(i, N)
  #define GENT(Z,N,_) BOOST_PP_CAT(T, N)
  #define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);

  #define RPC_INTERFACE_GENERATOR(Z,N,FNAME_AND_CALL) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  void  BOOST_PP_TUPLE_ELEM(3,0,FNAME_AND_CALL) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, senders.size()); \
    BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,FNAME_AND_CALL),N) \
        <F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(senders[target],  BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL), target, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \
  
  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (remote_call, dc_impl::remote_call_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (fast_remote_call,dc_impl::remote_call_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (control_call, dc_impl::remote_call_issue, FAST_CALL | CONTROL_PACKET) )
 

  #define REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
    BOOST_PP_TUPLE_ELEM(3,0,ARGS) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, senders.size()); \
    return BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,ARGS),N) \
        <F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(senders[target],  BOOST_PP_TUPLE_ELEM(3,2,ARGS), target, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \

  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
 // BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type remote_request, dc_impl::remote_request_issue, 0) )
   BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type remote_request, dc_impl::remote_request_issue, STANDARD_CALL) )

  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type fast_remote_request, dc_impl::remote_request_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type control_request, dc_impl::remote_request_issue, FAST_CALL | CONTROL_PACKET) )
 

  
  #undef RPC_INTERFACE_GENERATOR
  #undef REQUEST_INTERFACE_GENERATOR
  #undef GENARC
  #undef GENT
  #undef GENI
  #undef GENARGS
  
  /**
  Immediately calls the function described by the data
  inside the buffer. This should not be called directly.
  */
  void exec_function_call(procid_t source, unsigned char packet_type_mask, std::istream &istrm);
  
  
  
  /**
  Performs a deferred function call using the information
  inside the buffer. This function will take over ownership of 
  the buffer and will free it when done
  */
  void deferred_function_call(procid_t source, unsigned char packet_type_mask, 
                              char* buf, size_t len);
  

  /**
  This is called by the function handler threads
  */
  void fcallhandler_loop();
  
  inline void inc_calls_sent() {
    global_calls_sent.inc();
  }

  inline void inc_calls_received() {
    global_calls_received.inc();
  }

  inline size_t calls_sent() const {
    return global_calls_sent.value;
  }

  size_t calls_received() const {
    return global_calls_received.value;
  }

  size_t bytes_sent() const {
    size_t ret = 0;
    for (size_t i = 0;i < senders.size(); ++i) ret += senders[i]->bytes_sent();
    return ret;
  }  
  
  size_t bytes_received() const {
    size_t ret = 0;
    for (size_t i = 0;i < receivers.size(); ++i) ret += receivers[i]->bytes_received();
    return ret;
  }  

  /**
    Instantiates a find_dispatch with the right arguments,
    and store the dispatch function in the hash map.
  */
  template <typename F, F f>
  void register_rpc(std::string c) {
    portable_dispatch_request_map[c] = (dc_impl::dispatch_type)
              dc_impl::portable_detail::find_dispatcher<F,        // function type
                              FRESULT,                            // result
                              boost::function_traits<               
                                    typename boost::remove_pointer<F>::type
                                                    >::arity ,   // number of arguments
                              f,                                    // function itself
                              typename dc_impl::is_rpc_call<F>::type  // whether it is an RPC style call
                              >::dispatch_request_fn();
                              
    portable_dispatch_call_map[c] = (dc_impl::dispatch_type)
              dc_impl::portable_detail::find_dispatcher<F,        // function type
                              FRESULT,                            // result
                              boost::function_traits<               
                                    typename boost::remove_pointer<F>::type
                                                    >::arity ,   // number of arguments
                              f,                                    // function itself
                              typename dc_impl::is_rpc_call<F>::type  // whether it is an RPC style call
                              >::dispatch_call_fn();
  }


  inline size_t register_object(void* v, dc_impl::dc_dist_object_base *rmiinstance) {
    ASSERT_NE(v, (void*)NULL);
    registered_objects.push_back(v);
    registered_rmi_instance.push_back(rmiinstance);
    return registered_objects.size() - 1;
  }


  inline void* get_registered_object(size_t id) {
    while(id >= num_registered_objects()) sched_yield();
    ASSERT_NE(registered_objects[id], (void*)NULL);
    return registered_objects[id];
  }

  inline dc_impl::dc_dist_object_base* get_rmi_instance(size_t id) {
    while(id >= num_registered_objects()) sched_yield();
    ASSERT_NE(registered_rmi_instance[id], (void*)NULL);
    return registered_rmi_instance[id];
  }  
  inline void clear_registered_object(size_t id) {
    registered_objects[id] = (void*)NULL;
    registered_rmi_instance[id] = NULL;
  }
  
  
  dc_services& services();
  
  /**
   This comm barrier is not a true "barrier" but is
   essentially a sequentialization point. It guarantees that
   all calls from this machine to the target machine performed
   before the comm_barrier() call are completed before any call
   sent after the comm barrier() call.
  */
  void comm_barrier(procid_t targetmachine);
  
  /**
    This is a convenience function which broadcasts a comm_barrier()
    \note having all machines call the comm barrier does not guarantee
    that all calls have been processed. Basically 'p' local barriers
    do not result in a global barrier.
  */
  void comm_barrier();

  /**
  This is a blocking send_to. It send an object T to the target 
  machine, but waits for the target machine to call recv_from
  before returning. Functionally similar to MPI's matched sending/receiving
  */
  template <typename T>
  void send_to(procid_t target, T& t, bool control = false);
  
  /**
  A blocking recv_from. Must be matched with a send_to call from the
  target before both source and target resumes.
  */
  template <typename T>
  void recv_from(procid_t source, T& t, bool control = false);
  
  /**
  When a process leaves this barrier, it is guaranteed that
   - all processes have called full_barrier()
   - all remote calls/requested have been evaluated and completed
  */
  void full_barrier();
  
};




}

#define REGISTER_RPC(dc, f) dc.register_rpc<typeof(f)*, f>(std::string(BOOST_PP_STRINGIZE(f))) 

#include <graphlab/rpc/function_arg_types_undef.hpp>
#include <graphlab/rpc/function_call_dispatch.hpp>
#include <graphlab/rpc/request_dispatch.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/dc_services.hpp>

/*
Implementations of a couple of template functions in DC
*/
namespace graphlab {
template <typename T>
void distributed_control::send_to(procid_t target, T& t, bool control) {
  services().rmi_instance().send_to(target, t, control);
}

template <typename T>
void distributed_control::recv_from(procid_t source, T& t, bool control) {
  services().rmi_instance().recv_from(source, t, control);
}

} // namespace graphlab
#endif
