#ifndef DC_HPP
#define DC_HPP
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
      function_call_block(procid_t source, 
                          char* data, size_t len): 
                          source(source), data(data), len(len) {}
      procid_t source;
      char* data;
      size_t len;
    };
  private:
   /// initialize receiver threads. private form of the constructor
   void init(const std::vector<std::string> &machines,
             const std::string &initstring,
             procid_t curmachineid,
             size_t numhandlerthreads,
             dc_comm_type commtype = TCP_COMM);
   
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

  /// For convenience, we provide a instance of dc_services 
  dc_services* distributed_services;

  /// ID of the local machine
  procid_t localprocid;
  /// Number of machines
  procid_t localnumprocs;
  
  /// the callback given to the comms class. Called when data is inbound
  friend void dc_recv_callback(void* tag, procid_t src, const char* buf, size_t len);
  
  
  /// the callback given to the comms class. Called when data is inbound
  template <typename T> friend class dc_dist_object;
  
  
  /// disable the operator= by placing it in private 
  distributed_control& operator=(const distributed_control& dc) { return *this; }


  std::map<std::string, std::string> parse_options(std::string initstring);
  
  
  
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
                      size_t numhandlerthreads = 8,
                      dc_comm_type commtype = TCP_COMM) {
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
  void exec_function_call(procid_t source, std::istream &istrm);
  
  
  
  /**
  Performs a deferred function call using the information
  inside the buffer. This function will take over ownership of 
  the buffer and will free it when done
  */
  void deferred_function_call(procid_t source, char* buf, size_t len);
  

  /**
  This is called by the function handler threads
  */
  void fcallhandler_loop();
  
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


  inline size_t register_object(void* v) {
    ASSERT_NE(v, (void*)NULL);
    registered_objects.push_back(v);
    return registered_objects.size() - 1;
  }


  inline void* get_registered_object(size_t id) {
    ASSERT_LT(id, registered_objects.size());
    ASSERT_NE(registered_objects[id], (void*)NULL);
    return registered_objects[id];
  }
  
  inline void clear_registered_object(size_t id) {
    registered_objects[id] = (void*)NULL;
  }
  
  
  dc_services& services();
  
  void comm_barrier(procid_t targetmachine);
  void comm_barrier();
};




}

#define REGISTER_RPC(dc, f) dc.register_rpc<typeof(f)*, f>(std::string(BOOST_PP_STRINGIZE(f))) 

#include <graphlab/rpc/function_arg_types_undef.hpp>
#include <graphlab/rpc/function_call_dispatch.hpp>
#include <graphlab/rpc/request_dispatch.hpp>
#include <graphlab/rpc/dc_services.hpp>
#endif
