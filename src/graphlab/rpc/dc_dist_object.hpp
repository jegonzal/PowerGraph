#ifndef DC_DIST_OBJECT_HPP
#define DC_DIST_OBJECT_HPP
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object_base.hpp>
#include <graphlab/rpc/object_request_issue.hpp>
#include <graphlab/rpc/object_call_issue.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>
/**
Provides capabilities for distributed objects
Your class should either inherit this, or instantiate it
before any distributed object call.
The requirement for using the distributed object is that 
all machines must construct the distributed objects in the same 
order. And, no distributed object calls should be make until
it is guaranteed that all machines have constructed their respective
distributed objects.
*/
namespace graphlab {



template <typename T>
class dc_dist_object : public dc_impl::dc_dist_object_base{
 private:
  distributed_control &dc_;
  size_t obj_id;
  T* owner;
  atomic<size_t> callsreceived;
  atomic<size_t> callssent;
  // make operator= private
  dc_dist_object<T>& operator=(const dc_dist_object<T> &d) {return *this;}
  friend class distributed_control;
 public:
  //internal stuff which should not be used.
  void inc_calls_received() {
    callsreceived.inc();
  }
  void inc_calls_sent() {
    callssent.inc();
  }

 public:
  dc_dist_object(distributed_control &dc_, T* owner):dc_(dc_),owner(owner) {
    obj_id = dc_.register_object(owner, this);
  }
  
  size_t calls_received() const {
    return callsreceived.value;
  }

  size_t calls_sent() const {
    return callssent.value;
  }
  
    
  distributed_control& dc() {
    return dc_;
  }

  const distributed_control& dc() const {
    return dc_;
  }
  
  inline procid_t procid() {
    return dc_.procid();
  }

  inline procid_t numprocs() {
    return dc_.numprocs();
  }

  
  /**
   This comm barrier is not a true "barrier" but is
   essentially a sequentialization point. It guarantees that
   all calls from this machine to the target machine performed
   before the comm_barrier() call are completed before any call
   sent after the comm barrier() call.
  */
  inline void comm_barrier(procid_t targetmachine) {
    return dc_.comm_barrier(targetmachine);
  }
  /**
    This is a convenience function which broadcasts a comm_barrier()
    \note having all machines call the comm barrier does not guarantee
    that all calls have been processed. Basically 'p' local barriers
    do not result in a global barrier.
  */
  inline void comm_barrier() {
    return dc_.comm_barrier();
  }

  inline dc_services& services() {
    return dc_.services();
  }

    /**
  This generates the interface functions for the standard calls, basic calls, and fast calls
  */
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(T, N) BOOST_PP_CAT(i, N)
  #define GENI(Z,N,_) BOOST_PP_CAT(i, N)
  #define GENT(Z,N,_) BOOST_PP_CAT(T, N)
  #define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);

  #define RPC_INTERFACE_GENERATOR(Z,N,FNAME_AND_CALL) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  void  BOOST_PP_TUPLE_ELEM(3,0,FNAME_AND_CALL) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL) & CONTROL_PACKET) == 0) inc_calls_sent(); \
    BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,FNAME_AND_CALL),N) \
        <T, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL), target,obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \
  
  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (remote_call, dc_impl::object_call_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (fast_remote_call,dc_impl::object_call_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (control_call,dc_impl::object_call_issue, FAST_CALL | CONTROL_PACKET) )
 

  #define REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
    BOOST_PP_TUPLE_ELEM(3,0,ARGS) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,ARGS) & CONTROL_PACKET) == 0) inc_calls_sent(); \
    return BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,ARGS),N) \
        <T, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,ARGS), target,obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \

  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type remote_request, dc_impl::object_request_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type fast_remote_request, dc_impl::object_request_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type control_request, dc_impl::object_request_issue, FAST_CALL | CONTROL_PACKET) )
 


  #undef RPC_INTERFACE_GENERATOR
  #undef REQUEST_INTERFACE_GENERATOR
  #undef GENARC
  #undef GENT
  #undef GENI
  #undef GENARGS
  

};


#include <graphlab/rpc/mem_function_arg_types_undef.hpp>

}// namespace graphlab
#endif
