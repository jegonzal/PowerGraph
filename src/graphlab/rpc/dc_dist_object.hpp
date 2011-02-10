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


#define DC_SERVICES_BARRIER_BRANCH_FACTOR 128

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
  size_t control_obj_id;  // object id of this object
  T* owner;
  atomic<size_t> callsreceived;
  atomic<size_t> callssent;
  using dc_impl::dc_dist_object_base::recv_froms;
  // make operator= private
  dc_dist_object<T>& operator=(const dc_dist_object<T> &d) {return *this;}
  friend class distributed_control;
  
  
  // ------- Sense reversing barrier data ----------
  /// The next value of the barrier. either +1 or -1
  int barrier_sense;
  /// When this flag == the current barrier value. The barrier is complete
  int barrier_release;
  /** when barrier sense is 1, barrier clears when 
   * child_barrier_counter == numchild. When barrier sense is -1, barrier
   * clears when child_barrier_counter == 0;
   */
  atomic<int> child_barrier_counter;
  /// condition variable and mutex protecting the barrier variables
  conditional barrier_cond;
  mutex barrier_mut;
  size_t parent;  /// parent node
  size_t childbase; /// id of my first child
  size_t numchild;  /// number of children




  /**
    The child calls this function in the parent once the child enters the barrier
  */
  void __child_to_parent_barrier_trigger(procid_t source) {
    barrier_mut.lock();
    // assert childbase <= source <= childbase + DC_SERVICES_BARRIER_BRANCH_FACTOR
    ASSERT_GE(source, childbase);
    ASSERT_LT(source, childbase + DC_SERVICES_BARRIER_BRANCH_FACTOR);
    child_barrier_counter.inc(barrier_sense);
    barrier_cond.signal();
    barrier_mut.unlock();
  } 
  
  /**
    This is on the downward pass of the barrier. The parent calls this function
    to release all the children's barriers
  */
  void __parent_to_child_barrier_release(int releaseval) {
    // send the release downwards
    // get my largest child
    for (size_t i = 0;i < numchild; ++i) {
      internal_control_call(childbase + i,
                            &dc_dist_object<T>::__parent_to_child_barrier_release,
                            releaseval);
  
    }
    barrier_mut.lock();
    barrier_release = releaseval;    
    barrier_cond.signal();
    barrier_mut.unlock();
  }
  
  
  
  /**
    Used to implement the match send/recv 
    Sender call this function in the receiver which stores
    and sets the condition variable
  */
  void block_and_wait_for_recv(size_t src,
                             size_t objid,
                             std::string& str,
                             size_t tag) {
    recv_froms[src].lock.lock();
    recv_froms[src].data = str;
    recv_froms[src].tag = tag;
    recv_froms[src].hasdata = true;
    recv_froms[src].cond.signal();
    recv_froms[src].lock.unlock();
  }

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
    control_obj_id = dc_.register_object(this, this);
    
    recv_froms.resize(dc_.numprocs());
    
    // reset the child barrier values
    child_barrier_counter.value = 0;
    barrier_sense = 1;
    barrier_release = -1;

  
    // compute my children
    childbase = size_t(dc_.procid()) * DC_SERVICES_BARRIER_BRANCH_FACTOR + 1;
    if (childbase >= dc_.numprocs()) {
      numchild = 0;
    }
    else {
      size_t maxchild = std::min<size_t>(dc_.numprocs(), 
                                         childbase + DC_SERVICES_BARRIER_BRANCH_FACTOR);
      numchild = maxchild - childbase;
    }
  
    parent =  (dc_.procid() - 1) / DC_SERVICES_BARRIER_BRANCH_FACTOR   ;

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

  /**
    This returns the set of services for the parent DC.
  */
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
  
  // Now generate the interface functions which allow me to call this dc_dist_object directly

  #define RPC_INTERFACE_GENERATOR(Z,N,FNAME_AND_CALL) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  void  BOOST_PP_TUPLE_ELEM(3,0,FNAME_AND_CALL) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL) & CONTROL_PACKET) == 0) inc_calls_sent(); \
    BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,FNAME_AND_CALL),N) \
        <dc_dist_object<T>, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL), target,control_obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \
  
  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (internal_call,dc_impl::object_call_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (internal_fast_call,dc_impl::object_call_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (internal_control_call,dc_impl::object_call_issue, FAST_CALL | CONTROL_PACKET) )
 

  #define REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
    BOOST_PP_TUPLE_ELEM(3,0,ARGS) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,ARGS) & CONTROL_PACKET) == 0) inc_calls_sent(); \
    return BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,ARGS),N) \
        <dc_dist_object<T>, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,ARGS), target,control_obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \

  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type internal_request, dc_impl::object_request_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type internal_fast_request, dc_impl::object_request_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type internal_control_request, dc_impl::object_request_issue, FAST_CALL | CONTROL_PACKET) )
 

  #undef RPC_INTERFACE_GENERATOR
  #undef REQUEST_INTERFACE_GENERATOR
  #undef GENARC
  #undef GENT
  #undef GENI
  #undef GENARGS
  



  /**
  This is a blocking send_to. It send an object T to the target 
  machine, but waits for the target machine to call recv_from
  before returning. Functionally similar to MPI's matched sending/receiving
  */
  template <typename U>
  void send_to(procid_t target, U& t, bool control = false) {
    std::stringstream strm;
    oarchive oarc(strm);
    oarc << t;
    dc_impl::reply_ret_type rt(REQUEST_WAIT_METHOD);
    // I shouldn't use a request to block here since 
    // that will take up a thread on the remote side
    // so I simulate a request here.
    size_t rtptr = reinterpret_cast<size_t>(&rt);
    if (control == false) {
      internal_call(target, &dc_dist_object<T>::block_and_wait_for_recv, 
                     procid(),obj_id, strm.str(), rtptr);
    }
    else {
      internal_control_call(target, &dc_dist_object<T>::block_and_wait_for_recv, 
                  procid(), obj_id, strm.str(), rtptr);
    }
    // wait for reply
    rt.wait();
    
    if (control == false) inc_calls_received();
  }
  
  
  /**
  A blocking recv_from. Must be matched with a send_to call from the
  target before both source and target resumes.
  */
  template <typename U>
  void recv_from(procid_t source, U& t, bool control = false) {
    // wait on the condition variable until I have data
    dc_impl::recv_from_struct &recvstruct = recv_froms[source];
    recvstruct.lock.lock();
    while (recvstruct.hasdata == false) {
      recvstruct.cond.wait(recvstruct.lock);
    }
    
    // got the data. deserialize it
    std::stringstream strm(recvstruct.data);
    iarchive iarc(strm);
    iarc >> t;
    // clear the data
    std::string("").swap(recvstruct.data);
    // remember the tag so we can unlock it before the remote call
    size_t tag = recvstruct.tag;
    // clear the has data flag
    recvstruct.hasdata = false;
    // unlock
    recvstruct.lock.unlock();
    if (control == false) {
      // remote call to release the sender. Use an empty blob
      dc_.fast_remote_call(source, reply_increment_counter, tag, dc_impl::blob());
      // I have to increment the calls sent manually here
      // since the matched send/recv calls do not go through the 
      // typical object calls. It goes through the DC, but I also want to charge
      // it to this object
      inc_calls_sent();
    }
    else {
      dc_.control_call(source, reply_increment_counter, tag, dc_impl::blob());
    }
  }


  
  /**
  A local barrier on this object
  */
  void barrier() {
    // upward message
    char barrier_val = barrier_sense;      
    barrier_mut.lock();
    // wait for all children to be done
    while(1) {
      if ((barrier_sense == -1 && child_barrier_counter.value == 0) || 
          (barrier_sense == 1 && child_barrier_counter.value == (int)(numchild))) {
        // flip the barrier sense
        barrier_sense = -barrier_sense;
        // call child to parent in parent
        barrier_mut.unlock();
        if (procid() != 0) {
          internal_control_call(parent, 
                           &dc_dist_object<T>::__child_to_parent_barrier_trigger,
                           procid());
        }
        break;
      }
      barrier_cond.wait(barrier_mut);
    }
    
    
    logger(LOG_INFO, "barrier phase 1 complete");
    // I am root. send the barrier release downwards
    if (procid() == 0) {
      barrier_release = barrier_val;
  
      for (size_t i = 0;i < numchild; ++i) {
        internal_control_call(childbase + i,
                             &dc_dist_object<T>::__parent_to_child_barrier_release,
                             barrier_val);
      
      }
    }
    // wait for the downward message releasing the barrier
    barrier_mut.lock();
    while(1) {
      if (barrier_release == barrier_val) break;
      barrier_cond.wait(barrier_mut);
    }
    barrier_mut.unlock();
  
    logger(LOG_INFO, "barrier phase 2 complete");
  }

};


#include <graphlab/rpc/mem_function_arg_types_undef.hpp>

}// namespace graphlab
#endif
