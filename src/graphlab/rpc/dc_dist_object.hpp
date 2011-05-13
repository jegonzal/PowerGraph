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

#include <graphlab/rpc/dc.hpp>

#ifndef DC_DIST_OBJECT_HPP
#define DC_DIST_OBJECT_HPP
#include <vector>
#include <string>
#include <set>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_dist_object_base.hpp>
#include <graphlab/rpc/object_request_issue.hpp>
#include <graphlab/rpc/object_call_issue.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#include <graphlab/util/charstream.hpp>
#include <boost/preprocessor.hpp>
#include <graphlab/macros_def.hpp>

#define BARRIER_BRANCH_FACTOR 128


namespace graphlab {


/**
\ingroup rpc
Provides capabilities for distributed objects Your class should either
inherit this, or instantiate it before any distributed object call.
The requirement for using the distributed object is that all machines
must construct the distributed objects in the same order. And, no
distributed object calls should be make until it is guaranteed that
all machines have constructed their respective distributed objects.

This class also acts as a single "context" spanning multiple machines.
For instance, the barrier implemented here is fully localized to within 
a particular instance of this object. That is, multiple instances of the object
can issue multiple independent barriers.
The dc_services() object is a thin wrapper around the dc_dist_object.

This class implements several MPI-like primitive ops such as 
barrier, gather, broadcast, etc. These operations are not particular optimized
and can be quite inefficient.
*/
template <typename T>
class dc_dist_object : public dc_impl::dc_dist_object_base{
 private:
  distributed_control &dc_;
  size_t obj_id;
  size_t control_obj_id;  // object id of this object
  T* owner;
  bool calltracking;
  std::vector<atomic<size_t> > callsreceived;
  std::vector<atomic<size_t> > callssent;
  std::vector<atomic<size_t> > bytessent;
  // make operator= private
  dc_dist_object<T>& operator=(const dc_dist_object<T> &d) {return *this;}
  friend class distributed_control;
  
  
  


 public:
  /// Should not be used by the user
  void inc_calls_received(procid_t p) {
    if (!full_barrier_in_effect) {
        callsreceived[p].inc();
    }
    else {
      //check the proc I just incremented.
      // If I just exceeded the required size, I need
      // to decrement the full barrier counter
      if (callsreceived[p].inc() == calls_to_receive[p]) {
        // if it was me who set the bit
        if (procs_complete.set_bit(p) == false) {
          // then decrement the incomplete count.
          // if it was me to decreased it to 0
          // lock and signal
          full_barrier_lock.lock();
          if (num_proc_recvs_incomplete.dec() == 0) {
            full_barrier_cond.signal();
          }
          full_barrier_lock.unlock();
        }
      }
    }
  }
  
  /// Should not be used by the user
  void inc_calls_sent(procid_t p) {
    callssent[p].inc();
  }
  
  /// Should not be used by the user
  void inc_bytes_sent(procid_t p, size_t bytes) {
    bytessent[p].inc(bytes);
  }

 public:
  dc_dist_object(distributed_control &dc_, T* owner, bool calltracking = false):
                                dc_(dc_),owner(owner),calltracking(calltracking) {
    callssent.resize(dc_.numprocs());
    callsreceived.resize(dc_.numprocs());
    bytessent.resize(dc_.numprocs());
    //------ Initialize the matched send/recv ------
    recv_froms.resize(dc_.numprocs());
    
    //------ Initialize the gatherer ------
    gather_receive.resize(dc_.numprocs());

    
    //------- Initialize the Barrier ----------
    child_barrier_counter.value = 0;
    barrier_sense = 1;
    barrier_release = -1;

  
    // compute my children
    childbase = size_t(dc_.procid()) * BARRIER_BRANCH_FACTOR + 1;
    if (childbase >= dc_.numprocs()) {
      numchild = 0;
    }
    else {
      size_t maxchild = std::min<size_t>(dc_.numprocs(), 
                                         childbase + BARRIER_BRANCH_FACTOR);
      numchild = (procid_t)(maxchild - childbase);
    }
  
    parent =  (procid_t)((dc_.procid() - 1) / BARRIER_BRANCH_FACTOR)   ;

    //-------- Initialize all gather --------------
    ab_child_barrier_counter.value = 0;
    ab_barrier_sense = 1;
    ab_barrier_release = -1;

    
    //-------- Initialize the full barrier ---------
    
    full_barrier_in_effect = false;
    procs_complete.resize(dc_.numprocs());
    
    // register
    obj_id = dc_.register_object(owner, this);
    control_obj_id = dc_.register_object(this, this);
  }
  
  /// The number of function calls received by this object
  size_t calls_received() const {
    size_t ctr = 0;
    for (size_t i = 0;i < numprocs(); ++i) {
      ctr += callsreceived[i].value;
    }
    return ctr;
  }

  /// The number of function calls send from this object
  size_t calls_sent() const {
    size_t ctr = 0;
    for (size_t i = 0;i < numprocs(); ++i) {
      ctr += callssent[i].value;
    }
    return ctr;
  }
  
  size_t bytes_sent() const {
    size_t ctr = 0;
    for (size_t i = 0;i < numprocs(); ++i) {
      ctr += bytessent[i].value;
    }
    return ctr;
  }
  
  /// A reference to the underlying dc
  distributed_control& dc() {
    return dc_;
  }

  /// A reference to the underlying dc
  const distributed_control& dc() const {
    return dc_;
  }
  
  /// The current process ID
  inline procid_t procid() const {
    return dc_.procid();
  }

  /// The number of processes in the distributed program.
  inline procid_t numprocs() const {
    return dc_.numprocs();
  }

  
  /**
   This comm barrier is not a true "barrier" but is
   essentially a sequentialization point. It guarantees that
   all calls from this machine to the target machine performed
   before the comm_barrier() call are completed before any call
   sent after the comm barrier() call.
   
    \note This affects the global context
  */
  inline void comm_barrier(procid_t targetmachine) {
    return dc_.comm_barrier(targetmachine);
  }
  /**
    This is a convenience function which broadcasts a comm_barrier()
    \note having all machines call the comm barrier does not guarantee
    that all calls have been processed. Basically 'p' local barriers
    do not result in a global barrier.
    
    \note This affects the global context
  */
  inline void comm_barrier() {
    return dc_.comm_barrier();
  }

  /**
    This returns the set of services for the parent DC.
    This is deprecated. Use dc() to get access to the global context
  */
   __attribute__((__deprecated__)) inline dc_services& services() {
    return dc_.services();
  }

    /*
  This generates the interface functions for the standard calls, basic calls, and fast calls
  The function looks like this:
  \code
  template<typename F , typename T0> void remote_call (procid_t target, F remote_function , T0 i0 )
  {
      ASSERT_LT(target, dc_.senders.size());
      if ((STANDARD_CALL & CONTROL_PACKET) == 0) inc_calls_sent(target);
      dc_impl::object_call_issue1 <T, F , T0> ::exec(dc_.senders[target], 
                                                      STANDARD_CALL, 
                                                      target,obj_id, 
                                                      remote_function , 
                                                      i0 );
  }

  The argument to the RPC_INTERFACE_GENERATOR are:
    - the name of the rpc call ("remote_call" in the first one)
    - the name of the issueing processor ("object_call_issue")
    - The flags to set on the call ("STANDARD_CALL")
    
    The call can be issued with
    rmi.remote_call(target,
                    &object_type::function_name,
                    arg1,
                    arg2...)
  \endcode
  */
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(T, N) BOOST_PP_CAT(i, N)
  #define GENI(Z,N,_) BOOST_PP_CAT(i, N)
  #define GENT(Z,N,_) BOOST_PP_CAT(T, N)
  #define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);

  #define RPC_INTERFACE_GENERATOR(Z,N,FNAME_AND_CALL) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  void  BOOST_PP_TUPLE_ELEM(3,0,FNAME_AND_CALL) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL) & CONTROL_PACKET) == 0) inc_calls_sent(target); \
    BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,FNAME_AND_CALL),N) \
        <T, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(this, dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL), target,obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \
  
  /*
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (remote_call, dc_impl::object_call_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (fast_remote_call,dc_impl::object_call_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (control_call,dc_impl::object_call_issue, (FAST_CALL | CONTROL_PACKET)) )
 

  /*
  The generation procedure for requests are the same. The only difference is that the function
  name has to be changed a little to be identify the return type of the function,
  (typename dc_impl::function_ret_type<__GLRPC_FRESULT>) and the issuing processor is object_request_issue.
  
    The call can be issued with
    \code
    ret = rmi.remote_request(target,
                              &object_type::function_name,
                              arg1,
                              arg2...)
    \endcode
  */
  #define REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
    BOOST_PP_TUPLE_ELEM(3,0,ARGS) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,ARGS) & CONTROL_PACKET) == 0) inc_calls_sent(target); \
    return BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,ARGS),N) \
        <T, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(this, dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,ARGS), target,obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \

  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<__GLRPC_FRESULT>::type remote_request, dc_impl::object_request_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<__GLRPC_FRESULT>::type fast_remote_request, dc_impl::object_request_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<__GLRPC_FRESULT>::type control_request, dc_impl::object_request_issue, (FAST_CALL | CONTROL_PACKET)) )
 


  #undef RPC_INTERFACE_GENERATOR
  #undef REQUEST_INTERFACE_GENERATOR
  
  /* Now generate the interface functions which allow me to call this dc_dist_object directly
  The internal calls are similar to the ones above. The only difference is that is that instead of
  'obj_id', the parameter passed to the issue processor is "control_obj_id" which identifies the
  current RMI class.
  */
  #define RPC_INTERFACE_GENERATOR(Z,N,FNAME_AND_CALL) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  void  BOOST_PP_TUPLE_ELEM(3,0,FNAME_AND_CALL) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL) & CONTROL_PACKET) == 0) inc_calls_sent(target); \
    BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,FNAME_AND_CALL),N) \
        <dc_dist_object<T>, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(this, dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL), target,control_obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \
  
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (internal_call,dc_impl::object_call_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (internal_fast_call,dc_impl::object_call_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (internal_control_call,dc_impl::object_call_issue, (FAST_CALL | CONTROL_PACKET)) )
 

  #define REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
    BOOST_PP_TUPLE_ELEM(3,0,ARGS) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, dc_.senders.size()); \
    if ((BOOST_PP_TUPLE_ELEM(3,2,ARGS) & CONTROL_PACKET) == 0) inc_calls_sent(target); \
    return BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,ARGS),N) \
        <dc_dist_object<T>, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(this, dc_.senders[target],  BOOST_PP_TUPLE_ELEM(3,2,ARGS), target,control_obj_id, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \

  /*
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<__GLRPC_FRESULT>::type internal_request, dc_impl::object_request_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<__GLRPC_FRESULT>::type internal_fast_request, dc_impl::object_request_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<__GLRPC_FRESULT>::type internal_control_request, dc_impl::object_request_issue, (FAST_CALL | CONTROL_PACKET)) )
 

  #undef RPC_INTERFACE_GENERATOR
  #undef REQUEST_INTERFACE_GENERATOR
  #undef GENARC
  #undef GENT
  #undef GENI
  #undef GENARGS
  


/*****************************************************************************
                      Implementation of matched send_to / recv_from
 *****************************************************************************/


 private:
  std::vector<dc_impl::recv_from_struct> recv_froms;
  
  void block_and_wait_for_recv(size_t src,
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
    strm.flush();
    dc_impl::reply_ret_type rt(REQUEST_WAIT_METHOD);
    // I shouldn't use a request to block here since 
    // that will take up a thread on the remote side
    // so I simulate a request here.
    size_t rtptr = reinterpret_cast<size_t>(&rt);
    if (control == false) {
      internal_call(target, &dc_dist_object<T>::block_and_wait_for_recv, 
                     procid(), strm.str(), rtptr);
    }
    else {
      internal_control_call(target, &dc_dist_object<T>::block_and_wait_for_recv, 
                  procid(), strm.str(), rtptr);
    }
    // wait for reply
    rt.wait();
    
    if (control == false) inc_calls_received(target);
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
      inc_calls_sent(source);
    }
    else {
      dc_.control_call(source, reply_increment_counter, tag, dc_impl::blob());
    }
  }





/*****************************************************************************
                      Implementation of Broadcast
 *****************************************************************************/
  
private:

  std::string broadcast_receive;

  void set_broadcast_receive(const std::string &s) {
    broadcast_receive = s;
  }


 public:
 
  /**
     This function allows one machine to broadcasts a variable to all machines.

     The originator calls broadcast with data provided in 
     in 'data' and originator set to true. 
     All other callers call with originator set to false.

     The originator will then return 'data'. All other machines
     will receive the originator's transmission in the "data" parameter.

     This call is guaranteed to have barrier-like behavior. That is to say,
     this call will block until all machines enter the broadcast function.

     \note Behavior is undefined if more than one machine calls broadcast
     with originator set to true.

     \note Behavior is undefined if multiple threads on the same machine
     call broadcast simultaneously. If multiple-thread broadcast is necessary,
     each thread should use its own instance of the services class.
  */
  template <typename U>
  void broadcast(U& data, bool originator, bool control = false) { 
    if (originator) {
      // construct the data stream
      std::stringstream strm;
      oarchive oarc(strm);
      oarc << data;
      strm.flush();
      broadcast_receive = strm.str();
      if (control == false) {
        for (size_t i = 0;i < numprocs(); ++i) {
          if (i != procid()) {
            internal_request(i,
                            &dc_dist_object<T>::set_broadcast_receive,
                            broadcast_receive);
          }
        }
      }
      else {
        for (size_t i = 0;i < numprocs(); ++i) {
          if (i != procid()) {
            internal_control_request(i,
                                    &dc_dist_object<T>::set_broadcast_receive,
                                    broadcast_receive);
          }
        }
      }
    }
  
    // by the time originator gets here, all machines
    // will have received the data due to the broadcast_receive
    // set a barrier here.
    barrier();
  
    // all machines will now deserialize the data
    if (!originator) {
      std::stringstream strm(broadcast_receive);
      iarchive iarc(strm);
      iarc >> data;
    }
  }


/*****************************************************************************
      Implementation of Gather, all_gather and gather_partition
 *****************************************************************************/

 private:
  std::vector<std::string> gather_receive;
  atomic<size_t> gatherid;
  
  void set_gather_receive(procid_t source, const std::string &s, size_t gid) {
    while(gatherid.value != gid) sched_yield();
    gather_receive[source] = s;
  }
 public:
  /**
   * Collects information contributed by each machine onto 
   * one machine.
   * data must be of length data[numprocs].
   * My data is stored in data[dc.procid()].
   * when function returns, machine sendto will have the complete vector
   * where data[i] is the data contributed by machine i.
   * All machines must have the same parameter for "sendto"
   */
  template <typename U>
  void gather(std::vector<U>& data, procid_t sendto, bool control = false) {
    // if not root
    if (sendto != procid()) {
      std::stringstream strm( std::ios::out | std::ios::binary );
      oarchive oarc(strm);
      oarc << data[procid()];
      strm.flush();
      if (control == false) {
        internal_request(sendto,
                        &dc_dist_object<T>::set_gather_receive,
                        procid(),
                        strm.str(),
                        gatherid.value);
      }
      else {
        internal_control_request(sendto,
                                  &dc_dist_object<T>::set_gather_receive,
                                  procid(),
                                  strm.str(),
                                  gatherid.value);
      }
    }
    barrier();
    if (sendto == procid()) {
      // if I am the receiver
      for (procid_t i = 0; i < numprocs(); ++i) {
        if (i != procid()) {
          // receiving only from others
          std::stringstream strm(gather_receive[i], 
                                 std::ios::in | std::ios::binary);
          assert(strm.good());
          iarchive iarc(strm);
          iarc >> data[i];
        }
      }
    }
    gatherid.inc();
  
  }

/********************************************************************
             Implementation of all gather
*********************************************************************/



 private:
  // ------- Sense reversing barrier data ----------
  /// The next value of the barrier. either +1 or -1
  int ab_barrier_sense;
  /// When this flag == the current barrier value. The barrier is complete
  int ab_barrier_release;
  /** when barrier sense is 1, barrier clears when
   * child_barrier_counter == numchild. When barrier sense is -1, barrier
   * clears when child_barrier_counter == 0;
   */
  atomic<int> ab_child_barrier_counter;
  /// condition variable and mutex protecting the barrier variables
  conditional ab_barrier_cond;
  mutex ab_barrier_mut;
  std::string ab_children_data[BARRIER_BRANCH_FACTOR];
  std::string ab_alldata;
  
  /**
    The child calls this function in the parent once the child enters the barrier
  */
  void __ab_child_to_parent_barrier_trigger(procid_t source, std::string collect) {
    ab_barrier_mut.lock();
    // assert childbase <= source <= childbase + BARRIER_BRANCH_FACTOR
    ASSERT_GE(source, childbase);
    ASSERT_LT(source, childbase + BARRIER_BRANCH_FACTOR);
    ab_children_data[source - childbase] = collect;
    ab_child_barrier_counter.inc(ab_barrier_sense);
    ab_barrier_cond.signal();
    ab_barrier_mut.unlock();
  }

  /**
    This is on the downward pass of the barrier. The parent calls this function
    to release all the children's barriers
  */
  void __ab_parent_to_child_barrier_release(int releaseval, 
                                            std::string allstrings,
                                            int use_control_calls) {
    // send the release downwards
    // get my largest child
    ab_alldata = allstrings;
    for (procid_t i = 0;i < numchild; ++i) {
      if (use_control_calls) {
        internal_control_call((procid_t)(childbase + i),
                              &dc_dist_object<T>::__ab_parent_to_child_barrier_release,
                              releaseval,
                              ab_alldata,
                              use_control_calls);
      }
      else {
        internal_call((procid_t)(childbase + i),
                      &dc_dist_object<T>::__ab_parent_to_child_barrier_release,
                      releaseval,
                      ab_alldata,
                      use_control_calls);
      }
    }
    ab_barrier_mut.lock();
    ab_barrier_release = releaseval;
    ab_barrier_cond.signal();
    ab_barrier_mut.unlock();
  }


 public:
   
  /**
   * Each machine creates a vector 'data' with size equivalent to the number of machines.
   * Each machine then fills the entry data[procid()] with information that it 
   * wishes to communicate.
   * After calling all_gather(), all machines will return with identical
   * vectors 'data', where data[i] contains the information machine i stored.
   */
  template <typename U>
  void all_gather(std::vector<U>& data, bool control = false) {
    if (numprocs() == 1) return;
    // get the string representation of the data
    charstream strm(128);
    oarchive oarc(strm);
    oarc << data[procid()];
    strm.flush();
    // upward message
    int ab_barrier_val = ab_barrier_sense;
    ab_barrier_mut.lock();
    // wait for all children to be done
    while(1) {
      if ((ab_barrier_sense == -1 && ab_child_barrier_counter.value == 0) ||
          (ab_barrier_sense == 1 && ab_child_barrier_counter.value == (int)(numchild))) {
        // flip the barrier sense
        ab_barrier_sense = -ab_barrier_sense;
        // call child to parent in parent
        ab_barrier_mut.unlock();
        if (procid() != 0) {
          // collect all my children data
          charstream strstrm(128);
          oarchive oarc2(strstrm);
          oarc2 << std::string(strm->c_str(), strm->size());
          for (procid_t i = 0;i < numchild; ++i) {
            strstrm.write(ab_children_data[i].c_str(), ab_children_data[i].length());
          }
          strstrm.flush();
          if (control) {
            internal_control_call(parent,
                            &dc_dist_object<T>::__ab_child_to_parent_barrier_trigger,
                            procid(),
                            std::string(strstrm->c_str(), strstrm->size()));
          }
          else {
            internal_call(parent,
                          &dc_dist_object<T>::__ab_child_to_parent_barrier_trigger,
                          procid(),
                          std::string(strstrm->c_str(), strstrm->size()));
          }
        }
        break;
      }
      ab_barrier_cond.wait(ab_barrier_mut);
    }


    //logger(LOG_DEBUG, "barrier phase 1 complete");
    // I am root. send the barrier release downwards
    if (procid() == 0) {
      ab_barrier_release = ab_barrier_val;
      // build the downward data
      charstream strstrm(128);
      oarchive oarc2(strstrm);
      oarc2 << std::string(strm->c_str(), strm->size());
      for (procid_t i = 0;i < numchild; ++i) {
        strstrm.write(ab_children_data[i].c_str(), ab_children_data[i].length());
      }
      strstrm.flush();
      ab_alldata = std::string(strstrm->c_str(), strstrm->size());
      for (procid_t i = 0;i < numchild; ++i) {
        internal_control_call((procid_t)(childbase + i),
                             &dc_dist_object<T>::__ab_parent_to_child_barrier_release,
                             ab_barrier_val,
                             ab_alldata,
                             (int)control);

      }
    }
    // wait for the downward message releasing the barrier
    ab_barrier_mut.lock();
    while(1) {
      if (ab_barrier_release == ab_barrier_val) break;
      ab_barrier_cond.wait(ab_barrier_mut);
    }
    // read the collected data and release the lock
    std::string local_ab_alldata = ab_alldata;
    ab_barrier_mut.unlock();

    //logger(LOG_DEBUG, "barrier phase 2 complete");
    // now the data is a DFS search of a heap
    // I need to unpack it
    size_t heappos = 0;
    std::stringstream istrm(local_ab_alldata);
    iarchive iarc(istrm);

    for (size_t i = 0;i < numprocs(); ++i) {
      std::string s;
      iarc >> s;

      std::stringstream strm2(s);
      iarchive iarc2(strm2);
      iarc2 >> data[heappos];
      
      if (i + 1 == numprocs()) break;
      // advance heappos
      // leftbranch
      bool lefttraverseblock = false;
      while (1) {
        // can we continue going deaper down the left?
        size_t leftbranch = heappos * BARRIER_BRANCH_FACTOR + 1;
        if (lefttraverseblock == false && leftbranch < numprocs()) {
          heappos = leftbranch;
          break;
        }
        // ok. can't go down the left
        bool this_is_a_right_branch = (((heappos - 1) % BARRIER_BRANCH_FACTOR) == BARRIER_BRANCH_FACTOR - 1);
        // if we are a left branch, go to sibling
        if (this_is_a_right_branch == false) {
          size_t sibling = heappos + 1;
          if (sibling < numprocs()) {
            heappos = sibling;
            break;
          }
        }

        // we have finished this subtree, go back up to parent
        // and block the depth traversal on the next round
        // unless heappos is 0

        heappos = (heappos - 1) / BARRIER_BRANCH_FACTOR;
        lefttraverseblock = true;
        continue;
        // go to sibling
      }
      
    }
  }
  



////////////////////////////////////////////////////////////////////////////

  /**
   * This function is takes a vector of local elements T which must
   * be comparable and constructs a vector of length numprocs where
   * each element is a subset of the local contribution from that
   * machine and the union of all elements in the union of all local
   * contributions and all entries are unique:
   *
   * Usage: Each process reads the files that are stored locally and
   * wants to know which subset of local files to read even when
   * multiple processes see the same files.
   */
  template <typename U>
  void gather_partition(const std::vector<U>& local_contribution,
                        std::vector< std::vector<U> >& ret_partition,
                        bool control = false) {
    typedef std::set<U> set_type;

    // Compute the elements on each machine
    std::vector< std::set<U> > cpu2elems(numprocs());
    cpu2elems[procid()].insert(local_contribution.begin(), 
                                   local_contribution.end());

    gather(cpu2elems, 0);
    // Construct the "balanced" partitioning
    if(procid() == 0) {
      ret_partition.clear();
      ret_partition.resize(numprocs());
      // Construct the union
      std::set<U> unassigned_elems;
      foreach(const set_type& set, cpu2elems) 
        unassigned_elems.insert(set.begin(), set.end());
      // Assign elements to each of the machines      
      for(procid_t cpuid = 0; !unassigned_elems.empty(); 
          cpuid = (cpuid + 1) % cpu2elems.size()) {
        // while there are things left to be assigned to this cpu
        while( !cpu2elems[cpuid].empty() ) {
          // Get the next element and remove it
          U elem = *(cpu2elems[cpuid].begin());
          cpu2elems[cpuid].erase(cpu2elems[cpuid].begin());
          // if the next element on this cpu is not yet assigned then
          // assign it to this cpu
          if(unassigned_elems.count(elem) != 0) {
            unassigned_elems.erase(elem);
            ret_partition[cpuid].push_back(elem);
            break;
          }
        
        } // end of while loop
      } // end of loop over cpus
      assert(unassigned_elems.empty());
    }
    // Scatter the result
    broadcast(ret_partition, procid() == 0, control);    
  } // end of gather_partition






/*****************************************************************************
                      Implementation of Barrier
 *****************************************************************************/



 private:
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
  procid_t parent;  /// parent node
  size_t childbase; /// id of my first child
  procid_t numchild;  /// number of children




  /**
    The child calls this function in the parent once the child enters the barrier
  */
  void __child_to_parent_barrier_trigger(procid_t source) {
    barrier_mut.lock();
    // assert childbase <= source <= childbase + BARRIER_BRANCH_FACTOR
    ASSERT_GE(source, childbase);
    ASSERT_LT(source, childbase + BARRIER_BRANCH_FACTOR);
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
    for (procid_t i = 0;i < numchild; ++i) {
      internal_control_call((procid_t)(childbase + i),
                            &dc_dist_object<T>::__parent_to_child_barrier_release,
                            releaseval);
  
    }
    barrier_mut.lock();
    barrier_release = releaseval;    
    barrier_cond.signal();
    barrier_mut.unlock();
  }
  

 public:
  /**
    A regular barrier equivalent to MPI_Barrier.
    A machine entering this barrier will wait until every machine 
    reaches this barrier before continuing. Only one thread from each machine
    should call the barrier.
    
    \see full_barrier
    */
  void barrier() {
    // upward message
    int barrier_val = barrier_sense;      
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
    
    
    //logger(LOG_DEBUG, "barrier phase 1 complete");
    // I am root. send the barrier release downwards
    if (procid() == 0) {
      barrier_release = barrier_val;
  
      for (procid_t i = 0;i < numchild; ++i) {
        internal_control_call((procid_t)(childbase + i),
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
  
    //logger(LOG_DEBUG, "barrier phase 2 complete");
  }
  
  
 /*****************************************************************************
                      Implementation of Full Barrier
*****************************************************************************/
 private:
  mutex full_barrier_lock;
  conditional full_barrier_cond;
  std::vector<size_t> calls_to_receive;
  // used to inform the counter that the full barrier
  // is in effect and all modifications to the calls_recv
  // counter will need to lock and signal
  bool full_barrier_in_effect;
  
  /** number of 'source' processor counts which have
  not achieved the right recv count */
  atomic<size_t> num_proc_recvs_incomplete; 
                                      
  /// Marked as 1 if the proc is complete
  dense_bitset procs_complete;

 public:  
  /**
  Similar to the barrier(), but provides additional guarantees that 
  all RMI calls issued prior to this barrier are completed before
  returning. 
  
  \note This function could return prematurely if
  other threads are still issuing function calls since we
  cannot differentiate between calls issued before the barrier
  and calls issued while the barrier is being evaluated.
  Therefore, when used in a multithreaded scenario, the user must ensure
  that all other threads which may perform operations using this object
  are stopped before the full barrier is initated.
  
  \see barrier
  */
  void full_barrier() {
    // gather a sum of all the calls issued to machine 0
    std::vector<size_t> calls_sent_to_target(numprocs(), 0);
    for (size_t i = 0;i < numprocs(); ++i) {
      calls_sent_to_target[i] = callssent[i].value;
    }
    
    // tell node 0 how many calls there are
    std::vector<std::vector<size_t> > all_calls_sent(numprocs());
    all_calls_sent[procid()] = calls_sent_to_target;
    all_gather(all_calls_sent, true);
    
    // get the number of calls I am supposed to receive from each machine
    calls_to_receive.clear(); calls_to_receive.resize(numprocs(), 0);
    for (size_t i = 0;i < numprocs(); ++i) {
      calls_to_receive[i] += all_calls_sent[i][procid()];
//      std::cout << "Expecting " << calls_to_receive[i] << " calls from " << i << std::endl;
    }
    // clear the counters
    num_proc_recvs_incomplete.value = numprocs();
    procs_complete.clear();
    // activate the full barrier
    full_barrier_in_effect = true;
    // begin one pass to set all which are already completed
    for (size_t i = 0;i < numprocs(); ++i) {
      if (callsreceived[i].value >= calls_to_receive[i]) {
        if (procs_complete.set_bit((uint32_t)i) == false) {
          num_proc_recvs_incomplete.dec();
        }
      }
    }
    
    full_barrier_lock.lock();
    while (num_proc_recvs_incomplete.value > 0) full_barrier_cond.wait(full_barrier_lock);
    full_barrier_lock.unlock();
    full_barrier_in_effect = false;
//     for (size_t i = 0; i < numprocs(); ++i) {
//       std::cout << "Received " << global_calls_received[i].value << " from " << i << std::endl;
//     }
    barrier();
  }
  
 /* --------------------  Implementation of Gather Statistics -----------------*/ 
 private:
  struct collected_statistics {
    size_t callssent;
    size_t bytessent;
    collected_statistics(): callssent(0), bytessent(0) { }
    void save(oarchive &oarc) const {
      oarc << callssent << bytessent;
    }
    void load(iarchive &iarc) {
      iarc >> callssent >> bytessent;
    }
  };
 public:
  /** Gather RPC statistics. All machines must call 
   this function at the same time. However, only proc 0 will
   return values */
  std::map<std::string, size_t> gather_statistics() {
    std::map<std::string, size_t> ret;

    std::vector<collected_statistics> stats(numprocs());
    stats[procid()].callssent = calls_sent();
    stats[procid()].bytessent = bytes_sent();
    logstream(LOG_INFO) << procid() << ": calls_sent: ";
    for (size_t i = 0;i < numprocs(); ++i) {
      logstream(LOG_INFO) << callssent[i].value << ", ";
    }
    logstream(LOG_INFO) << std::endl;
    logstream(LOG_INFO) << procid() << ": calls_recv: ";
    for (size_t i = 0;i < numprocs(); ++i) {
      logstream(LOG_INFO) << callsreceived[i].value << ", ";
    }
    logstream(LOG_INFO) << std::endl;


    gather(stats, 0, true);
    if (procid() == 0) {
      collected_statistics cs;
      for (size_t i = 0;i < numprocs(); ++i) {
        cs.callssent += stats[i].callssent;
        cs.bytessent += stats[i].bytessent;
      }
      ret["total_calls_sent"] = cs.callssent;
      ret["total_bytes_sent"] = cs.bytessent;
    }
    return ret;
  }
};

#include <graphlab/macros_undef.hpp>
#include <graphlab/rpc/mem_function_arg_types_undef.hpp>
#undef BARRIER_BRANCH_FACTOR
}// namespace graphlab
#endif

