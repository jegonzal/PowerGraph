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
#include <boost/preprocessor.hpp>
#include <graphlab/macros_def.hpp>

#define BARRIER_BRANCH_FACTOR 128


namespace graphlab {


/**
Provides capabilities for distributed objects Your class should either
inherit this, or instantiate it before any distributed object call.
The requirement for using the distributed object is that all machines
must construct the distributed objects in the same order. And, no
distributed object calls should be make until it is guaranteed that
all machines have constructed their respective distributed objects.

This class also acts as a single "xontext" spanning multiple machines.
For instance, the barrier implemented here is fully localized to within 
a particular instance of this object. That is, multiple instances of the object
can issue multiple independent barriers.
The dc_services() object is a thin wrapper around the dc_dist_object.
*/
template <typename T>
class dc_dist_object : public dc_impl::dc_dist_object_base{
 private:
  distributed_control &dc_;
  size_t obj_id;
  size_t control_obj_id;  // object id of this object
  T* owner;
  bool calltracking;
  atomic<size_t> callsreceived;
  atomic<size_t> callssent;
  // make operator= private
  dc_dist_object<T>& operator=(const dc_dist_object<T> &d) {return *this;}
  friend class distributed_control;
  
  
  


 public:
  //internal stuff which should not be used.
  void inc_calls_received() {
    if (full_barrier_in_effect == false) {
      callsreceived.inc();
    }
    else {
      // use the more costly option
      full_barrier_lock.lock();
      callsreceived.inc();
      full_barrier_cond.signal();
      full_barrier_lock.unlock();
    }
  }
  void inc_calls_sent() {
    callssent.inc();
  }

 public:
  dc_dist_object(distributed_control &dc_, T* owner, bool calltracking = false):
                                dc_(dc_),owner(owner),calltracking(calltracking) {
    obj_id = dc_.register_object(owner, this);
    control_obj_id = dc_.register_object(this, this);
    
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
      numchild = maxchild - childbase;
    }
  
    parent =  (dc_.procid() - 1) / BARRIER_BRANCH_FACTOR   ;

    //-------- Initialize the full barrier ---------
    full_barrier_in_effect = false;
    full_barrier_curid = 0;
    full_barrier_released = false;

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
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (control_call,dc_impl::object_call_issue, (FAST_CALL | CONTROL_PACKET)) )
 

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
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type control_request, dc_impl::object_request_issue, (FAST_CALL | CONTROL_PACKET)) )
 


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
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (internal_control_call,dc_impl::object_call_issue, (FAST_CALL | CONTROL_PACKET)) )
 

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
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type internal_control_request, dc_impl::object_request_issue, (FAST_CALL | CONTROL_PACKET)) )
 

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

  void set_gather_receive(procid_t source, const std::string &s) {
    gather_receive[source] = s;
  }
 public:
  /**
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
      if (control == false) {
        internal_request(sendto,
                        &dc_dist_object<T>::set_gather_receive,
                        procid(),
                        strm.str());
      }
      else {
        internal_control_request(sendto,
                                  &dc_dist_object<T>::set_gather_receive,
                                  procid(),
                                  strm.str());
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
  
  }

  /**
   * data must be of length data[numprocs].
   * My data is stored in data[dc.procid()]
   * when function returns, everyone will have the same data vector
   * where data[i] is the data contributed by machine i.
   */
  template <typename U>
  void all_gather(std::vector<U>& data, bool control = false) {
    gather(data, 0, control);
    broadcast(data, procid() == 0, control);
  }



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
  size_t parent;  /// parent node
  size_t childbase; /// id of my first child
  size_t numchild;  /// number of children




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
  

 public:
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
    
    
    logger(LOG_DEBUG, "barrier phase 1 complete");
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
  
    logger(LOG_DEBUG, "barrier phase 2 complete");
  }
  
  
 /*****************************************************************************
                      Implementation of Full Barrier
*****************************************************************************/
 private:
  bool full_barrier_released;
  mutex full_barrier_lock;
  conditional full_barrier_cond;
  atomic<size_t> all_recv_count;
  size_t all_send_count;
  size_t full_barrier_curid; // to protect against fast repeated calls to full_barrier
  
  // used to inform the counter that the full barrier
  // is in effect and all modifications to the calls_recv
  // counter will need to lock and signal
  bool full_barrier_in_effect;
  
  void release_full_barrier(size_t id) {
    if (id != full_barrier_curid) return;
    full_barrier_lock.lock();
    full_barrier_released = true;
    full_barrier_in_effect = false;
    full_barrier_cond.signal();
    full_barrier_lock.unlock();
  }
  
  void full_barrier_add_to_recv(size_t id, size_t r) {
      if (id != full_barrier_curid) return;
    // we want the previous value of the atom
    // so we can find the first time it crosses
    // the send counter, and avoid multiple releases
    size_t prevval = all_recv_count.inc_ret_last(r);
    if (prevval < all_send_count && prevval + r >= all_send_count) {
      // release myself
      release_full_barrier(full_barrier_curid);
      // release everyone
      for (size_t i = 0;i < numprocs(); ++i) {
        if (i != procid()) {
          internal_control_call(i,
                               &dc_dist_object<T>::release_full_barrier,
                               full_barrier_curid);
        }
      }
    }
  }
 public:  
  /**
  This barrier ensures globally across all machines that
  all calls issued prior to this barrier are completed before
  returning. This function could return prematurely if
  other threads are still issuing function calls since we
  cannot differentiate between calls issued before the barrier
  and calls issued while the barrier is being evaluated.
  */
  void full_barrier() {
    // gather a sum of all the calls issued to machine 0
    size_t local_num_calls_sent = calls_sent();
    
    // tell node 0 how many calls there are
    std::vector<size_t> all_calls_sent(numprocs());
    all_calls_sent[procid()] = local_num_calls_sent;
    gather(all_calls_sent, 0, true);
    
    // proc 0 computes the total number of calls sent    
    all_send_count = 0;
    all_recv_count.value = 0;
    if (procid() == 0) {
      for (size_t i = 0;i < all_calls_sent.size(); ++i) {
        all_send_count += all_calls_sent[i];
      }
    }
    // issue a barrier to make sure everyone stops here
    // while node 0 prepares the counters.
    barrier();
    // ok. now we basically keep recomunicating with
    // node 0 the number of calls we have received so far
    // until node 0 releases the barrier
    full_barrier_lock.lock();
    full_barrier_in_effect = true;
    size_t last_communicated_recv_count = 0;
    
    
    while(full_barrier_released == false) {
      while (calls_received() != last_communicated_recv_count) {
        size_t nextval = calls_received();
        assert(nextval > last_communicated_recv_count);
        // unlock for a while to issue the RPC call
        full_barrier_lock.unlock();
        if (procid() == 0) {
          full_barrier_add_to_recv(full_barrier_curid, 
                                   nextval - last_communicated_recv_count);
        }
        else {
          internal_control_call(0,
                             &dc_dist_object<T>::full_barrier_add_to_recv,
                             full_barrier_curid, 
                             nextval - last_communicated_recv_count);
        }
        last_communicated_recv_count = nextval;
        full_barrier_lock.lock();
        // now there could be a race here because I released
        // the lock to issue the calls.
        // I need to check again before I let the condition
        // variable take over
        // the inner while loop will check the counting case
        // but I need to check the exterior case
      }
      if (full_barrier_released) break;
      full_barrier_cond.wait(full_barrier_lock);    
    }
    full_barrier_curid++;
    full_barrier_lock.unlock();
  }
};

#include <graphlab/macros_undef.hpp>
#include <graphlab/rpc/mem_function_arg_types_undef.hpp>
#undef BARRIER_BRANCH_FACTOR
}// namespace graphlab
#endif
