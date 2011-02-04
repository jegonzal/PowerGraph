#ifndef DC_SERVICES_HPP
#define DC_SERVICES_HPP
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>

#define DC_SERVICES_BARRIER_BRANCH_FACTOR 128
namespace graphlab {

class dc_services {
 private:
  dc_dist_object<dc_services> rpc;
  
  // Sense reversing barrier
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
  
  std::string broadcast_receive;

  std::vector<std::string> gather_receive;
  // set the waiting flag
  void __child_to_parent_barrier_trigger(procid_t source);
  
  void __parent_to_child_barrier_release(int releaseval);
  
  
  void set_broadcast_receive(const std::string &s) {
    broadcast_receive = s;
  }

  void set_gather_receive(procid_t source, const std::string &s) {
    gather_receive[source] = s;
  }
  
 public:
  dc_services(distributed_control &dc):rpc(dc, this) { 
    // reset the child barrier values
    child_barrier_counter.value = 0;
    barrier_sense = 1;
    barrier_release = -1;

    // compute my children
    childbase = size_t(rpc.dc().procid()) * DC_SERVICES_BARRIER_BRANCH_FACTOR + 1;
    if (childbase >= rpc.dc().numprocs()) {
      numchild = 0;
    }
    else {
      size_t maxchild = std::min<size_t>(rpc.dc().numprocs(), 
                               childbase + DC_SERVICES_BARRIER_BRANCH_FACTOR);
      numchild = maxchild - childbase;
    }
    
    parent =  (rpc.dc().procid() - 1) / DC_SERVICES_BARRIER_BRANCH_FACTOR   ;
    logger(LOG_INFO, "%d %d %d", childbase, numchild, parent);
  }
 
  /**
   * tree-reduction based sense-reversing barrier with a branching factor of 
   * DC_SERVICES_BARRIER_BRANCH_FACTOR
   */
  void barrier();
  
  
  inline void comm_barrier(procid_t targetmachine) {
    rpc.dc().comm_barrier(targetmachine);
  }
  
  inline void comm_barrier() {
    rpc.dc().comm_barrier();
  }
  
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
  template <typename T>
  void broadcast(T& data, bool originator) { 
    if (originator) {
      // construct the data stream
      std::stringstream strm;
      oarchive oarc(strm);
      oarc << data;
      broadcast_receive = strm.str();
      
      for (size_t i = 0;i < rpc.numprocs(); ++i) {
        if (i != rpc.procid()) {
          rpc.fast_remote_request(i,
                                  &dc_services::set_broadcast_receive,
                                  broadcast_receive);
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

  /**
   * data must be of length data[numprocs].
   * My data is stored in data[dc.procid()].
   * when function returns, machine sendto will have the complete vector
   * where data[i] is the data contributed by machine i.
   * All machines must have the same parameter for "sendto"
   */
  template <typename T>
  void gather(std::vector<T>& data, procid_t sendto) {
    if (sendto != rpc.procid()) {
      std::stringstream strm;
      oarchive oarc(strm);
      oarc << data;
      gather_receive.resize(rpc.numprocs());
      rpc.fast_remote_request(sendto,
                              &dc_services::set_gather_receive,
                              rpc.procid(),
                              strm.str());
    }
    barrier();

    if (sendto == rpc.procid()) {
      for (procid_t i = 0; i < rpc.numprocs(); ++i) {
        std::stringstream strm(gather_receive[i]);
        iarchive iarc(strm);
        iarc >> data[i];
      }
    }
    
  }

  /**
   * data must be of length data[numprocs].
   * My data is stored in data[dc.procid()]
   * when function returns, everyone will have the same data vector
   * where data[i] is the data contributed by machine i.
   */
  template <typename T>
  void all_gather(std::vector<T>& data) {
    gather(data, 0);
    broadcast(data, rpc.procid() == 0);
  }
};

}
#endif
