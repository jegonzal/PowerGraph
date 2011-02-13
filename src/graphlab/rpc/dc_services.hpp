#include <graphlab/rpc/dc_dist_object.hpp>
#ifndef GRAPHLAB_DC_SERVICES_HPP
#define GRAPHLAB_DC_SERVICES_HPP
#include <graphlab/parallel/pthread_tools.hpp>



#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
    A thin wrapper around the dc_dist_object.
    Provides a non-templated reference to a "context"
  */
  class dc_services {
  private:
    dc_dist_object<dc_services> rmi;
  
  public:
    dc_services(distributed_control &dc):rmi(dc, this) {  }
    
    dc_dist_object<dc_services>& rmi_instance() {
      return rmi;
    }

    const dc_dist_object<dc_services>& rmi_instance() const {
      return rmi;
    }
  
    inline void comm_barrier(procid_t targetmachine) {
      rmi.comm_barrier(targetmachine);
    }
  
    inline void comm_barrier() {
      rmi.comm_barrier();
    }
  
    template <typename U>
    inline void send_to(procid_t target, U& t, bool control = false) {
      rmi.send_to(target, t, control);
    }
    
    template <typename U>
    inline void recv_from(procid_t source, U& t, bool control = false) {
      rmi.recv_from(source, t, control);
    }

    template <typename U>
    inline void broadcast(U& data, bool originator, bool control = false) { 
      rmi.broadcast(data, originator, control);
    }

    template <typename U>
    inline void gather(std::vector<U>& data, procid_t sendto, bool control = false) {
      rmi.gather(data, sendto, control);
    }

    template <typename U>
    inline void all_gather(std::vector<U>& data, bool control = false) {
      rmi.all_gather(data, control);
    }

    template <typename U>
    inline void gather_partition(const std::vector<U>& local_contribution,
                          std::vector< std::vector<U> >& ret_partition,
                          bool control = false) {
      rmi.gather_partition(local_contribution, ret_partition, control);
    }
    

    inline void barrier() {
      rmi.barrier();
    }
    
    inline void full_barrier() {
      rmi.full_barrier();
    }
  
 

  };


} // end of namespace graphlab


#include <graphlab/macros_undef.hpp>
#endif
