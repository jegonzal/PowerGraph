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

#include <graphlab/rpc/dc_dist_object.hpp>
#ifndef GRAPHLAB_DC_SERVICES_HPP
#define GRAPHLAB_DC_SERVICES_HPP
#include <graphlab/parallel/pthread_tools.hpp>



#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
    \ingroup rpc
    Creates a new context for MPI-like global global operations.
    Where all machines create an instance of dc_services at the same time,
    operations performed by the new dc_services instance will not interfere
    and will run in parallel with other contexts. 
    i.e. If I have two distributed dc_services instances, one instance can 
    perform a barrier while another instance performs a broadcast() at the same 
    time.
    
    \note Only simple algorithms for the MPI collective operations (barrier, broadcast, etc)
    are implemented. Significant work is necessary to improve the performance of the collectives.
  */
  class dc_services {
  private:
    dc_dist_object<dc_services> rmi;
  
  public:
    dc_services(distributed_control &dc):rmi(dc, this) {  }
    
    /// Returns the underlying dc_dist_object 
    dc_dist_object<dc_services>& rmi_instance() {
      return rmi;
    }

    /// Returns the underlying dc_dist_object 
    const dc_dist_object<dc_services>& rmi_instance() const {
      return rmi;
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
      rmi.comm_barrier(targetmachine);
    }

  /**
    This is a convenience function which broadcasts a comm_barrier()
    \note having all machines call the comm barrier does not guarantee
    that all calls have been processed. Basically 'p' local barriers
    do not result in a global barrier.
    
    \note This affects the global context
  */
    inline void comm_barrier() {
      rmi.comm_barrier();
    }
    
    /**
    This is a blocking send_to. It send an object T to the target 
    machine, but waits for the target machine to call recv_from
    before returning. Functionally similar to MPI's matched sending/receiving
    */
    template <typename U>
    inline void send_to(procid_t target, U& t, bool control = false) {
      rmi.send_to(target, t, control);
    }
    
    /**
    A blocking recv_from. Must be matched with a send_to call from the
    target before both source and target resumes.
    */
    template <typename U>
    inline void recv_from(procid_t source, U& t, bool control = false) {
      rmi.recv_from(source, t, control);
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
    template <typename U>
    inline void broadcast(U& data, bool originator, bool control = false) { 
      rmi.broadcast(data, originator, control);
    }

  /**
   * data must be of length data[numprocs].
   * My data is stored in data[dc.procid()].
   * when function returns, machine sendto will have the complete vector
   * where data[i] is the data contributed by machine i.
   * All machines must have the same parameter for "sendto"
   */
    template <typename U>
    inline void gather(std::vector<U>& data, procid_t sendto, bool control = false) {
      rmi.gather(data, sendto, control);
    }

  /**
   * Each machine creates a vector 'data' with size equivalent to the number of machines.
   * Each machine then fills the entry data[procid()] with information that it 
   * wishes to communicate.
   * After calling all_gather(), all machines will return with identical
   * vectors 'data', where data[i] contains the information machine i stored.
   */
    template <typename U>
    inline void all_gather(std::vector<U>& data, bool control = false) {
      rmi.all_gather(data, control);
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
    inline void gather_partition(const std::vector<U>& local_contribution,
                          std::vector< std::vector<U> >& ret_partition,
                          bool control = false) {
      rmi.gather_partition(local_contribution, ret_partition, control);
    }
    
    /**
    A regular barrier equivalent to MPI_Barrier.
    A thread machine entering this barrier will wait until one thread on each 
    machines enter this barrier.
    
    \see full_barrier
    */
    inline void barrier() {
      rmi.barrier();
    }
    
    
    /**
  This barrier ensures globally across all machines that
  all calls issued prior to this barrier are completed before
  returning. This function could return prematurely if
  other threads are still issuing function calls since we
  cannot differentiate between calls issued before the barrier
  and calls issued while the barrier is being evaluated.
  
  Therefore, when used in a multithreaded scenario, the user must ensure
  that all other threads which may perform operations using this object
  are stopped before the full barrier is initated.
  
  \see barrier
  */
    inline void full_barrier() {
      rmi.full_barrier();
    }
  
 

  };


} // end of namespace graphlab


#include <graphlab/macros_undef.hpp>
#endif

