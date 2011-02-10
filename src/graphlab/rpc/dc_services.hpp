#ifndef GRAPHLAB_DC_SERVICES_HPP
#define GRAPHLAB_DC_SERVICES_HPP
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  class dc_services {
  private:
    dc_dist_object<dc_services> rpc;
 
    std::string broadcast_receive;

    std::vector<std::string> gather_receive;

  
    void set_broadcast_receive(const std::string &s) {
      broadcast_receive = s;
    }

    void set_gather_receive(procid_t source, const std::string &s) {
      gather_receive[source] = s;
    }
  
  public:
    dc_services(distributed_control &dc):rpc(dc, this) { 
      // initialize gathers
      gather_receive.resize(rpc.numprocs());

    }
    
    dc_dist_object<dc_services>& rmi_instance() {
      return rpc;
    }
 
    /**
     * tree-reduction based sense-reversing barrier with a branching factor of 
     * DC_SERVICES_BARRIER_BRANCH_FACTOR
     */
    inline void barrier() {
      rpc.barrier();
    }
  
  
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
      // if not root
      if (sendto != rpc.procid()) {
        std::stringstream strm( std::ios::out | std::ios::binary );
        oarchive oarc(strm);
        oarc << data[rpc.procid()];
        rpc.fast_remote_request(sendto,
                                &dc_services::set_gather_receive,
                                rpc.procid(),
                                strm.str());
      }
      barrier();
      if (sendto == rpc.procid()) {
        // if I am the receiver
        for (procid_t i = 0; i < rpc.numprocs(); ++i) {
          if (i != rpc.procid()) {
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
    template <typename T>
    void all_gather(std::vector<T>& data) {
      gather(data, 0);
      broadcast(data, rpc.procid() == 0);
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
    template <typename T>
    void gather_partition(const std::vector<T>& local_contribution,
                          std::vector< std::vector<T> >& ret_partition) {
      typedef std::set<T> set_type;

      // Compute the elements on each machine
      std::vector< std::set<T> > cpu2elems(rpc.numprocs());
      cpu2elems[rpc.procid()].insert(local_contribution.begin(), 
                                     local_contribution.end());

      gather(cpu2elems, 0);
      // Construct the "balanced" partitioning
      if(rpc.procid() == 0) {
        ret_partition.clear();
        ret_partition.resize(rpc.numprocs());
        // Construct the union
        std::set< T > unassigned_elems;
        foreach(const set_type& set, cpu2elems) 
          unassigned_elems.insert(set.begin(), set.end());
        // Assign elements to each of the machines      
        for(procid_t cpuid = 0; !unassigned_elems.empty(); 
            cpuid = (cpuid + 1) % cpu2elems.size()) {
          // while there are things left to be assigned to this cpu
          while( !cpu2elems[cpuid].empty() ) {
            // Get the next element and remove it
            T elem = *(cpu2elems[cpuid].begin());
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
      broadcast(ret_partition, rpc.procid() == 0);    
    } // end of gather_partition






  };


} // end of namespace graphlab


#include <graphlab/macros_undef.hpp>
#endif
