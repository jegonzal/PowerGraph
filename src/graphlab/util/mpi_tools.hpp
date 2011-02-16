#ifndef GRAPHLAB_MPI_TOOLS
#define GRAPHLAB_MPI_TOOLS

#include <sys/types.h>
#include <ifaddrs.h>
#include <netinet/in.h>

#include <mpi.h>

#include <vector>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/charstream.hpp>






#include <graphlab/macros_def.hpp>

namespace graphlab {
  namespace mpi_tools {



    inline void init(int& argc, char**& argv) {
      int error = MPI_Init(&argc, &argv);
      assert(error == MPI_SUCCESS);
    }

    inline void finalize() {
      int error = MPI_Finalize();
      assert(error == MPI_SUCCESS);
    }

    inline size_t rank() {
      int mpi_rank(-1);
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      assert(mpi_rank >= 0);
      return size_t(mpi_rank);
    }

    inline size_t size() {
      int mpi_size(-1);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      assert(mpi_size >= 0);
      return size_t(mpi_size);
    }

    

    template<typename T>
    void all_gather(const T& elem, std::vector<T>& results) {
      // Get the mpi rank and size
      size_t mpi_size(size());
      if(results.size() != mpi_size) results.resize(mpi_size);

      // Serialize the local map
      graphlab::charstream cstrm(128);
      graphlab::oarchive oarc(cstrm);
      oarc << elem;
      cstrm.flush();
      char* send_buffer = cstrm->c_str();
      int send_buffer_size = cstrm->size();
      assert(send_buffer_size >= 0);

      // compute the sizes
      std::vector<int> recv_sizes(mpi_size, -1);
      // Compute the sizes
      int error = MPI_Allgather(&send_buffer_size,  // Send buffer
                                1,                  // send count
                                MPI_INT,            // send type
                                &(recv_sizes[0]),  // recvbuffer
                                1,                  // recvcount
                                MPI_INT,           // recvtype
                                MPI_COMM_WORLD);  
      assert(error == MPI_SUCCESS);
      for(size_t i = 0; i < recv_sizes.size(); ++i)
        assert(recv_sizes[i] >= 0);
  

      // Construct offsets
      std::vector<int> recv_offsets(recv_sizes);
      int sum = 0, tmp = 0;
      for(size_t i = 0; i < recv_offsets.size(); ++i) {
        tmp = recv_offsets[i]; recv_offsets[i] = sum; sum += tmp; 
      }

      // if necessary realloac recv_buffer
      std::vector<char> recv_buffer(sum);
     
      // recv all the maps 
      error = MPI_Allgatherv(send_buffer,         // send buffer
                             send_buffer_size,    // how much to send
                             MPI_BYTE,            // send type
                             &(recv_buffer[0]),   // recv buffer
                             &(recv_sizes[0]),    // amount to recv
                                                  // for each cpuess
                             &(recv_offsets[0]),  // where to place data
                             MPI_BYTE,
                             MPI_COMM_WORLD);
      assert(error == MPI_SUCCESS);
      // Update the local map
      namespace bio = boost::iostreams;
      typedef bio::stream<bio::array_source> icharstream;
      icharstream strm(&(recv_buffer[0]), recv_buffer.size());
      graphlab::iarchive iarc(strm);
      for(size_t i = 0; i < results.size(); ++i) {
        iarc >> results[i];
      }  
    }







    /**
     * called on the root.  must be matched with gather(const T& elem);
     */
    template<typename T>
    void gather(size_t root, const T& elem) {
       // Get the mpi rank and size
      assert(root < size_t(std::numeric_limits<int>::max()));
      int mpi_root(root);

      // Serialize the local map
      graphlab::charstream cstrm(128);
      graphlab::oarchive oarc(cstrm);
      oarc << elem;
      cstrm.flush();
      char* send_buffer = cstrm->c_str();
      int send_buffer_size = cstrm->size();
      assert(send_buffer_size >= 0);

      // compute the sizes
      // Compute the sizes
      int error = MPI_Gather(&send_buffer_size,  // Send buffer
                             1,                  // send count
                             MPI_INT,            // send type
                             NULL,               // recvbuffer
                             1,                  // recvcount
                             MPI_INT,           // recvtype
                             mpi_root,          // root rank
                             MPI_COMM_WORLD);  
      assert(error == MPI_SUCCESS); 

     
      // recv all the maps 
      error = MPI_Gatherv(send_buffer,         // send buffer
                          send_buffer_size,    // how much to send
                          MPI_BYTE,            // send type
                          NULL,                // recv buffer
                          NULL,                // amount to recv
                                               // for each cpuess
                          NULL,                // where to place data
                          MPI_BYTE,
                          mpi_root,            // root rank
                          MPI_COMM_WORLD);
      assert(error == MPI_SUCCESS);
    } // end of gather






    /**
     * called on the root.  must be matched with gather(const T& elem);
     */
    template<typename T>
    void gather(const T& elem, std::vector<T>& results) {
      // Get the mpi rank and size
      size_t mpi_size(size());
      int mpi_rank(rank());
      if(results.size() != mpi_size) results.resize(mpi_size);

      // Serialize the local map
      graphlab::charstream cstrm(128);
      graphlab::oarchive oarc(cstrm);
      oarc << elem;
      cstrm.flush();
      char* send_buffer = cstrm->c_str();
      int send_buffer_size = cstrm->size();
      assert(send_buffer_size >= 0);

      // compute the sizes
      std::vector<int> recv_sizes(mpi_size, -1);
      // Compute the sizes
      int error = MPI_Gather(&send_buffer_size,  // Send buffer
                             1,                  // send count
                             MPI_INT,            // send type
                             &(recv_sizes[0]),  // recvbuffer
                             1,                  // recvcount
                             MPI_INT,           // recvtype
                             mpi_rank,          // root rank
                             MPI_COMM_WORLD);  
      assert(error == MPI_SUCCESS);
      for(size_t i = 0; i < recv_sizes.size(); ++i)
        assert(recv_sizes[i] >= 0);
  

      // Construct offsets
      std::vector<int> recv_offsets(recv_sizes);
      int sum = 0, tmp = 0;
      for(size_t i = 0; i < recv_offsets.size(); ++i) {
        tmp = recv_offsets[i]; recv_offsets[i] = sum; sum += tmp; 
      }

      // if necessary realloac recv_buffer
      std::vector<char> recv_buffer(sum);
     
      // recv all the maps 
      error = MPI_Gatherv(send_buffer,         // send buffer
                          send_buffer_size,    // how much to send
                          MPI_BYTE,            // send type
                          &(recv_buffer[0]),   // recv buffer
                          &(recv_sizes[0]),    // amount to recv
                                               // for each cpuess
                          &(recv_offsets[0]),  // where to place data
                          MPI_BYTE,
                          mpi_rank,            // root rank
                          MPI_COMM_WORLD);
      assert(error == MPI_SUCCESS);
      // Update the local map
      namespace bio = boost::iostreams;
      typedef bio::stream<bio::array_source> icharstream;
      icharstream strm(&(recv_buffer[0]), recv_buffer.size());
      graphlab::iarchive iarc(strm);
      for(size_t i = 0; i < results.size(); ++i) {
        iarc >> results[i];
      }  
    } // end of gether



    template<typename T>
    void send(const T& elem, const size_t id, const int tag = 0) {
      // Get the mpi rank and size
      assert(id < size());
      // Serialize the local map
      graphlab::charstream cstrm(128);
      graphlab::oarchive oarc(cstrm);
      oarc << elem;
      cstrm.flush();
      char* send_buffer = cstrm->c_str();
      int send_buffer_size = cstrm->size();
      assert(send_buffer_size >= 0);

      int dest(id);
      // send the size
      int error = MPI_Send(&send_buffer_size,  // Send buffer
                           1,                  // send count
                           MPI_INT,            // send type
                           dest,               // destination
                           tag,                  // tag
                           MPI_COMM_WORLD);
      assert(error == MPI_SUCCESS);
  
      // send the actual content
      error = MPI_Send(send_buffer,         // send buffer
                       send_buffer_size,    // how much to send
                       MPI_BYTE,            // send type
                       dest,
                       tag,
                       MPI_COMM_WORLD);
      assert(error == MPI_SUCCESS);
    } // end of send



    template<typename T>
    void recv(T& elem, const size_t id, const int tag = 0) {
      // Get the mpi rank and size
      assert(id < size());

      int recv_buffer_size(-1);
      int dest(id);
      MPI_Status status;
      // recv the size
      int error = MPI_Recv(&recv_buffer_size,  
                           1,                  
                           MPI_INT,            
                           dest,               
                           tag,                
                           MPI_COMM_WORLD,
                           &status);
      assert(error == MPI_SUCCESS);
      assert(recv_buffer_size > 0);

      std::vector<char> recv_buffer(recv_buffer_size);
      // recv the actual content
      error = MPI_Recv(&(recv_buffer[0]),
                       recv_buffer_size, 
                       MPI_BYTE,         
                       dest,
                       tag,
                       MPI_COMM_WORLD, 
                       &status);
      assert(error == MPI_SUCCESS);
      // deserialize
      // Update the local map
      namespace bio = boost::iostreams;
      typedef bio::stream<bio::array_source> icharstream;
      icharstream strm(&(recv_buffer[0]), recv_buffer.size());
      graphlab::iarchive iarc(strm);
      iarc >> elem;
    }



    uint32_t get_local_ip() {
      uint32_t ip(-1);
      // code adapted from
      struct ifaddrs * ifAddrStruct = NULL;
      getifaddrs(&ifAddrStruct);
      struct ifaddrs * firstifaddr = ifAddrStruct;
      ASSERT_NE(ifAddrStruct, NULL);
      while (ifAddrStruct != NULL) {
        if (ifAddrStruct->ifa_addr != NULL && 
            ifAddrStruct->ifa_addr->sa_family == AF_INET) {
          char* tmpAddrPtr = NULL;
          // check it is IP4 and not lo0.
          tmpAddrPtr = (char*)&((struct sockaddr_in *)ifAddrStruct->ifa_addr)->sin_addr;
          ASSERT_NE(tmpAddrPtr, NULL);
          if (tmpAddrPtr[0] != 127) {
            memcpy(&ip, tmpAddrPtr, 4);
            break;
          }
          //break;
        }
        ifAddrStruct=ifAddrStruct->ifa_next;
      }
      freeifaddrs(firstifaddr);
      return ip;
    }



    void get_master_ranks(std::set<size_t>& master_ranks) {
      uint32_t local_ip = get_local_ip();
      std::vector<uint32_t> all_ips;
      all_gather(local_ip, all_ips);
      std::set<uint32_t> visited_ips;
      master_ranks.clear();
      for(size_t i = 0; i < all_ips.size(); ++i) {
        if(visited_ips.count(all_ips[i]) == 0) {
          visited_ips.insert(all_ips[i]);
          master_ranks.insert(i);
        }
      }
    }





  }; // end of namespace mpi tools
}; //end of graphlab namespace
#include <graphlab/macros_undef.hpp>
#endif
