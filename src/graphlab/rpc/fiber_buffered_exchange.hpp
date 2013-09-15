/**
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef GRAPHLAB_FIBER_BUFFERED_EXCHANGE_HPP
#define GRAPHLAB_FIBER_BUFFERED_EXCHANGE_HPP

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/fiber_control.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/util/mpi_tools.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * \ingroup rpc
   * \internal
   */
  template<typename T>
  class fiber_buffered_exchange {
  public:
    typedef std::vector<T> buffer_type;

    struct buffer_record {
      procid_t proc;
      buffer_type buffer;
      buffer_record() : proc(-1)  { }
    }; // end of buffer record
    typedef std::vector<buffer_record> recv_buffer_type;
  private:



    /** The rpc interface for this class */
    mutable dc_dist_object<fiber_buffered_exchange> rpc;

    std::vector<std::vector< buffer_record> > recv_buffers;


    struct send_record {
      oarchive* oarc;
      size_t numinserts;
    };

    std::vector<std::vector<send_record> > send_buffers;
    const size_t max_buffer_size;


    // typedef boost::function<void (const T& tref)> handler_type;
    // handler_type recv_handler;

  public:
    fiber_buffered_exchange(distributed_control& dc,
                      const size_t max_buffer_size = DEFAULT_BUFFERED_EXCHANGE_SIZE) :
      rpc(dc, this),
      max_buffer_size(max_buffer_size) {
       send_buffers.resize(fiber_control::get_instance().num_workers());
       recv_buffers.resize(fiber_control::get_instance().num_workers());
       for (size_t i = 0;i < send_buffers.size(); ++i) {
         send_buffers[i].resize(dc.numprocs());
         for (size_t j = 0;j < send_buffers[i].size(); ++j) {
           send_buffers[i][j].oarc = NULL;
           send_buffers[i][j].numinserts = 0;
         }
       }
       rpc.barrier();
      }


    ~fiber_buffered_exchange() {
      // clear the send buffers
      for (size_t i = 0;i < send_buffers.size(); ++i) {
        for (size_t j = 0;j < send_buffers[i].size(); ++j) {
          if (send_buffers[i][j].oarc) rpc.split_call_cancel(send_buffers[i][j].oarc);
        }
      }
    }
    // fiber_buffered_exchange(distributed_control& dc, handler_type recv_handler,
    //                   size_t buffer_size = 1000) :
    // rpc(dc, this), send_buffers(dc.numprocs()), send_locks(dc.numprocs()),
    // max_buffer_size(buffer_size), recv_handler(recv_handler) { rpc.barrier(); }


    void send(const procid_t proc, const T& value) {
      size_t wid = fiber_control::get_worker_id();
      if (send_buffers[wid][proc].oarc == NULL) {
        send_buffers[wid][proc].oarc = rpc.split_call_begin(&fiber_buffered_exchange::rpc_recv);
        // write a header
        (*send_buffers[wid][proc].oarc) << rpc.procid();
        send_buffers[wid][proc].numinserts = 0;
      }

      (*(send_buffers[wid][proc].oarc)) << value;
      ++send_buffers[wid][proc].numinserts;


      if(send_buffers[wid][proc].oarc->off >= max_buffer_size) {
        flush_buffer(wid, proc);
      }
    } // end of send

    void flush_buffer(size_t wid, procid_t proc) {
      if(send_buffers[wid][proc].oarc) {
        // write the length at the end of the buffere are returning
        (*send_buffers[wid][proc].oarc) << (size_t)(send_buffers[wid][proc].numinserts);
        rpc.split_call_end(proc, send_buffers[wid][proc].oarc);
//         logstream(LOG_DEBUG) << rpc.procid() << ": Sending exchange of length " 
//                              << send_buffers[wid][proc].oarc->off << " to " 
//                              << proc << std::endl;
        send_buffers[wid][proc].oarc = NULL;
        send_buffers[wid][proc].numinserts = 0;
      }
    }

    void partial_flush() {
      for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
        flush_buffer(fiber_control::get_worker_id(), proc);
      }
    }

    void flush() {
      for(size_t i = 0; i < send_buffers.size(); ++i) {
        for (size_t j = 0;j < send_buffers[i].size(); ++j) {
          flush_buffer(i,j);
        }
      }
      rpc.dc().flush();
      rpc.full_barrier();
    } // end of flush


    bool recv(std::vector<buffer_record>& ret_buffer,
              const bool self_buffer = true) {
      fiber_control::fast_yield();
      ret_buffer.clear();
      bool success = false;
      if (self_buffer) {
        // get from my own buffer
        size_t wid = fiber_control::get_worker_id();
        if(!recv_buffers[wid].empty()) {
          success = true;
          std::swap(ret_buffer, recv_buffers[wid]);
        }
      } else {
        for (size_t i = 0;i < recv_buffers.size(); ++i) {
          if(!recv_buffers[i].empty()) {
            success = true;
            std::swap(ret_buffer, recv_buffers[i]);
            break;
          }
        }
      }
      return success;
    } // end of recv



    /**
     * Returns the number of elements to recv
     */
    size_t size() const {
      size_t count = 0;
      for (size_t i = 0;i < recv_buffers.size(); ++i) {
        count += recv_buffers[i].size();
      }
      return count;
    } // end of size

    bool empty() const { 
      for (size_t i = 0;i < recv_buffers.size(); ++i) {
        if (recv_buffers[i].size() > 0) return false;
      }
      return true;
    }

    void clear() { }

    void barrier() { rpc.barrier(); }
  private:
    void rpc_recv(size_t len, wild_pointer w) {
      buffer_type tmp;
      iarchive iarc(reinterpret_cast<const char*>(w.ptr), len);
      // first desrialize the source process
      procid_t src_proc; iarc >> src_proc;
//       logstream(LOG_DEBUG) << rpc.procid() << ": Receiving exchange of length "
//                            << len << " from " << src_proc << std::endl;
      // create an iarchive which just points to the last size_t bytes
      // to get the number of elements
      iarchive numel_iarc(reinterpret_cast<const char*>(w.ptr) + len - sizeof(size_t),
                          sizeof(size_t));
      size_t numel; numel_iarc >> numel;
      //std::cout << "Receiving: " << numel << "\n";
      tmp.resize(numel);
      for (size_t i = 0;i < numel; ++i) {
        iarc >> tmp[i];
      }

      size_t wid = fiber_control::get_worker_id();
      recv_buffers[wid].push_back(buffer_record());
      buffer_record& rec = recv_buffers[wid].back();
      rec.proc = src_proc;
      rec.buffer.swap(tmp);
    } // end of rpc rcv



  }; // end of buffered exchange


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>

#endif


