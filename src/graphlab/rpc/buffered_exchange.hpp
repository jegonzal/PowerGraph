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


#ifndef GRAPHLAB_BUFFERED_EXCHANGE_HPP
#define GRAPHLAB_BUFFERED_EXCHANGE_HPP

#include <graphlab/parallel/pthread_tools.hpp>
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
  class buffered_exchange {
  public:
    typedef std::vector<T> buffer_type;

  private:
    struct buffer_record {
      procid_t proc;
      buffer_type buffer;
      buffer_record() : proc(-1)  { }
    }; // end of buffer record



    /** The rpc interface for this class */
    mutable dc_dist_object<buffered_exchange> rpc;

    std::deque< buffer_record > recv_buffers;
    mutex recv_lock;


    struct send_record {
      oarchive* oarc;
      size_t numinserts;
    };

    std::vector<send_record> send_buffers;
    std::vector< mutex >  send_locks;
    const size_t num_threads;
    const size_t max_buffer_size;


    // typedef boost::function<void (const T& tref)> handler_type;
    // handler_type recv_handler;

  public:
    buffered_exchange(distributed_control& dc,
                      const size_t num_threads = 1,
                      const size_t max_buffer_size = 1024 * 1024 /* 1MB */) :
      rpc(dc, this),
      send_buffers(num_threads *  dc.numprocs()),
      send_locks(num_threads *  dc.numprocs()),
      num_threads(num_threads),
      max_buffer_size(max_buffer_size) {
       //
       for (size_t i = 0;i < send_buffers.size(); ++i) {
         // initialize the split call
         send_buffers[i].oarc = rpc.split_call_begin(&buffered_exchange::rpc_recv);
         send_buffers[i].numinserts = 0;
         // begin by writing the src proc.
         (*(send_buffers[i].oarc)) << rpc.procid();
       }
       rpc.barrier();
      }


    ~buffered_exchange() {
      // clear the send buffers
      for (size_t i = 0;i < send_buffers.size(); ++i) {
        rpc.split_call_cancel(send_buffers[i].oarc);
      }
    }
    // buffered_exchange(distributed_control& dc, handler_type recv_handler,
    //                   size_t buffer_size = 1000) :
    // rpc(dc, this), send_buffers(dc.numprocs()), send_locks(dc.numprocs()),
    // max_buffer_size(buffer_size), recv_handler(recv_handler) { rpc.barrier(); }


    void send(const procid_t proc, const T& value, const size_t thread_id = 0) {
      ASSERT_LT(proc, rpc.numprocs());
      ASSERT_LT(thread_id, num_threads);
      const size_t index = thread_id * rpc.numprocs() + proc;
      ASSERT_LT(index, send_locks.size());
      send_locks[index].lock();

      (*(send_buffers[index].oarc)) << value;
      ++send_buffers[index].numinserts;

      if(send_buffers[index].oarc->off >= max_buffer_size) {
        oarchive* prevarc = swap_buffer(index);
        send_locks[index].unlock();
        // complete the send
        rpc.split_call_end(proc, prevarc);
      } else {
        send_locks[index].unlock();
      }
    } // end of send


    void partial_flush(size_t thread_id) {
      for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
        const size_t index = thread_id * rpc.numprocs() + proc;
        ASSERT_LT(proc, rpc.numprocs());
        if (send_buffers[index].numinserts > 0) {
          send_locks[index].lock();
          oarchive* prevarc = swap_buffer(index);
          send_locks[index].unlock();
          // complete the send
          rpc.split_call_end(proc, prevarc);
        }
      }
    }

    void flush() {
      for(size_t i = 0; i < send_buffers.size(); ++i) {
        const procid_t proc = i % rpc.numprocs();
        ASSERT_LT(proc, rpc.numprocs());
        send_locks[i].lock();
        if (send_buffers[i].numinserts > 0) {
          oarchive* prevarc = swap_buffer(i);
          // complete the send
          rpc.split_call_end(proc, prevarc);
        }
        send_locks[i].unlock();
      }
      rpc.full_barrier();
    } // end of flush


    bool recv(procid_t& ret_proc, buffer_type& ret_buffer,
              const bool try_lock = false) {
      dc_impl::blob read_buffer;
      bool has_lock = false;
      if(try_lock) {
        if (recv_buffers.empty()) return false;
        has_lock = recv_lock.try_lock();
      } else {
        recv_lock.lock();
        has_lock = true;
      }
      bool success = false;
      if(has_lock) {
        if(!recv_buffers.empty()) {
          success = true;
          buffer_record& rec =  recv_buffers.front();
          // read the record
          ret_proc = rec.proc;
          ret_buffer.swap(rec.buffer);
          ASSERT_LT(ret_proc, rpc.numprocs());
          recv_buffers.pop_front();
        }
        recv_lock.unlock();
      }

      return success;
    } // end of recv



    /**
     * Returns the number of elements to recv
     */
    size_t size() const {
      typedef typename std::deque< buffer_record >::const_iterator iterator;
      recv_lock.lock();
      size_t count = 0;
      foreach(const buffer_record& rec, recv_buffers) {
        count += rec.buffer.size();
      }
      recv_lock.unlock();
      return count;
    } // end of size

    bool empty() const { return recv_buffers.empty(); }

    void clear() { }

    void barrier() { rpc.barrier(); }
  private:
    void rpc_recv(size_t len, wild_pointer w) {
      buffer_type tmp;
      iarchive iarc(reinterpret_cast<const char*>(w.ptr), len);
      // first desrialize the source process
      procid_t src_proc; iarc >> src_proc;
      ASSERT_LT(src_proc, rpc.numprocs());
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

      recv_lock.lock();
      recv_buffers.push_back(buffer_record());
      buffer_record& rec = recv_buffers.back();
      rec.proc = src_proc;
      rec.buffer.swap(tmp);
      recv_lock.unlock();
    } // end of rpc rcv


    // create a new buffer for send_buffer[index], returning the old buffer
    oarchive* swap_buffer(size_t index) {
      oarchive* swaparc = rpc.split_call_begin(&buffered_exchange::rpc_recv);
      swaparc->expand_buf(max_buffer_size * 1.2);
      std::swap(send_buffers[index].oarc, swaparc);
      // write the length at the end of the buffere are returning
      (*swaparc) << (size_t)(send_buffers[index].numinserts);

      //std::cout << "Sending : " << (send_buffers[index].numinserts)<< "\n";
      // reset the insertion count
      send_buffers[index].numinserts = 0;
      // write the current procid into the new buffer
      (*(send_buffers[index].oarc)) << rpc.procid();
      return swaparc;
    }


  }; // end of buffered exchange


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>

#endif


