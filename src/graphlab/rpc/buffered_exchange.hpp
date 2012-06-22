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
      size_t numel;
      dc_impl::blob buffer;
      buffer_record() : proc(-1), numel(0) { }
    }; // end of buffer record



    /** The rpc interface for this class */
    mutable dc_dist_object<buffered_exchange> rpc;

    std::deque< buffer_record > recv_buffers;
    mutex recv_lock;


    struct send_record {
      send_record():buffer(128),numinserts(0){}
      // need a fake copy constructor here
      // just so I can make a vector of these
      send_record(const send_record& ):buffer(128), numinserts(0) { }
      charstream buffer;
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
                      const size_t max_buffer_size = 1000000) : 
      rpc(dc, this), 
      send_buffers(num_threads *  dc.numprocs()), 
      send_locks(num_threads *  dc.numprocs()),
      num_threads(num_threads),
      max_buffer_size(max_buffer_size) { 
       rpc.barrier(); 
      }



    ~buffered_exchange() { 
     foreach(buffer_record& buf, recv_buffers) {
        buf.buffer.free();
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

      oarchive oarc(send_buffers[index].buffer);
      ++send_buffers[index].numinserts;
      oarc << value;
      if(send_buffers[index].buffer->size() > max_buffer_size) {
        send_buffers[index].buffer.flush();
        graphlab::dc_impl::blob b(send_buffers[index].buffer->str, 
                                  send_buffers[index].buffer->len);
        if(proc == rpc.procid()) {
          rpc_recv(proc, send_buffers[index].numinserts, b);
          // here the blob is transimtted directly
        } else {
          rpc.remote_call(proc, &buffered_exchange::rpc_recv,
                          rpc.procid(), send_buffers[index].numinserts, b);
          // here I need to free the blob
          b.free();
        }
        send_buffers[index].buffer->relinquish();
        send_buffers[index].numinserts = 0;
      }
      send_locks[index].unlock();
    } // end of send


    void flush() {
      for(size_t i = 0; i < send_buffers.size(); ++i) {
        const procid_t proc = i % rpc.numprocs();
        ASSERT_LT(proc, rpc.numprocs());
        send_locks[i].lock();
        send_buffers[i].buffer.flush();
        graphlab::dc_impl::blob b(send_buffers[i].buffer->str, 
                                  send_buffers[i].buffer->len);
        if(proc == rpc.procid()) {
          rpc_recv(proc, send_buffers[i].numinserts, b);
        } else {
          rpc.remote_call(proc, &buffered_exchange::rpc_recv,
                          rpc.procid(),send_buffers[i].numinserts, b);
          b.free();
        }
        send_buffers[i].buffer->relinquish();
        send_buffers[i].numinserts = 0;
        send_locks[i].unlock();
      }
      rpc.full_barrier();
    } // end of flush


    bool recv(procid_t& ret_proc, buffer_type& ret_buffer, 
              const bool try_lock = false) {
      dc_impl::blob read_buffer;
      size_t numel = 0;
      bool has_lock = false;
      if(try_lock) {
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
          read_buffer = rec.buffer;
          numel = rec.numel;
          ASSERT_LT(ret_proc, rpc.numprocs());
          recv_buffers.pop_front();
        }
        recv_lock.unlock();
      }

     ret_buffer.clear();
      if (success) {
        // deserialize tmp_buffer
        boost::iostreams::stream<boost::iostreams::array_source> strm(read_buffer.c, read_buffer.len);
        iarchive iarc(strm);
        ret_buffer.resize(numel);
        for (size_t i = 0;i < numel; ++i) {
          iarc >> ret_buffer[i];
        }
        read_buffer.free();
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
        count += rec.buffer.len;
      }
      recv_lock.unlock();
      return count;
    } // end of size

    bool empty() const { return recv_buffers.empty(); }

    void clear() {
      foreach(buffer_record& buf, recv_buffers) {
        buf.buffer.free();
      }
    }

  private:
    void rpc_recv(procid_t src_proc, size_t numel, dc_impl::blob& buffer) {
      recv_lock.lock();
      recv_buffers.push_back(buffer_record());
      buffer_record& rec = recv_buffers.back();
      rec.proc = src_proc;
      rec.numel = numel;
      rec.buffer = buffer;
      recv_lock.unlock();
    } // end of rpc rcv

  }; // end of buffered exchange


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>

#endif


