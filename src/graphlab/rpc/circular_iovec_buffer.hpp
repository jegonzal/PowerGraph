/*  
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


#ifndef GRAPHLAB_RPC_CIRCULAR_IOVEC_BUFFER_HPP
#define GRAPHLAB_RPC_CIRCULAR_IOVEC_BUFFER_HPP
#include <vector>
#include <sys/socket.h>

namespace graphlab{
namespace dc_impl {
 
/**
 * \ingroup rpc
 * \internal
 * A circular buffer which maintains a parallel sequence of iovecs.
 * One sequence is basic iovecs
 * The other sequence is used for storing the original unomidifed pointers
 * This is minimally checked. length must be a power of 2
 */
struct circular_iovec_buffer {
  inline circular_iovec_buffer(size_t len = 4096) {
    v.resize(4096);
    parallel_v.resize(4096);
    head = 0;
    tail = 0;
    numel = 0;
  }
  
  inline bool empty() const {
    return numel == 0;
  }
  
  size_t size() const {
    return numel;
  }
  
  /**
   * Writes an entry into the buffer, resizing the buffer if necessary.
   * This buffer will take over all iovec pointers and free them when done
   */
  inline void write(const iovec &entry) {
    if (numel == v.size()) {
      std::vector<struct iovec> newv(v.size() * 2);
      std::vector<struct iovec> new_parallel_v(v.size() * 2);
      size_t newi = 0;
      // copy to the new vector
      if (head < tail) {
        // head is before tail. just copy to tail and I am done
        while(head < tail) {
          newv[newi] = v[head];
          new_parallel_v[newi] = parallel_v[head];
          ++newi; ++head;
        }
      }
      else {
        // otherwise there is a loop around
        while(head < numel) {
           newv[newi] = v[head];
           new_parallel_v[newi] = parallel_v[head];
          ++newi; ++head;
        }
        head = 0;
        while(head < tail) {
          newv[newi] = v[head];
          new_parallel_v[newi] = parallel_v[head];
          ++newi; ++head;
        }
      }
      v.swap(newv);
      parallel_v.swap(new_parallel_v);
      head = 0;
      tail = newi;
    }
    
    v[tail] = entry;
    parallel_v[tail] = entry;
    tail = (tail + 1) & (v.size() - 1); ++numel;
  }

  
  /**
   * Erases a single iovec from the head and free the pointer
   */
  inline void erase_from_head_and_free() {
    free(v[head].iov_base);
    head = (head + 1) & (v.size() - 1);
    --numel;
  }

  /**
   * Fills a msghdr for unsent data.
   */
  void fill_msghdr(struct msghdr& data) {
    data.msg_iov = &(parallel_v[head]);
    if (head < tail) {
      data.msg_iovlen = tail - head;
    }
    else {
      data.msg_iovlen = v.size() - head;
    }
    data.msg_iovlen = std::min<size_t>(IOV_MAX, data.msg_iovlen);
  }
  
  /**
   * Advances the head as if some amount of data was sent.
   */
  void sent(size_t len) {
    while(len > 0) {
      size_t curv_sent_len = std::min(len, parallel_v[head].iov_len);
      parallel_v[head].iov_len -= curv_sent_len;
      parallel_v[head].iov_base = (char*)(parallel_v[head].iov_base) + curv_sent_len;
      len -= curv_sent_len;
      if (parallel_v[head].iov_len == 0) {
        erase_from_head_and_free();
      }
    }
  }

  std::vector<struct iovec> v;
  std::vector<struct iovec> parallel_v;
  size_t head;
  size_t tail;
  size_t numel;
};
  
}
}

#endif