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


#include <iostream>
#include <boost/iostreams/stream.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_buffered_stream_send2.hpp>
#include <graphlab/util/branch_hints.hpp>
namespace graphlab {
namespace dc_impl {

  void dc_buffered_stream_send2::flush() {
    comm->trigger_send_timeout(target, true);
  }

  void dc_buffered_stream_send2::flush_soon() {
    comm->trigger_send_timeout(target, true);
  }



  inline size_t dc_buffered_stream_send2::bytes_sent() {
    size_t ret = total_bytes_sent;
    lock.lock();
    for (size_t i = 0;i < send_buffers.size(); ++i) {
      ret += send_buffers[i]->get_bytes_sent(target); 
    }
    lock.unlock();
    return ret;
  }

  void dc_buffered_stream_send2::write_to_buffer(char* c, size_t len)  {
    lock.lock();
    additional_flush_buffers.push_back(std::make_pair(c, len));
    lock.unlock();
  }

  void dc_buffered_stream_send2::register_send_buffer(thread_local_buffer* buffer) {
    lock.lock();
    send_buffers.push_back(buffer);
    to_send.resize(send_buffers.size());
    lock.unlock();
  }

  void dc_buffered_stream_send2::unregister_send_buffer(thread_local_buffer* buffer) {
    lock.lock();
    for (size_t i = 0;i < send_buffers.size(); ++i) {
      if (send_buffers[i] == buffer) {
        total_bytes_sent.inc(send_buffers[i]->get_bytes_sent(target));
        send_buffers.erase(send_buffers.begin() + i);
        break;
      }
    }
    to_send.resize(send_buffers.size());
    lock.unlock();
  }

  dc_buffered_stream_send2::~dc_buffered_stream_send2() {
    // unregister all the buffers.
    std::vector<thread_local_buffer*> all_buffers;
    for (size_t i = 0; i < all_buffers.size(); ++i) {
      unregister_send_buffer(all_buffers[i]);
    }
  }

  size_t dc_buffered_stream_send2::get_outgoing_data(circular_iovec_buffer& outdata) {
    lock.lock();

    size_t sendlen = 0;
    for (size_t i = 0;i < send_buffers.size(); ++i) {
      to_send[i] = send_buffers[i]->extract(target);
      for (size_t j = 0;j < to_send[i].size(); ++j) {
        sendlen += to_send[i][j].second;
      }
    }
    for (size_t i = 0;i < additional_flush_buffers.size(); ++i) {
      sendlen += additional_flush_buffers[i].second;
    }
    // nothing to send
    if (sendlen == 0) {
      lock.unlock();
      return 0;
    }

    // ok. I have something to send.
    // construct the block msg header
//     block_header_type* blockheader = new block_header_type;
//     (*blockheader) = sendlen;
//     iovec blockheader_iovec;
//     blockheader_iovec.iov_base = reinterpret_cast<void*>(blockheader);
//     blockheader_iovec.iov_len = sizeof(block_header_type);
//     outdata.write(blockheader_iovec);

    for (size_t i = 0;i < send_buffers.size(); ++i) {
      for (size_t j = 0;j < to_send[i].size(); ++j) {
        if (to_send[i][j].second > 0) {
          iovec sendvec;
          sendvec.iov_base = to_send[i][j].first;
          sendvec.iov_len = to_send[i][j].second;
          outdata.write(sendvec);
        }
      }
    }

    for (size_t i = 0;i < additional_flush_buffers.size(); ++i) {
      iovec sendvec;
      sendvec.iov_base = additional_flush_buffers[i].first;
      sendvec.iov_len = additional_flush_buffers[i].second;
      outdata.write(sendvec);
    }
    additional_flush_buffers.clear();
    lock.unlock();
    return sendlen + sizeof(block_header_type);
  }
} // namespace dc_impl
} // namespace graphlab


