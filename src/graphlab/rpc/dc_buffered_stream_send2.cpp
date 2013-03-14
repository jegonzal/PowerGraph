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

  void dc_buffered_stream_send2::send_data(procid_t target,
                                           unsigned char packet_type_mask,
                                           char* data, size_t len) {
    size_t actual_data_len = len - sizeof(packet_hdr) - sizeof(size_t);

    if ((packet_type_mask & CONTROL_PACKET) == 0) {
      if (packet_type_mask & (STANDARD_CALL)) {
        dc->inc_calls_sent(target);
      }
    }
    writebuffer_totallen.inc(actual_data_len);

    // build the packet header
    packet_hdr* hdr = reinterpret_cast<packet_hdr*>(data + sizeof(size_t));
    memset(hdr, 0, sizeof(packet_hdr));

    hdr->len = len - sizeof(packet_hdr);
    hdr->src = dc->procid();
    hdr->sequentialization_key = dc->get_sequentialization_key();
    hdr->packet_type_mask = packet_type_mask;

    sendqueue.enqueue(data);

    size_t sqsize = approx_send_queue_size;
    approx_send_queue_size = sqsize + 1;
    if (sqsize == 256) comm->trigger_send_timeout(target, false);
    else if ((packet_type_mask &
            (CONTROL_PACKET | WAIT_FOR_REPLY | REPLY_PACKET))) {
      comm->trigger_send_timeout(target, true);
    }
  }

  void dc_buffered_stream_send2::flush() {
    comm->trigger_send_timeout(target, true);
  }

  void dc_buffered_stream_send2::copy_and_send_data(procid_t target,
                                          unsigned char packet_type_mask,
                                          char* data, size_t len) {
    char* c = (char*)malloc(sizeof(packet_hdr) + len);
    memcpy(c + sizeof(size_t) + sizeof(packet_hdr), data, len);
    send_data(target, packet_type_mask, c, len + sizeof(size_t) + sizeof(packet_hdr));
  }


  size_t dc_buffered_stream_send2::get_outgoing_data(circular_iovec_buffer& outdata) {
    // fast exit if no buffer
    if (writebuffer_totallen.value == 0) return 0;

    char* sendqueue_head = sendqueue.dequeue_all();
    if (sendqueue_head == NULL) return 0;

    approx_send_queue_size = 0;
    size_t real_send_len = 0;

    // construct the block msg header
    block_header_type* blockheader = new block_header_type;
    // now I don't really know what is the size of it yet.
    // create a block header iovec
    iovec blockheader_iovec;
    blockheader_iovec.iov_base = reinterpret_cast<void*>(blockheader);
    blockheader_iovec.iov_len = sizeof(block_header_type);
    outdata.write(blockheader_iovec);

    while(!sendqueue.end_of_dequeue_list(sendqueue_head)) {
      iovec tosend, tofree;
      tofree.iov_base = sendqueue_head;
      tosend.iov_base = sendqueue_head + sizeof(size_t);
      // I need to read the length
      packet_hdr* hdr = reinterpret_cast<packet_hdr*>(sendqueue_head + sizeof(size_t));
      tosend.iov_len = hdr->len + sizeof(packet_hdr);
      tofree.iov_len = tosend.iov_len + sizeof(size_t);
      outdata.write(tosend, tofree);
      real_send_len += tosend.iov_len;
      // advance to the next list item
      while(__unlikely__(inplace_lf_queue::get_next(sendqueue_head) == NULL)) {
        asm volatile("pause\n": : :"memory");
      }
      sendqueue_head = inplace_lf_queue::get_next(sendqueue_head);
    }
    (*blockheader) = real_send_len;
    writebuffer_totallen.dec(real_send_len);
    return real_send_len;
  }
} // namespace dc_impl
} // namespace graphlab


