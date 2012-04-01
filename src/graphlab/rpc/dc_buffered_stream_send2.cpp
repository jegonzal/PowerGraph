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


#include <iostream>
#include <boost/iostreams/stream.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_buffered_stream_send2.hpp>

namespace graphlab {
namespace dc_impl {

  void dc_buffered_stream_send2::send_data(procid_t target,
                                           unsigned char packet_type_mask,
                                           char* data, size_t len) {
    if ((packet_type_mask & CONTROL_PACKET) == 0) {
      if (packet_type_mask & (STANDARD_CALL)) {
        dc->inc_calls_sent(target);
      }
      bytessent.inc(len - sizeof(packet_hdr));
    }
    
    // build the packet header
    packet_hdr* hdr = reinterpret_cast<packet_hdr*>(data);
    memset(hdr, 0, sizeof(packet_hdr));

    hdr->len = len - sizeof(packet_hdr);
    hdr->src = dc->procid();
    hdr->sequentialization_key = dc->get_sequentialization_key();
    hdr->packet_type_mask = packet_type_mask;
    iovec msg;
    msg.iov_base = data;
    msg.iov_len = len;
    bool trigger = false;
    while(1) {
      size_t curid;
      while(1) {
        curid = bufid;
        int32_t cref = buffer[curid].ref_count;
        if (cref < 0 || 
            !atomic_compare_and_swap(buffer[curid].ref_count, cref, cref + 1)) continue;

        if (curid != bufid) {
          __sync_fetch_and_sub(&(buffer[curid].ref_count), 1);
        }
        else {
          break;
        }
      }
      // ok, we have a reference count into curid, we can write to it
      size_t insertloc = buffer[curid].numel.inc_ret_last();
      // ooops out of buffer room. release the reference count, flush and retry
      if (insertloc >= buffer[curid].buf.size()) {
        __sync_fetch_and_sub(&(buffer[curid].ref_count), 1);
        usleep(1);
        continue;
      }
      buffer[curid].buf[insertloc] = msg;
      buffer[curid].numbytes.inc(len);    
      trigger = ((writebuffer_totallen.inc_ret_last(len)) == 0);
      // decrement the reference count
      __sync_fetch_and_sub(&(buffer[curid].ref_count), 1);
      break;
    }
    
    if (trigger) comm->trigger_send_timeout(target);
  }

  void dc_buffered_stream_send2::flush() {
    while(writebuffer_totallen.value) usleep(1);
  }

  void dc_buffered_stream_send2::copy_and_send_data(procid_t target,
                                          unsigned char packet_type_mask,
                                          char* data, size_t len) {
    char* c = (char*)malloc(sizeof(packet_hdr) + len);
    memcpy(c + sizeof(packet_hdr), data, len);
    send_data(target, packet_type_mask, c, len + sizeof(packet_hdr));
  }


  bool dc_buffered_stream_send2::get_outgoing_data(std::vector<iovec>& outdata,
                                                   size_t& numel) {
    if (writebuffer_totallen.value == 0) return false;
    
    // swap the buffer
    size_t curid = bufid;
    bufid = !bufid;
    // decrement the reference count
    __sync_fetch_and_sub(&(buffer[curid].ref_count), 1);
    // wait till the reference count is negative
    while(buffer[curid].ref_count >= 0) usleep(1);
    
    // ok now we have exclusive access to the buffer
    size_t sendlen = buffer[curid].numbytes;
    if (sendlen > 0) {
      numel = std::min((size_t)(buffer[curid].numel.value), buffer[curid].buf.size());
      bool buffull = (numel == buffer[curid].buf.size());
      std::vector<iovec> &sendbuffer = buffer[curid].buf;
      
      writebuffer_totallen.dec(sendlen);    
      block_header_type* blockheader = new block_header_type;
      (*blockheader) = sendlen;
      
      // fill the first msg block
      sendbuffer[0].iov_base = reinterpret_cast<void*>(blockheader);
      sendbuffer[0].iov_len = sizeof(block_header_type);
      // give the buffer away
      outdata.swap(sendbuffer);
      // reset the buffer;
      buffer[curid].numbytes = 0;
      buffer[curid].numel = 1;

      if (buffull) {
        sendbuffer.resize(2 * numel);
      }
      else {
        sendbuffer.resize(numel);
      }
      __sync_fetch_and_add(&(buffer[curid].ref_count), 1);
      return true;
    }
    else {
      // reset the buffer;
      buffer[curid].numbytes = 0;
      buffer[curid].numel = 1;
      __sync_fetch_and_add(&(buffer[curid].ref_count), 1);
      return false;
    }
  }
} // namespace dc_impl
} // namespace graphlab


