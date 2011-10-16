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
#include <graphlab/rpc/dc_buffered_stream_send_expqueue.hpp>

namespace graphlab {
namespace dc_impl {

void dc_buffered_stream_send_expqueue::send_data(procid_t target_, 
                unsigned char packet_type_mask,
                std::istream &istrm,
                size_t len) {
  ASSERT_EQ(target, target_);
  if (len != size_t(-1)) {
    char cbuffer[len];
    while(len > 0 && istrm.good()) {
      size_t l = istrm.readsome(cbuffer, len);
      len -= l;
    }
    send_data(target, packet_type_mask, cbuffer, len);
  }
  else {
    // annoying. we have to compute the length of the stream
    // allocate a 128byte block first.
    // \todo: This can be optimized. Though, I don't think this
    //        code path is even used.
    size_t len = 0;
    size_t cursize = 128;
    char* data = (char*)malloc(128);
    // while the stream is good. read stuff
    // len is the current length of the contents
    // cursize is the max length of the array
    // when we run out of space in the array, we double the size.
    while (istrm.good()) {
      len += istrm.readsome(data+len, cursize-len);
      if (cursize - len == 0) {
        cursize *= 2;
        data = (char*)realloc(data, cursize);
      }
    }
    send_data(target, packet_type_mask, data, len);
    free(data);
  }
}

void dc_buffered_stream_send_expqueue::send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len) {
  if ((packet_type_mask & CONTROL_PACKET) == 0) {
    if (packet_type_mask & (FAST_CALL | STANDARD_CALL)) {
      dc->inc_calls_sent(target);
    }
    bytessent.inc(len);
  }

  // build the packet header
  packet_hdr hdr;
  memset(&hdr, 0, sizeof(packet_hdr));
  
  hdr.len = len;
  hdr.src = dc->procid(); 
  hdr.sequentialization_key = dc->get_sequentialization_key();
  hdr.packet_type_mask = packet_type_mask;
  
  std::streamsize numbytes_needed = sizeof(packet_hdr) + len;
  expqueue_entry eentry;
  eentry.len = numbytes_needed;
  eentry.c = (char*)malloc(numbytes_needed);
  memcpy(eentry.c, &hdr, sizeof(packet_hdr));
  memcpy(eentry.c + sizeof(packet_hdr), data, len);
  sendqueue.enqueue_conditional_signal(eentry, 32);
}


void dc_buffered_stream_send_expqueue::write_combining_send(std::deque<expqueue_entry> &stuff) {
  // fill new entries until I fill up the length
  char sendcombining[combine_upper_threshold];
  size_t length = 0;
  // now send the stuff. stick them into the combining buffer.
  while(!stuff.empty()){  
    expqueue_entry ent = stuff.front();
    stuff.pop_front();
    // if length with the new entry exceeds my combining buffer length
    // dump the combining buffer first and clear the buffer
    if (length + ent.len > combine_upper_threshold) {
      comm->send(target, sendcombining, length);
      length = 0;
    }

    // if there is room in the combining buffer, write into the combining buffer
    if (length + ent.len <= combine_upper_threshold) {
      memcpy(sendcombining + length, ent.c , ent.len);
      length += ent.len;
      free(ent.c);
    }
    else {
      // length should be 0 here
      // too long for the combining buffer
      // send it seperately
      comm->send(target, ent.c, ent.len); 
      free(ent.c);
    }
  }
  
  if (length > 0) {
    comm->send(target, sendcombining, length);
  }
}

void dc_buffered_stream_send_expqueue::send_loop() {
  float t = lowres_time_seconds(); 
  while (1) {
    bool ret = sendqueue.timed_wait_for_data(100000, 32);
    if (ret) {
      std::deque<expqueue_entry> stuff_to_send;
      sendqueue.swap(stuff_to_send);
      if (lowres_time_seconds() - t > 10) {
        t = lowres_time_seconds();
        std::cout << dc->procid() << "->" << target << " send buffer = " << stuff_to_send.size() << std::endl;
      }
      write_combining_send(stuff_to_send);
    }
    else {
      break;
    }
  }
}

void dc_buffered_stream_send_expqueue::shutdown() {
  sendqueue.stop_blocking();
  thr.join();
}

} // namespace dc_impl
} // namespace graphlab


