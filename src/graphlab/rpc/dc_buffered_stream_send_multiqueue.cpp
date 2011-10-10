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
#include <graphlab/rpc/dc_buffered_stream_send_multiqueue.hpp>

namespace graphlab {
namespace dc_impl {


void dc_buffered_stream_send_multiqueue::send_data(procid_t target, 
                            unsigned char packet_type_mask,
                            std::istream &istrm,
                            size_t len) {
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

void dc_buffered_stream_send_multiqueue::send_data(procid_t target, 
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
  
  
  // whats the thread id?
  size_t targetqueue = thread::thread_id() % queues[target].size();
  
  buffer_spec& targetbuffer = queues[target][targetqueue];
  targetbuffer.lock.lock();
  targetbuffer.buf.write(reinterpret_cast<char*>(&(hdr)), sizeof(packet_hdr));
  targetbuffer.buf.write(data, len);
  targetbuffer.len += len + sizeof(packet_hdr);
  targetbuffer.lock.unlock();
}


void dc_buffered_stream_send_multiqueue::send_loop(size_t sockrangelow, 
                                        size_t sockrangehigh) {
  while (!done) {
    for (size_t i = sockrangelow; i < sockrangehigh; ++i) {
      for (size_t j = 0;j < queues[i].size(); ++j) {
        buffer_spec& sourcebuffer = queues[i][j];
        if (sourcebuffer.len > 0) {
          sourcebuffer.lock.lock();
          while(1) {
            char* c = NULL;
            std::streamsize len = sourcebuffer.buf.introspective_read(c, 512*1024);
            if (len > 0) {
              comm->send(i, c, len);
            }
            else {
              break;
            }
          }
          sourcebuffer.len = 0;
          sourcebuffer.lock.unlock();
        }
      }
    }
    sched_yield();
  }
}

void dc_buffered_stream_send_multiqueue::shutdown() {
  done = true;
  while (1) {
    try {
      pool.join();
      break;
    }
    catch(const char*c) {
      std::cout << c << std::endl;
    }
  }
}



} // namespace dc_impl
} // namespace graphlab
