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
#include <graphlab/rpc/dc_stream_send.hpp>

namespace graphlab {
namespace dc_impl {
void dc_stream_send::send_data(procid_t target, 
                unsigned char packet_type_mask,
                std::istream &istrm,
                size_t len) {

  if (len != size_t(-1)) {
    if ((packet_type_mask & CONTROL_PACKET) == 0) {
      bytessent.inc(len);
    }
    // build the packet header
    packet_hdr hdr;
    memset(&hdr, 0, sizeof(packet_hdr));
    hdr.len = len;
    hdr.src = dc->procid(); 
    hdr.packet_type_mask = packet_type_mask;
    lock.lock();
    comm->send(target, 
                reinterpret_cast<char*>(&(hdr)),
                sizeof(packet_hdr));
    char cbuffer[10240];
    size_t l = 0;
    while(istrm.good()) {
      l = istrm.readsome(cbuffer, 10240);
      if (l == 0) break;  // end of buffer if readsome returns 0
      comm->send(target, cbuffer, l);
    }
    comm->flush(target);
    lock.unlock();
  }
  else {
    // annoying. we have to compute the length of the stream
    // allocate a 128byte block first
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

void dc_stream_send::send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len) {
  if ((packet_type_mask & CONTROL_PACKET) == 0) {
    if (packet_type_mask & (FAST_CALL | STANDARD_CALL)) {
      dc->inc_calls_sent(target);
    }
    bytessent.inc(len);
  }
  packet_hdr hdr;
  memset(&hdr, 0, sizeof(packet_hdr));
  hdr.len = len;
  hdr.src = dc->procid(); 
  hdr.sequentialization_key = dc->get_sequentialization_key();
  hdr.packet_type_mask = packet_type_mask;
  lock.lock();
 /* comm->send(target, 
              reinterpret_cast<char*>(&(hdr)),
              sizeof(packet_hdr));
  comm->send(target, 
             data,
             len); */
  comm->send2(target, reinterpret_cast<char*>(&hdr),sizeof(packet_hdr),
                     data,len);
  comm->flush(target);
  lock.unlock();
}

void dc_stream_send::shutdown() { }

} // namespace dc_impl
} // namespace graphlab

