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

  void dc_buffered_stream_send2::send_data(procid_t target_, 
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

  void dc_buffered_stream_send2::send_data(procid_t target, 
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

    lock.lock();
    writebuffer.write(reinterpret_cast<char*>(&hdr), sizeof(packet_hdr));
    writebuffer.write(data, len);
    if (writebuffer.size() >= wait_count_bytes) cond.signal();
    lock.unlock();
  }


  void dc_buffered_stream_send2::send_loop() {
    size_t last_sent = 0;
    graphlab::timer timer;
    timer.start();
    unsigned long long last_time = rdtsc();
    //const double nano2second = 1000*1000*1000;
    //const double second_wait = nanosecond_wait / nano2second;
    unsigned long long second_wait = estimate_ticks_per_second() / 1000;
    std::cout << second_wait << std::endl;
    while (1) {
      lock.lock();
      cond.timedwait_ns(lock, nanosecond_wait);
      if (!done && writebuffer.size() > 0) {
        sendbuffer.swap(writebuffer);
        lock.unlock();
        last_sent += writebuffer.size();
        comm->send(target, sendbuffer.str, sendbuffer.len);
        // shrink if we are not using much buffer
        if (sendbuffer.len < sendbuffer.buffer_size / 2 
            && sendbuffer.buffer_size > 10240) {
          sendbuffer.clear(sendbuffer.buffer_size / 2);
        }
        else {
          sendbuffer.clear();
        }
      } else {
        lock.unlock();
      }
      if (done) {
        break;
      }
      unsigned long long curtime = rdtsc();
      if(curtime - last_time >= second_wait) {
        wait_count_bytes = (wait_count_bytes + last_sent)/2;
        last_sent = 0;
        last_time = curtime;
      }

    }
  }

  void dc_buffered_stream_send2::shutdown() {
    lock.lock();
    done = true;
    cond.signal();
    lock.unlock();
    thr.join();
  }
  
  size_t dc_buffered_stream_send2::set_option(std::string opt, 
                                             size_t val) {
    size_t prevval = 0;
    if (opt == "nanosecond_wait") {
      prevval = nanosecond_wait;
      nanosecond_wait = val;
    }
    else if (opt == "wait_count_bytes") {
      prevval = wait_count_bytes;
      wait_count_bytes = val;
    }
    return prevval;
  }
  
} // namespace dc_impl
} // namespace graphlab


