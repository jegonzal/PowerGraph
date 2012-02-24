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

    size_t prevwbufsize = writebuffer.len;
    writebuffer.write(reinterpret_cast<char*>(&hdr), sizeof(packet_hdr));
    writebuffer.write(data, len);

    if (prevwbufsize == 0 && 
        writebuffer.len >= wait_count_bytes &&
        sendlock.try_lock()) {
      comm->send(target, writebuffer.str, writebuffer.len);
      sendlock.unlock();
      writebuffer.clear();
    }
    else if (prevwbufsize == 0 ||
        writebuffer.len >= wait_count_bytes) {
      cond.signal();
    }
    lock.unlock();
  }


  void dc_buffered_stream_send2::send_loop() {
    size_t last_sent = 0;
    graphlab::timer timer;
    timer.start();
    unsigned long long prevtime = rdtsc();
    //const double nano2second = 1000*1000*1000;
    //const double second_wait = nanosecond_wait / nano2second;
    unsigned long long second_wait = estimate_ticks_per_second() / 1000;

    lock.lock();
    while (1) {
      if (writebuffer.len > 0) {
        sendbuffer.swap(writebuffer);
        sendlock.lock();
        lock.unlock();

        wait_count_bytes = (wait_count_bytes + sendbuffer.len)/2;
        wait_count_bytes += (wait_count_bytes == 0);
        comm->send(target, sendbuffer.str, sendbuffer.len);
        sendlock.unlock();
        // shrink if we are not using much buffer
        if (sendbuffer.len < sendbuffer.buffer_size / 2 
            && sendbuffer.buffer_size > 10240) {
          sendbuffer.clear(sendbuffer.buffer_size / 2);
        }
        else {
          sendbuffer.clear();
        }
        lock.lock();
      } else {
        prevtime = rdtsc();
        while(writebuffer.len < wait_count_bytes &&
              prevtime + second_wait > rdtsc() && 
              !done) {
          if(writebuffer.len == 0) cond.wait(lock);
          else cond.timedwait_ns(lock, nanosecond_wait);
        //  std::cout << prevtime << " " << second_wait << " " << nexttime << " " << writebuffer.len << "\n";
        }

      }
      if (done) {
        break;
      }
    }
    lock.unlock();
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


