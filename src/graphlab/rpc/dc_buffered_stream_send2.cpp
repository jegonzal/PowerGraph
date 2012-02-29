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

  bool dc_buffered_stream_send2::adaptive_send_decision() {
    /** basically for each call, you have a decision.
     *  1: send it immediately, (either by waking up the sender, 
     *                           or sending it in-line)
     *  2: buffer it for sending.
     * 
     * I would like to send immediately if the callrate is low
     * (as compared to the thread wake up time.)
     * I would like to buffer it otherwise. 
     */
    return writebuffer.len >= buffer_length_trigger;
  }

  void dc_buffered_stream_send2::send_data(procid_t target, 
                                          unsigned char packet_type_mask,
                                          char* data, size_t len) {
    if ((packet_type_mask & CONTROL_PACKET) == 0) {
      if (packet_type_mask & (STANDARD_CALL)) {
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
    bool send_decision = adaptive_send_decision();
    /*if (0 && sendlock.try_lock()) {
      // try to immediately send if we have exceeded the threshold 
      // already nd we can acquire the lock
      sendbuffer.swap(writebuffer);
      lock.unlock();

      // fill in the chunk header with the length of the chunk
      *reinterpret_cast<block_header_type*>(sendbuffer.str) = (block_header_type)(sendbuffer.len - sizeof(block_header_type));
      comm->send(target, sendbuffer.str, sendbuffer.len);
      // shrink if we are not using much buffer
      if (sendbuffer.len < sendbuffer.buffer_size / 2 
          && sendbuffer.buffer_size > 10240) {
        sendbuffer.clear(sendbuffer.buffer_size / 2);
      }
      else {
        sendbuffer.clear();
      }
      char bufpad[sizeof(block_header_type)];
      sendbuffer.write(bufpad, sizeof(block_header_type));
      sendlock.unlock();
    }*/
    if (prevwbufsize == sizeof(block_header_type) || send_decision) {
      sendlock.lock();
      flush_flag = send_decision;
      cond.signal();
      sendlock.unlock();
      lock.unlock();
    }
    else {
      lock.unlock();
    }
  }


  void dc_buffered_stream_send2::send_loop() {
    graphlab::timer timer;
    timer.start();
    //const double nano2second = 1000*1000*1000;
    //const double second_wait = nanosecond_wait / nano2second;
    int wakeuptimes = 0;
    size_t sendlength = 0;
    lock.lock();
    while (1) {
      if (writebuffer.len > sizeof(block_header_type)) {
        sendbuffer.swap(writebuffer);
        lock.unlock();
        // fill in the chunk header with the length of the chunk
        *reinterpret_cast<block_header_type*>(sendbuffer.str) = (block_header_type)(sendbuffer.len - sizeof(block_header_type));
        comm->send(target, sendbuffer.str, sendbuffer.len);
        wakeuptimes++;
        sendlength += sendbuffer.len;
        if (wakeuptimes & 4) {
          sendlength /= wakeuptimes;
          buffer_length_trigger = (buffer_length_trigger + sendlength) / 2;
          buffer_length_trigger = std::min(buffer_length_trigger, max_buffer_length);
          buffer_length_trigger += (buffer_length_trigger == 0);
          sendlength = 0; wakeuptimes = 0;
        }
        // shrink if we are not using much buffer
        if (sendbuffer.len < sendbuffer.buffer_size / 2 
            && sendbuffer.buffer_size > 10240) {
          sendbuffer.clear(sendbuffer.buffer_size / 2);
        }
        else {
          sendbuffer.clear();
        }
        char bufpad[sizeof(block_header_type)];
        sendbuffer.write(bufpad, sizeof(block_header_type));
    
        lock.lock();
      } else {
        // sleep for 1 ms or up till we get wait_count_bytes
        lock.unlock();
        while(!flush_flag &&
              !done) {
          if (return_signal) {
            return_signal = false;
            flush_return_cond.signal();
          }
          if(writebuffer.len == sizeof(block_header_type)) {
            sendlock.lock();
            cond.wait(sendlock);
            sendlock.unlock();
          }
          else if (writebuffer.len < buffer_length_trigger) {
            my_sleep_ms(1);
            break;
          }
        }
        lock.lock();
        flush_flag = false;
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
  
  void dc_buffered_stream_send2::flush() {
    lock.lock();
    flush_flag = true;
    return_signal = true;
    cond.signal();
    while (return_signal) flush_return_cond.wait(lock);
    lock.unlock();
  }
  
  size_t dc_buffered_stream_send2::set_option(std::string opt, 
                                             size_t val) {
    size_t prevval = 0;
    if (opt == "nanosecond_wait") {
      prevval = nanosecond_wait;
      nanosecond_wait = val;
    }
    else if (opt == "max_buffer_length") {
      prevval = max_buffer_length;
      max_buffer_length = val;
    }
    return prevval;
  }
  
} // namespace dc_impl
} // namespace graphlab


