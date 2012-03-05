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

  void dc_buffered_stream_send2::copy_and_send_data(procid_t target, 
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

    if (prevwbufsize == sizeof(block_header_type) || send_decision) {
      buffer_empty_lock.lock();
      flush_flag = send_decision;
      cond.signal();
      buffer_empty_lock.unlock();
      lock.unlock();
    }
    else {
      lock.unlock();
    }
  }


  void dc_buffered_stream_send2::send_data(procid_t target,
                                          unsigned char packet_type_mask,
                                          char* data, size_t len) {
    copy_and_send_data(target, packet_type_mask, data, len);
    free(data);
  }
    
  void dc_buffered_stream_send2::send_loop() {
    graphlab::timer timer;
    timer.start();
    //const double nano2second = 1000*1000*1000;
    //const double second_wait = nanosecond_wait / nano2second;
    lock.lock();
    while (1) {
      if (writebuffer.len > sizeof(block_header_type)) {
        flush_locked();
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
            buffer_empty_lock.lock();
            cond.wait(buffer_empty_lock);
            buffer_empty_lock.unlock();
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
    flush_locked();
    lock.unlock();
  }

  void dc_buffered_stream_send2::flush_locked() {
    // if the writebuffer is empty, just return
    // note that "empty" means it must contain at leat the header byes.
    if (writebuffer.len == sizeof(block_header_type)) return;
    sendbuffer.swap(writebuffer);
    send_lock.lock();
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
    send_lock.unlock();
    lock.lock();
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


