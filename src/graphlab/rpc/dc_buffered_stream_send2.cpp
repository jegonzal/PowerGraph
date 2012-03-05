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

  inline bool dc_buffered_stream_send2::adaptive_send_decision() {
    /** basically for each call, you have a decision.
     *  1: send it immediately, (either by waking up the sender, 
     *                           or sending it in-line)
     *  2: buffer it for sending.
     * 
     * I would like to send immediately if the callrate is low
     * (as compared to the thread wake up time.)
     * I would like to buffer it otherwise. 
     */
    return writebuffer_totallen >= buffer_length_trigger;
  }

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
    lock.lock();
    writebuffer.push_back(msg);
    writebuffer_totallen += len;
    bool send_decision = adaptive_send_decision();
    // wake it up from cond sleep
    bool signal_decision = writebuffer.size() == 2;
    lock.unlock();
    // first insertion into buffer
    if (signal_decision || send_decision) {
      buffer_empty_lock.lock();
      cond.signal();
      buffer_empty_lock.unlock();
    }
  }


  void dc_buffered_stream_send2::copy_and_send_data(procid_t target,
                                          unsigned char packet_type_mask,
                                          char* data, size_t len) {
    char* c = (char*)malloc(sizeof(packet_hdr) + len);
    memcpy(c + sizeof(packet_hdr), data, len);
    send_data(target, packet_type_mask, c, len + sizeof(packet_hdr));
  }


  void dc_buffered_stream_send2::send_loop() {
    graphlab::timer timer;
    timer.start();
    //const double nano2second = 1000*1000*1000;
    //const double second_wait = nanosecond_wait / nano2second;
    lock.lock();
    while (1) {
      if (writebuffer_totallen > 0) {
        flush_locked();
      } else {
        // sleep for 1 ms or up till we get wait_count_bytes
        lock.unlock();
        buffer_empty_lock.lock();
        while(1) {
          if (!done) {
            if(writebuffer_totallen == 0) {
              cond.wait(buffer_empty_lock);
            }
            else if (writebuffer_totallen < buffer_length_trigger) {
              my_sleep_ms(1);
              break;
            }
            else {
              break;
            }
          }
          else {
            break;
          }
        }
        buffer_empty_lock.unlock();
        lock.lock();
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
    if (writebuffer_totallen == 0) return;
    sendbuffer.swap(writebuffer);
    size_t sendlen = writebuffer_totallen;
    writebuffer_totallen = 0;
    send_lock.lock();
    lock.unlock();
    block_header_type blockheader = sendlen;
    // fill the first msg block
    sendbuffer[0].iov_base = reinterpret_cast<void*>(&blockheader);
    sendbuffer[0].iov_len = sizeof(block_header_type);
    //remember what I just sent so that I can free it later. send_many
    // may modify the vector.
    std::vector<iovec> prevsend = sendbuffer;
    comm->send_many(target, sendbuffer);
    wakeuptimes++;
    sendlength += sendlen;
    if (wakeuptimes & 4) {
      sendlength /= wakeuptimes;
      buffer_length_trigger = (buffer_length_trigger + sendlength) / 2;
      buffer_length_trigger = std::min(buffer_length_trigger, max_buffer_length);
      buffer_length_trigger += (buffer_length_trigger == 0);
      sendlength = 0; wakeuptimes = 0;
    }
    sendbuffer.resize(1);
    send_lock.unlock();

    // now clear what I just sent. start from '1' to avoid the header
    for (size_t i = 1; i < prevsend.size(); ++i) {
      free(prevsend[i].iov_base);
    }
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


