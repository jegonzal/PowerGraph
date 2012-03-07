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
    return writebuffer_totallen.value >= buffer_length_trigger;
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

    bool send_decision = false;
    bool signal_decision = false;
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
      writebuffer_totallen.inc(len);
      send_decision = adaptive_send_decision();
      signal_decision = (insertloc == 1);
      // decrement the reference count
      __sync_fetch_and_sub(&(buffer[curid].ref_count), 1);
      break;
    }
    // wake it up from cond sleep
    // first insertion into buffer
    if (signal_decision || send_decision) {
      send_active_lock.lock(); 
        cond.signal();
        send_active_lock.unlock();
      
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
    while (1) {
      if (writebuffer_totallen.value > 0) {
        flush_impl();
      } else {
        // sleep for 1 ms or up till we get wait_count_bytes
        send_active_lock.lock();
        while(1) {
          if (!done) {
            if(writebuffer_totallen.value == 0) {
              cond.wait(send_active_lock);
            }
            else if (writebuffer_totallen.value < buffer_length_trigger) {
              send_active_lock.unlock();
              my_sleep_ms(1);
              send_active_lock.lock();
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
        send_active_lock.unlock();
      }
      if (done) {
        break;
      }
    }
  }

  void dc_buffered_stream_send2::shutdown() {
    send_active_lock.lock();
    done = true;
    cond.signal();
    send_active_lock.unlock();
    thr.join();
  }
  
  void dc_buffered_stream_send2::flush() {
    flush_impl();
  }

  void dc_buffered_stream_send2::flush_impl() {
    // if the writebuffer is empty, just return
    if (writebuffer_totallen.value == 0) return;
    send_lock.lock();
    // swap the buffer
    size_t curid = bufid;
    bufid = !bufid;
    // decrement the reference count
    __sync_fetch_and_sub(&(buffer[curid].ref_count), 1);
    // wait till the reference count is negative
    while(buffer[curid].ref_count >= 0) usleep(1);
    
    // ok. Now we have exclusive access to this buffer
    // take a reference for convenience
    size_t sendlen = buffer[curid].numbytes;
    if (sendlen > 0) {
      size_t numel = std::min((size_t)(buffer[curid].numel.value), buffer[curid].buf.size());
      std::vector<iovec> &sendbuffer = buffer[curid].buf;
      
      writebuffer_totallen.dec(sendlen);    
      sendlength += sendlen;
      block_header_type blockheader = sendlen;
      
      // fill the first msg block
      sendbuffer[0].iov_base = reinterpret_cast<void*>(&blockheader);
      sendbuffer[0].iov_len = sizeof(block_header_type);
      //remember what I just sent so that I can free it later. send_many
      // may modify the vector.
      std::vector<iovec> prevsend(sendbuffer.begin(), 
                                  sendbuffer.begin() + numel);
      comm->send_many(target, sendbuffer, numel);
      // reset the buffer;
      buffer[curid].numbytes = 0;
      buffer[curid].numel = 1;

      if (numel == sendbuffer.size()) {
         sendbuffer.resize(2 * numel);
//         std::cout << "r to " << sendbuffer.size() << std::endl;
      }

      __sync_fetch_and_add(&(buffer[curid].ref_count), 1);
      // now clear what I just sent. start from '1' to avoid the header
      for (size_t i = 1; i < prevsend.size(); ++i) {
        free(prevsend[i].iov_base);
      }
      
      wakeuptimes++;
      if (wakeuptimes & 4) {
        sendlength /= wakeuptimes;
        buffer_length_trigger = (buffer_length_trigger + sendlength) / 2;
        buffer_length_trigger = std::min(buffer_length_trigger, max_buffer_length);
        buffer_length_trigger += (buffer_length_trigger == 0);
        sendlength = 0; wakeuptimes = 0;
      }
    }
    else {
      // reset the buffer;
      buffer[curid].numbytes = 0;
      buffer[curid].numel = 1;
      __sync_fetch_and_add(&(buffer[curid].ref_count), 1);
    }

    send_lock.unlock();
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


