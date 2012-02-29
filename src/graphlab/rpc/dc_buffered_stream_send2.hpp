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


#ifndef DC_BUFFERED_STREAM_SEND2_HPP
#define DC_BUFFERED_STREAM_SEND2_HPP
#include <iostream>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/resizing_array_sink.hpp>
#include <graphlab/logger/logger.hpp>
namespace graphlab {
class distributed_control;

namespace dc_impl {


/**
   \ingroup rpc
Sender for the dc class.
  The job of the sender is to take as input data blocks of
  pieces which should be sent to a single destination socket.
  This can be thought of as a sending end of a multiplexor.
  This class performs buffered transmissions using an blocking 
  queue with one call per queue entry.
  A seperate thread is used to transmit queue entries. Rudimentary
  write combining is used to decrease transmission overhead.
  This is typically the best performing sender.
  
  This can be enabled by passing "buffered_queued_send=yes"
  in the distributed control initstring.
  
  dc_buffered_stream_send22 is similar, but does not perform write combining.
  
*/

class dc_buffered_stream_send2: public dc_send{
 public:
  dc_buffered_stream_send2(distributed_control* dc, 
                                   dc_comm_base *comm, 
                                   procid_t target) : 
                  dc(dc),  comm(comm), target(target), done(false), 
                  flush_flag(false), return_signal(false),
                  rtdsc_per_ms(estimate_ticks_per_second() / 1000) {
    char bufpad[sizeof(block_header_type)];
    writebuffer.write(bufpad, sizeof(block_header_type));
    sendbuffer.write(bufpad, sizeof(block_header_type));
    thr = launch_in_new_thread(boost::bind
                               (&dc_buffered_stream_send2::send_loop, 
                                this));
  }
  
  ~dc_buffered_stream_send2() {
  }
  

  inline bool channel_active(procid_t target) const {
    return comm->channel_active(target);
  }

  /**
   Called by the controller when there is data to send.
   if len is -1, the function has to compute the length by itself,
   or send the data from the stream directly. the strm is not copyable.
  */
  void send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 std::istream &istrm,
                 size_t len = size_t(-1));
                 
  /** Another possible interface the controller can
  call with when there is data to send. The caller has
  responsibility for freeing the pointer when this call returns*/
  void send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len);

  void send_loop();
  
  void flush();
  
  void shutdown();
  
  bool adaptive_send_decision();
  
  inline size_t bytes_sent() {
    return bytessent.value;
  }

  /**
   * Possible Options include 
   * nanosecond_wait: Maximum amount of time remote calls can age 
   *                  in the queue before the queue is flushed. (1000000)
   * wait_count_bytes: Maximum number of bytes in the buffer before
   *                   the queue is flushed. This number is self adjusting
   *                   with an exponential scaling rate and should not need
   *                   to be modified. (initial = 1024000)
   */
  size_t set_option(std::string opt, size_t val);

 private:
  /// pointer to the owner
  distributed_control* dc;
  dc_comm_base *comm;
  procid_t target;
  
  charstream_impl::resizing_array_sink<true> writebuffer;
  charstream_impl::resizing_array_sink<true> sendbuffer;
  
  mutex lock;
  mutex sendlock;
  conditional cond;

  thread thr;
  bool done;

  static atomic<size_t> callcount;
  atomic<size_t> bytessent; 
  
  bool flush_flag;
  bool return_signal;
  conditional flush_return_cond;
  
  static double calls_per_ms;
  
  size_t nanosecond_wait;
  static unsigned long long prevtime;
  unsigned long long rtdsc_per_ms;
  static mutex callcountmutex;
};



} // namespace dc_impl
} // namespace graphlab
#endif // DC_BUFFERED_STREAM_SEND_EXPQUEUE_HPP

