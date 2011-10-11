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


#ifndef DC_BUFFERED_STREAM_SEND_MULTIQUEUE_HPP
#define DC_BUFFERED_STREAM_SEND_MULTIQUEUE_HPP
#include <iostream>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/rpc/circular_char_buffer.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/non_blocking_queue.hpp>
#include <graphlab/logger/logger.hpp>
namespace graphlab {
class distributed_control;

namespace dc_impl {

struct multiqueue_entry{
  char* c;
  size_t len;
};


/**
   \ingroup rpc
Sender for the dc class.
  The job of the sender is to take as input data blocks of
  pieces which should be sent to a single destination socket.
  This can be thought of as a sending end of a multiplexor.
  This class performs buffered transmissions.
  That is, sends are relegated to an internal circular buffer, which is then
  passed to the communication classes on another thread.

  This implements a buffered sender and can be enabled by passing "buffered_send=yes"
  in the distributed control initstring.
*/

class dc_buffered_stream_send_multiqueue: public dc_send{
 public:
  dc_buffered_stream_send_multiqueue(distributed_control* dc, 
                                     dc_comm_base *comm, 
                                     size_t nummachines,
                                     size_t num_polling_threads): dc(dc), 
                                    comm(comm), done(false) { 
    sendqueues.resize(nummachines);

    // one threads for 2 sockets maximum
    
    for (size_t i = 0;i < num_polling_threads; ++i) {
      size_t sock_range_low = (i * nummachines) / num_polling_threads;
      size_t sock_range_high = ((i + 1) * nummachines) / num_polling_threads;
      pool.launch(boost::bind(&dc_buffered_stream_send_multiqueue::send_loop, 
                              this, sock_range_low, sock_range_high));
    }
  }
  
  ~dc_buffered_stream_send_multiqueue() {
    sendqueues.clear();
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

  void send_loop(size_t sockrangelow, size_t sockrangehigh);
  
  void write_combining_send(size_t target, 
                            non_blocking_queue<multiqueue_entry> &entry);
  
  void shutdown();
  
  inline size_t bytes_sent() {
    return bytessent.value;
  }

 private:
  struct buffer_spec{
    buffer_spec():buf(65536),len(0) { }
    circular_char_buffer buf;
    size_t len;
    mutex lock;
  };
  /// pointer to the owner
  distributed_control* dc;
  dc_comm_base *comm;
  atomic<size_t> bytessent;
  bool done;
  thread_group pool;

  //all queues in queues[i] go to machine i
  std::vector<non_blocking_queue<multiqueue_entry> > sendqueues;
};




} // namespace dc_impl
} // namespace graphlab
#endif

