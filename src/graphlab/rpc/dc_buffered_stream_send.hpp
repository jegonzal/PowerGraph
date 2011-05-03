/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DC_BUFFERED_STREAM_SEND_HPP
#define DC_BUFFERED_STREAM_SEND_HPP
#include <iostream>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/util/safe_circular_char_buffer.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/blocking_queue.hpp>
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
  This class performs buffered transmissions.
  That is, sends are relegated to an internal circular buffer, which is then
  passed to the communication classes on another thread.

  This implements a buffered sender and can be enabled by passing "buffered_send=yes"
  in the distributed control initstring.
*/

class dc_buffered_stream_send: public dc_send{
 public:
  dc_buffered_stream_send(distributed_control* dc, dc_comm_base *comm, procid_t target): dc(dc), 
                                    comm(comm), target(target), done(false) { 
    thr = launch_in_new_thread(boost::bind(&dc_buffered_stream_send::send_loop, 
                                      this));
  }
  
  ~dc_buffered_stream_send() {
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
  
  void shutdown();
  
  inline size_t bytes_sent() {
    return bytessent.value;
  }

 private:
  /// pointer to the owner
  distributed_control* dc;
  dc_comm_base *comm;
  procid_t target;
  safe_circular_char_buffer sendbuf;

  thread thr;
  bool done;
  atomic<size_t> bytessent;
  
  void send_till_empty();
};



} // namespace dc_impl
} // namespace graphlab
#endif

