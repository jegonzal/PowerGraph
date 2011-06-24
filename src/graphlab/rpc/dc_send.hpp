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


#ifndef DC_SEND_HPP
#define DC_SEND_HPP
#include <iostream>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
namespace graphlab {
namespace dc_impl {

/**
\ingroup rpc_internal
Base class of the data sending class.
This class forms the sending side of a "multiplexer"
send_data() will be called with a packet mask as well as a
character stream containing the contents of the packet.
This class must then transmit the data out of the associated 
comms.
*/
class dc_send{
 public:
  dc_send() { }
  virtual ~dc_send() { }
  /**
   Called by the controller when there is data to send.
   if len is -1, the function has to compute the length by itself,
   or send the data from the stream directly. the strm is not copyable.
  */
  virtual void send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 std::istream &istrm,
                 size_t len = size_t(-1)) = 0;
  /** Another possible interface the controller can
  call with when there is data to send. The caller has
  responsibility for freeing the pointer when this call returns*/
  virtual void send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len) = 0;

  /**
    Bytes sent must be incremented BEFORE the data is transmitted.
    Packets marked CONTROL_PACKET should not be counted
  */
  virtual size_t bytes_sent() = 0;
  

  /** returns true if the channel to the target
  machine is truly open. The dc_comm_base specification allows
  for lazy channels which are not created until it is used.
  For such implementations, this function should return true
  if the channel has been created, and false otherwise. Non-lazy
  implementations should return true all the time.
  The invariant to ensure is that this function must return true
  for a target machine ID if a packet has been sent from this machine
  to the target before this call.
  */
  virtual bool channel_active(procid_t target) const = 0;

  /**
   * Last call sent to any instance of dc_send.
   * If the sender multithreads, the sending thread must shut down.
   */
  virtual void shutdown() = 0;

};
  

} // namespace dc_impl
} // namespace graphlab
#endif

