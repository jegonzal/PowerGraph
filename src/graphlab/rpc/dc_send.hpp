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
#include <sys/types.h>
#include <sys/socket.h>

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
The class should accumulate the data in an iovec structure
and relinquish it on get_outgoing_data()
*/
class dc_send{
 public:
  dc_send() { }
  
  virtual ~dc_send() { }

  /** Called to send data to the target. The caller transfers control of
  the pointer. The caller MUST ensure that the data be prefixed
  with sizeof(packet_hdr) extra bytes at the start for placement of the
  packet header. This function must be reentrant. */
  virtual void send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len) = 0;

  /** Sends the data but without transferring control of the pointer.
   The function will make a copy of the data before sending it.
   Unlike send_data, no padding is necessary. This function must be reentrant. */
  virtual void copy_and_send_data(procid_t target,
                 unsigned char packet_type_mask,
                 char* data, size_t len) = 0;

  /**
    Bytes sent must be incremented BEFORE the data is transmitted.
    Packets marked CONTROL_PACKET should not be counted
  */
  virtual size_t bytes_sent() = 0;
  
  virtual void flush() = 0;
  
  virtual size_t set_option(std::string opt, size_t val) {
    return 0;
  }

  /**
   * Returns true if there is data, false otherwise. This function
   * must be reentrant, but it is guaranteed that only one thread will
   * call this function at anytime. 
   * numel may be less than outdata.size()
   */
  virtual bool get_outgoing_data(std::vector<iovec>& outdata, size_t &numel) = 0;

};
  

} // namespace dc_impl
} // namespace graphlab
#endif

