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

#ifndef DC_RECEIVE_HPP
#define DC_RECEIVE_HPP
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/parallel/atomic.hpp>
namespace graphlab {
namespace dc_impl {
  
/**
\ingroup rpc_internal
Base class of the data receiving class.
This class forms the receiving side of a "multiplexer"
Data entering from a single socket will be passed to this
function through the incoming_data function call.

This class must understand the packet header and issue the right
calls in the owning dc. 
*/
class dc_receive {
 public:
  dc_receive() { };
  virtual ~dc_receive() { };

  /**
   Called by the controller when there is data coming
   from the source
  */
  virtual void incoming_data(procid_t src, 
                     const char* buf, 
                     size_t len) = 0;
   
  /**
    If direct access support return true, 
    get buffer and commit buffer should be defined
  */
  virtual bool direct_access_support() = 0;
  
  /**
    gets a buffer. The buffer length is returned in retbuflength
    This will be used for receiving data.
    If get_buffer() or advance_buffer() is called,
    incoming_data will never be called.
  */
  virtual char* get_buffer(size_t& retbuflength) = 0;
  
  /**
    Commits a buffer obtained using get_buffer.
    c will be the result of a previous call to get_buffer() or advance_buffer()
    This function should commit a range of bytes starting of c,
    up to 'wrotelength' bytes. A new empty buffer should be returned
    and the size is returned in retbuflength
  */
  virtual char* advance_buffer(char* c, size_t wrotelength, 
                              size_t& retbuflength) = 0;
  
  /** called by the controller when a function
  call is completed */
  virtual void function_call_completed(unsigned char packettype) = 0;

  /**
  Bytes received must be updated before handing off 
  to the handlers. Bytes received do not count headers.
  If packet type is marked as CONTROL_PACKET, the packet is not counted
  */
  virtual size_t bytes_received() = 0;
  
  
  /**
   * Last call sent to any instance of dc_receive.
   * If the sender multithreads, the sending thread must shut down.
   */
  virtual void shutdown() = 0;
};


} // namespace dc_impl
} // namespace graphlab
#endif
