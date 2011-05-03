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

#ifndef DC_PACKET_MASK_HPP
#define DC_PACKET_MASK_HPP
namespace graphlab {
  // ---------  Packet header types --------------

  /** 
   * \ingroup rpc_internal 
   * 
   * Used for regular calls which go into a thread pool
   * for evaluation
   */
  const unsigned char STANDARD_CALL = 1;
  
  /** 
   * \ingroup rpc_internal 
   * 
   * If FAST_CALL is set, the RPC call will be 
  issued as early as possible, taking precedence over
  other non-fast-calls. Fast Calls also avoid additional
  buffer allocation. The actual function call,
  should be short, and its arguments should not require 
  much memory. */
  const unsigned char FAST_CALL = 2;

  /** \ingroup rpc_internal
   * 
   * If WAIT_FOR_REPLY is set, the function call's
  return will be passed back to the caller */
  const unsigned char WAIT_FOR_REPLY = 4;
  
  /** \ingroup rpc_internal
   * 
   * The BARRIER packet must have 0 data length
  and must be the only flag set. Upon receipt 
  of a barrier packet, the receiver set a barrier flag.
  All future packets sent from the same host will
  not be processed until all previous packets have
  completed evaluation. This is conceptually
  equivalent to a memory barrier. */
  const unsigned char BARRIER = 8; 
  
  /**
    \ingroup rpc_internal
   * 
    If control packet flag is set, this packet 
    does not increment any counters.
  */
  const unsigned char CONTROL_PACKET = 16; 
  
  /**
   * \ingroup rpc_internal
   * 
   * Used to identify that this packet was 
   * a reply to a previous request.
   */
  const unsigned char REPLY_PACKET = 32;
}
#endif
