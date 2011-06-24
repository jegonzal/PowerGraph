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

