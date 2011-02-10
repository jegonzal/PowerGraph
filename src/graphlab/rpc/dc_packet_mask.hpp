#ifndef DC_PACKET_MASK_HPP
#define DC_PACKET_MASK_HPP
namespace graphlab {
  // ---------  Packet header types --------------
 
  const unsigned char STANDARD_CALL = 1;
  
  /** If FAST_CALL is set, the RPC call will be 
  issued as early as possible, taking precedence over
  other non-fast-calls. Fast Calls also avoid additional
  buffer allocation. The actual function call,
  should be short, and its arguments should not require 
  much memory. */
  const unsigned char FAST_CALL = 2;

  /** If WAIT_FOR_REPLY is set, the function call's
  return will be passed back to the caller */
  const unsigned char WAIT_FOR_REPLY = 4;
  
  /** The BARRIER packet must have 0 data length
  and must be the only flag set. Upon receipt 
  of a barrier packet, the receiver set a barrier flag.
  All future packets sent from the same host will
  not be processed until all previous packets have
  completed evaluation. This is conceptually
  equivalent to a memory barrier. */
  const unsigned char BARRIER = 8; 
  
  /**
    If control packet flag is set, this packet 
    does not increment any counters.
  */
  const unsigned char CONTROL_PACKET = 16; 
}
#endif
