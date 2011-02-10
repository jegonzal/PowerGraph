#ifndef DC_RECEIVE_HPP
#define DC_RECEIVE_HPP
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/parallel/atomic.hpp>
namespace graphlab {
namespace dc_impl {
  
/**
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
  Calls received is incremented only after completion
  of a function handler. Control packets should not be counted.
  */
  virtual size_t calls_received() = 0;
  
  virtual void shutdown() = 0;
};


} // namespace dc_impl
} // namespace graphlab
#endif
