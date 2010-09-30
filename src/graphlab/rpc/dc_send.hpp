#ifndef DC_SEND_HPP
#define DC_SEND_HPP
#include <iostream>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
namespace graphlab {
namespace dc_impl {

/**
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
  call with when there is data to send. The data pointer's
  ownership is transfered here and is this class's responsibility
  to free the pointer when done. */
  virtual void send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len) = 0;

};
  

} // namespace dc_impl
} // namespace graphlab
#endif

