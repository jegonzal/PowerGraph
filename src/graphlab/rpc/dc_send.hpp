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

  virtual void shutdown() = 0;

};
  

} // namespace dc_impl
} // namespace graphlab
#endif

