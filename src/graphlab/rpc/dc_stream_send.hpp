#ifndef DC_STREAM_SEND_HPP
#define DC_STREAM_SEND_HPP
#include <iostream>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/blocking_queue.hpp>
#include <graphlab/logger/logger.hpp>
namespace graphlab {
class distributed_control;

namespace dc_impl {

/**
  Sender for the dc class.
  The job of the sender is to take as input data blocks of
  pieces which should be sent to a single destination socket.
  This can be thought of as a sending end of a multiplexor.
*/

class dc_stream_send: public dc_send{
 public:
  dc_stream_send(distributed_control* dc, dc_comm_base *comm): dc(dc), comm(comm){ 
 
  }
  
  ~dc_stream_send() {
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
  call with when there is data to send. The data pointer's
  ownership is transfered here and is this class's responsibility
  to free the pointer when done. */
  void send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len);

  void shutdown();
  inline size_t bytes_sent() {
    return bytessent.value;
  }

 private:
  mutex lock;
  
  /// pointer to the owner
  distributed_control* dc;
  dc_comm_base *comm;
  atomic<size_t> bytessent, callssent;

};



} // namespace dc_impl
} // namespace graphlab
#endif

