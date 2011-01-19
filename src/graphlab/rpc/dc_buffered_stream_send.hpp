#ifndef DC_BUFFERED_STREAM_SEND_HPP
#define DC_BUFFERED_STREAM_SEND_HPP
#include <iostream>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/rpc/circular_char_buffer.hpp>
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
  Essentially each "send call" here must end up being paired with a 
  corresponding "receive call" on the receiver end.
*/

class dc_buffered_stream_send: public dc_send{
 public:
  dc_buffered_stream_send(distributed_control* dc, dc_comm_base *comm): dc(dc), comm(comm), done(false){ 
    thr = launch_in_new_thread(boost::bind(&dc_buffered_stream_send::send_loop, 
                                      this));
  }
  
  ~dc_buffered_stream_send() {
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

  void send_loop();
  
  void shutdown();
  
 private:
  /// pointer to the owner
  distributed_control* dc;
  dc_comm_base *comm;

  mutex sendbuflock;
  circular_char_buffer sendbuf;
  conditional sendcond;

  thread thr;
  bool done;
};



} // namespace dc_impl
} // namespace graphlab
#endif

