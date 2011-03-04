#ifndef DC_STREAM_RECEIVE_LZ_HPP
#define DC_STREAM_RECEIVE_LZ_HPP
#include <boost/type_traits/is_base_of.hpp>
#include <graphlab/extern/quicklz/quicklz.hpp>
#include <graphlab/rpc/circular_char_buffer.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_receive.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/logger.hpp>

namespace graphlab {
class distributed_control;



namespace dc_impl {

/**
  Receiver processor for the dc class.
  The job of the receiver is to take as input a byte stream
  (as received from the socket) and cut it up into meaningful chunks.
  This can be thought of as a receiving end of a multiplexor.
*/
class dc_stream_receive_lz: public dc_receive{
 public:
  
  dc_stream_receive_lz(distributed_control* dc): 
                  buffer(10240),
                  barrier(false), dc(dc),
                  bytesreceived(0){ 
    state_decompress = (qlz_state_decompress *)malloc(sizeof(qlz_state_decompress));
    zbuffer = (char*)malloc(128*1024 + 400);
    decompbuffer = (char*)malloc(128*1024);
  }
  
  ~dc_stream_receive_lz() {
    free(state_decompress);
    free(zbuffer);
    free(decompbuffer);
  }

  /**
   Called by the controller when there is data coming
   from the source
  */
  void incoming_data(procid_t src, 
                     const char* buf, 
                     size_t len);
   
  /** called by the controller when a function
  call is completed */
  void function_call_completed(unsigned char packettype) ;
 private:
  /// the mutex protecting the buffer and the barrier 
  mutex bufferlock;
  
  /** the incoming data stream. This is protected
  by the bufferlock */
  circular_char_buffer buffer;

  /** number of rpc calls from this other processor
     which are in the deferred execution queue */
  atomic<size_t> pending_calls;
  
  /** whether a barrier has been issued. 
      this is protected by the bufferlock */
  bool barrier;
  
  /// pointer to the owner
  distributed_control* dc;

  size_t bytesreceived;
  atomic<size_t> compressed_bytesreceived;
  
  /**
    Reads the incoming buffer and processes, dispatching
    calls when enough bytes are received
  */
  void process_buffer(bool outsidelocked) ;

  size_t bytes_received();
  
  void shutdown();

  inline bool direct_access_support() {
    return false;
  }
  
  char* get_buffer(size_t& retbuflength);
  

  char* advance_buffer(char* c, size_t wrotelength, 
                              size_t& retbuflength);
                              
	qlz_state_decompress *state_decompress;
  char* zbuffer;
  size_t zbufferoffset;

  char* decompbuffer;
  
};


} // namespace dc_impl
} // namespace graphlab
#endif
