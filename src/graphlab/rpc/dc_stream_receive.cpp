
#include <iostream>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_stream_receive.hpp>

//#define DC_RECEIVE_DEBUG
namespace graphlab {
namespace dc_impl {

/**
  Called by the controller when there is data coming
  from the source
*/
void dc_stream_receive::incoming_data(procid_t src, 
                    const char* buf, 
                    size_t len) {
  bufferlock.lock();
  buffer.write(buf, len);
  bufferlock.unlock();
  process_buffer();
}
  
/** called by the controller when a function
call is completed */
void dc_stream_receive::function_call_completed() {
  size_t pending = pending_calls.dec();
  if (barrier && pending == 0) {
    bufferlock.lock();
    barrier = false;
    bufferlock.unlock();
    process_buffer();
  }
}
void dc_stream_receive::process_buffer() {
  // if barrier is set. we should not process anything
  if (barrier) return;
  if (bufferlock.try_lock()) {
    // only makes sense to process if we at least have
    // a header
    while (size_t(buffer.size()) >= sizeof(packet_hdr)) {
      // read the header
      packet_hdr hdr;
      buffer.peek((char*)(&hdr), sizeof(hdr));
      #ifdef DC_RECEIVE_DEBUG
      logstream(LOG_INFO) << "peeked packet header. Has length " 
                          << hdr.len << std::endl;         
      #endif
      //do we have enough to extract a single packet
      // if not, quit now!
      if (size_t(buffer.size()) < sizeof(packet_hdr) + hdr.len) break;


      if (hdr.packet_type_mask & BARRIER) {
        #ifdef DC_RECEIVE_DEBUG
        logstream(LOG_INFO) << "Comm barrier" << std::endl;
        #endif
        ASSERT_EQ(hdr.len, 0); // barrier packets cannot contain data
        // barrier only makes sense if we have incomplete calls
        barrier = pending_calls.value > 0;
        // ok. we do have incomplete calls. quit processing.
        if (barrier) break;
      }
      if (hdr.packet_type_mask & FAST_CALL) {
        // if it is a fast call, dispatch the function immediately
        buffer.skip(sizeof(packet_hdr));

        #ifdef DC_RECEIVE_DEBUG
        logstream(LOG_INFO) << "Is fast call" << std::endl;
        #endif
        boost::iostreams::stream<circular_char_buffer_source> strm(buffer,hdr.len);
        dc->exec_function_call(hdr.src, strm);
      }
      else {
        #ifdef DC_RECEIVE_DEBUG
        logstream(LOG_INFO) << "Is deferred call" << std::endl;
        #endif
        buffer.skip(sizeof(packet_hdr));
        // not a fast call. so read out the buffer
        char* tmpbuf = new char[hdr.len];
        buffer.read(tmpbuf, hdr.len);
        pending_calls.inc();
        dc->deferred_function_call(hdr.src, tmpbuf, hdr.len);
      }
    }
    bufferlock.unlock();
  }
}

void dc_stream_receive::shutdown() { }

} // namespace dc_impl
} // namespace graphlab
