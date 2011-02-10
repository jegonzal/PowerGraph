
#include <iostream>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_buffered_stream_receive.hpp>

//#define DC_RECEIVE_DEBUG
namespace graphlab {
namespace dc_impl {

/**
  Called by the controller when there is data coming
  from the source
*/
void dc_buffered_stream_receive::incoming_data(procid_t src,
                    const char* buf, 
                    size_t len) {
  bufferlock.lock();
  buffer.write(buf, len);
  recvcond.signal();
  bufferlock.unlock();
}
  
/** called by the controller when a function
call is completed */
void dc_buffered_stream_receive::function_call_completed(unsigned char packettype) {
  size_t pending = pending_calls.dec();
  if ((packettype & CONTROL_PACKET) == 0) callsreceived.inc();
  if (barrier && pending == 0) {
    bufferlock.lock();
    barrier = false;
    recvcond.signal();
    bufferlock.unlock();
  }
}
void dc_buffered_stream_receive::process_buffer() {
  // if barrier is set. we should not process anything
  if (barrier) return;
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

    buffer.skip(sizeof(packet_hdr));

    if ((hdr.packet_type_mask & CONTROL_PACKET) == 0) {
      bytesreceived += hdr.len;
    }

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
    else if (hdr.packet_type_mask & FAST_CALL) {
      // if it is a fast call, dispatch the function immediately
      #ifdef DC_RECEIVE_DEBUG
      logstream(LOG_INFO) << "Is fast call" << std::endl;
      #endif
      boost::iostreams::stream<circular_char_buffer_source> strm(buffer,hdr.len);
      dc->exec_function_call(hdr.src, hdr.packet_type_mask, strm);
      if ((hdr.packet_type_mask & CONTROL_PACKET) == 0) callsreceived.inc();
    }
    else if (hdr.packet_type_mask & STANDARD_CALL) {
      #ifdef DC_RECEIVE_DEBUG
      logstream(LOG_INFO) << "Is deferred call" << std::endl;
      #endif
      // not a fast call. so read out the buffer
      char* tmpbuf = new char[hdr.len];
      buffer.read(tmpbuf, hdr.len);
      pending_calls.inc();
      dc->deferred_function_call(hdr.src,hdr.packet_type_mask, tmpbuf, hdr.len);
    }
  }
}

void dc_buffered_stream_receive::receive_loop() {
  bufferlock.lock();
  while (!done) {
    process_buffer();
    recvcond.wait(bufferlock);
  }
  bufferlock.unlock();
}

size_t dc_buffered_stream_receive::bytes_received() {
  return bytesreceived;
}
size_t dc_buffered_stream_receive::calls_received() {
  return callsreceived.value;
}

void dc_buffered_stream_receive::shutdown() {
  done = true;
  recvcond.signal();
  thr.join();  
}

} // namespace dc_impl
} // namespace graphlab
