
#include <iostream>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_stream_receive_lz.hpp>

//#define DC_RECEIVE_DEBUG

namespace graphlab {
namespace dc_impl {



/**
  Called by the controller when there is data coming
  from the source
*/
void dc_stream_receive_lz::incoming_data(procid_t src, 
                    const char* buf, 
                    size_t len) {
  const size_t MAXBUFLEN = 128*1024 + 400;
  do{
    // copy into zbuffer if the offset is small.
    // since the zbuffer block is at least 9 bytes, if the 
    // current buffer is smaller than 9, it must be the beginning of the
    // next block. make sure we get 9 bytes first
    if (zbufferoffset < 9) {
      size_t copylen = std::min(9 - zbufferoffset, len);
      memcpy(zbuffer + zbufferoffset, buf, copylen);
      zbufferoffset += copylen;
      buf += copylen;
      len -= copylen;
      // if we still don't have enough data, we are done. quit
      if (zbufferoffset < 9) {
        ASSERT_EQ(len, 0);
        break;
      }
    }
    // we must have at least 9 bytes here
    //inspect the header
    size_t bodylen = qlz_size_compressed(zbuffer);
    ASSERT_LT(bodylen, MAXBUFLEN);
    // make sure the buffer contains the entire body
    if (zbufferoffset < bodylen) {
      size_t copylen = std::min((bodylen - zbufferoffset), len);
      memcpy(zbuffer + zbufferoffset, buf, copylen);
      zbufferoffset += copylen;
      buf += copylen;
      len -= copylen;
    }    
    // do we have enough data?
    if (bodylen <= zbufferoffset) {
      //yes!
      // decompress!
      size_t bytesout = qlz_decompress(zbuffer, decompbuffer, state_decompress);
      buffer.write(decompbuffer, bytesout);
      zbufferoffset -= bodylen;
      // by design, we only copied barely enough into the buffer.
      ASSERT_EQ(zbufferoffset, 0);
    }
    else {
      // insufficient data to decompress
      // wait for more. len MUST be zero at this point
      ASSERT_EQ(len, 0);
    }
  } while(len > 0);
  process_buffer(false);

  //std::cout << dc->procid() << ": From " << src << ": " << zstrm.total_in << " --> " << zstrm.total_out << std::endl; 
}
  
/** called by the controller when a function
call is completed */
void dc_stream_receive_lz::function_call_completed(unsigned char packettype) {
  size_t pending = pending_calls.dec();
  if (barrier && pending == 0) {
    bufferlock.lock();
    barrier = false;
    bufferlock.unlock();
    process_buffer(false);
  }
}
void dc_stream_receive_lz::process_buffer(bool outsidelocked) {
  // if barrier is set. we should not process anything
  if (barrier) return;
  if (outsidelocked || bufferlock.try_lock()) {
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
      else if ((hdr.packet_type_mask & FAST_CALL) || (hdr.packet_type_mask & REPLY_PACKET)) {
        // if it is a fast call, dispatch the function immediately. replies are also fast tracked
        #ifdef DC_RECEIVE_DEBUG
        if (hdr.packet_type_mask & REPLY_PACKET) logstream(LOG_INFO) << "Is reply" << std::endl;
        else logstream(LOG_INFO) << "Is fast call" << std::endl;
        #endif
        boost::iostreams::stream<circular_char_buffer_source> strm(buffer,hdr.len);
        dc->exec_function_call(hdr.src,hdr, strm);
      }
      else if (hdr.packet_type_mask & STANDARD_CALL) {
        #ifdef DC_RECEIVE_DEBUG
        logstream(LOG_INFO) << "Is deferred call" << std::endl;
        #endif
        // not a fast call. so read out the buffer
        char* tmpbuf = new char[hdr.len];
        buffer.read(tmpbuf, hdr.len);
        pending_calls.inc();
        dc->deferred_function_call(hdr.src,hdr, tmpbuf, hdr.len);
      }
    }
    if (!outsidelocked) bufferlock.unlock();
  }
}

char* dc_stream_receive_lz::get_buffer(size_t& retbuflength) {
  char* ret;
  bufferlock.lock();
  // get a write section
  retbuflength = buffer.introspective_write(ret);
  assert(retbuflength > 0);
  bufferlock.unlock();
  return ret;
}


char* dc_stream_receive_lz::advance_buffer(char* c, size_t wrotelength, 
                            size_t& retbuflength) {
  char* ret;
  bufferlock.lock();
  buffer.advance_write(wrotelength);

  process_buffer(true);
  
  retbuflength = buffer.introspective_write(ret);
  // if the writeable section is too small, sqeeze the buffer 
  // and try again
  if (retbuflength < 1024) {
    // realign the buffer if it is cheap to do so
    if (buffer.align_requires_alloc() == false) {
      buffer.align();
      retbuflength = buffer.introspective_write(ret);
    }
    // try again
    // if this is still too small
    if (retbuflength < 1024) {
      //reserve more capacity
      buffer.reserve(2 * buffer.reserved_size());
      retbuflength = buffer.introspective_write(ret);
    }
    
    // if still too small
    if (retbuflength < 1024) {
      //reserve more capacity
      buffer.reserve(2 * buffer.reserved_size());
      retbuflength = buffer.introspective_write(ret);    
    }
    // by design of the circular buffer, this is guaranteed to 
    // free up enough space
  }
  bufferlock.unlock();

  return ret;
}


size_t dc_stream_receive_lz::bytes_received() {
  return bytesreceived;
}
  
void dc_stream_receive_lz::shutdown() { }

} // namespace dc_impl
} // namespace graphlab
