/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_stream_receive_z.hpp>
#include <zlib.h>
//#define DC_RECEIVE_DEBUG

namespace graphlab {
namespace dc_impl {


/**
  Called by the controller when there is data coming
  from the source
*/
void dc_stream_receive_z::incoming_data(procid_t src, 
                    const char* buf, 
                    size_t len) {
  zstrm.avail_in = (uInt)len;
  zstrm.next_in = (Bytef*)buf;
 
  compressed_bytesreceived.inc(len); 
  do{
    zstrm.avail_out = 128*1024;
    zstrm.next_out = (Bytef*)zbuffer;
    ASSERT_NE(inflate(&zstrm, Z_SYNC_FLUSH), Z_STREAM_ERROR);
    bufferlock.lock();
    buffer.write(zbuffer, 128*1024 - zstrm.avail_out);
    bufferlock.unlock();
  } while(zstrm.avail_out == 0);
  process_buffer(false);

  //std::cout << dc->procid() << ": From " << src << ": " << zstrm.total_in << " --> " << zstrm.total_out << std::endl; 
}
  
/** called by the controller when a function
call is completed */
void dc_stream_receive_z::function_call_completed(unsigned char packettype) {
  size_t pending = pending_calls.dec();
  if (barrier && pending == 0) {
    bufferlock.lock();
    barrier = false;
    bufferlock.unlock();
    process_buffer(false);
  }
}
void dc_stream_receive_z::process_buffer(bool outsidelocked) {
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

char* dc_stream_receive_z::get_buffer(size_t& retbuflength) {
  char* ret;
  bufferlock.lock();
  // get a write section
  retbuflength = buffer.introspective_write(ret);
  assert(retbuflength > 0);
  bufferlock.unlock();
  return ret;
}


char* dc_stream_receive_z::advance_buffer(char* c, size_t wrotelength, 
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


size_t dc_stream_receive_z::bytes_received() {
  return bytesreceived;
}
  
void dc_stream_receive_z::shutdown() { }

} // namespace dc_impl
} // namespace graphlab
