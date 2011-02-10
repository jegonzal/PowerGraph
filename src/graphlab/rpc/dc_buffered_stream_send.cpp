#include <iostream>
#include <boost/iostreams/stream.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_buffered_stream_send.hpp>

namespace graphlab {
namespace dc_impl {

/* Used to differentiate elements
 * in the send buffer */
struct buffer_head{
  size_t target;
  size_t len;
};
  

void dc_buffered_stream_send::send_data(procid_t target, 
                unsigned char packet_type_mask,
                std::istream &istrm,
                size_t len) {
  if (len != size_t(-1)) {
    // build the packet header
    if ((packet_type_mask & CONTROL_PACKET) == 0) {
      if (packet_type_mask & (FAST_CALL | STANDARD_CALL)) callssent.inc();
      bytessent.inc(len);
    }
    packet_hdr hdr;
    hdr.len = len;
    hdr.src = dc->procid(); 
    hdr.packet_type_mask = packet_type_mask;
    
    buffer_head bufhead;
    bufhead.target = target;
    bufhead.len = sizeof(packet_hdr) + len;
    sendbuflock.lock();
    sendbuf.write(reinterpret_cast<char*>(&(bufhead)), sizeof(buffer_head));
    sendbuf.write(reinterpret_cast<char*>(&(hdr)), sizeof(packet_hdr));
    
    char cbuffer[10240];
    size_t l = 0;
    while(istrm.good()) {
      l = istrm.readsome(cbuffer, 10240);
      if (l == 0) break;  // 0 return is empty buffer
      sendbuf.write(cbuffer, l);
    }
    sendcond.signal();
    sendbuflock.unlock();
  }
  else {
    // annoying. we have to compute the length of the stream
    // allocate a 128byte block first.
    // \todo: This can be optimized. Though, I don't think this
    //        code path is even used.
    size_t len = 0;
    size_t cursize = 128;
    char* data = (char*)malloc(128);
    // while the stream is good. read stuff
    // len is the current length of the contents
    // cursize is the max length of the array
    // when we run out of space in the array, we double the size.
    while (istrm.good()) {
      len += istrm.readsome(data+len, cursize-len);
      if (cursize - len == 0) {
        cursize *= 2;
        data = (char*)realloc(data, cursize);
      }
    }
    send_data(target, packet_type_mask, data, len);
  }
}

void dc_buffered_stream_send::send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len) {
  if ((packet_type_mask & CONTROL_PACKET) == 0) {
    if (packet_type_mask & (FAST_CALL | STANDARD_CALL)) {
      dc->inc_calls_sent();
      callssent.inc();
    }
    bytessent.inc(len);
  }

  // build the packet header
  packet_hdr hdr;
  hdr.len = len;
  hdr.src = dc->procid(); 
  hdr.packet_type_mask = packet_type_mask;
  
  buffer_head bufhead;
  bufhead.target = target;
  bufhead.len = sizeof(packet_hdr) + len;

  // write into the send buffer
  sendbuflock.lock();
  sendbuf.write(reinterpret_cast<char*>(&(bufhead)), sizeof(buffer_head));
  sendbuf.write(reinterpret_cast<char*>(&(hdr)), sizeof(packet_hdr));
  sendbuf.write(data, len);
  sendcond.signal();
  sendbuflock.unlock();
  
  free(data);
}

void dc_buffered_stream_send::send_loop() {
  sendbuflock.lock();
  while (!done) {
    while (sendbuf.size() != 0) {
      char* c;
      // get the size and target of this send
      buffer_head bufhead;
      sendbuf.read(reinterpret_cast<char*>(&bufhead), sizeof(buffer_head));
      
      while(bufhead.len > 0) {
        std::streamsize len = sendbuf.introspective_read(c, bufhead.len);
        comm->send(bufhead.target, c,  len);
        bufhead.len -= len;
      }
    }
    sendcond.wait(sendbuflock);
  }
  sendbuflock.unlock();
}

void dc_buffered_stream_send::shutdown() {
  done = true;
  sendcond.signal();
  thr.join();
}

} // namespace dc_impl
} // namespace graphlab


