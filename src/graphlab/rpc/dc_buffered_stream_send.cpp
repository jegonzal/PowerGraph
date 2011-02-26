#include <iostream>
#include <boost/iostreams/stream.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_buffered_stream_send.hpp>

namespace graphlab {
namespace dc_impl {

void dc_buffered_stream_send::send_data(procid_t target_, 
                unsigned char packet_type_mask,
                std::istream &istrm,
                size_t len) {
  ASSERT_EQ(target, target_);
  if (len != size_t(-1)) {
    char cbuffer[len];
    while(len > 0 && istrm.good()) {
      size_t l = istrm.readsome(cbuffer, len);
      len -= l;
    }
    send_data(target, packet_type_mask, cbuffer, len);
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
    free(data);
  }
}

void dc_buffered_stream_send::send_data(procid_t target, 
                 unsigned char packet_type_mask,
                 char* data, size_t len) {
  if ((packet_type_mask & CONTROL_PACKET) == 0) {
    if (packet_type_mask & (FAST_CALL | STANDARD_CALL)) {
      dc->inc_calls_sent(target);
    }
    bytessent.inc(len);
  }

  // build the packet header
  packet_hdr hdr;
  memset(&hdr, 0, sizeof(packet_hdr));
  
  hdr.len = len;
  hdr.src = dc->procid(); 
  hdr.sequentialization_key = dc->get_sequentialization_key();
  hdr.packet_type_mask = packet_type_mask;
  
  std::streamsize numbytes_needed = sizeof(packet_hdr) + len;
  // write into the send buffer
  sendbuf.begin_critical_section();
  // is there room?
  if (numbytes_needed <= sendbuf.free_space()) {
    // buffer has room to hold it. everything is all right with this world
    sendbuf.write_unsafe(reinterpret_cast<char*>(&(hdr)), sizeof(packet_hdr));
    sendbuf.write_unsafe(data, len);
    sendbuf.end_critical_section_with_signal();
  }
  else if (numbytes_needed < sendbuf.reserved_size()) {
    // the buffer is large enough. Just not enough room
    while (sendbuf.free_space() < numbytes_needed) {
      if (sendbuf.reader_is_blocked()) {
        // sender is definitely inside a lock
        // I take over the responsibilty of sending...
        send_till_empty();
      }
      else {
        // sender may be still sending.
        // lets wait till it is done
        sched_yield();
      }
    }
    sendbuf.write_unsafe(reinterpret_cast<char*>(&(hdr)), sizeof(packet_hdr));
    sendbuf.write_unsafe(data, len);      
    sendbuf.end_critical_section_with_signal();

  }
  else {
    // the buffer is not large enough
    // wait for the send buffer to clear
    // the buffer is large enough. Just not enough room
    while (sendbuf.size() > 0) {
      if (sendbuf.reader_is_blocked()) {
        // sender is definitely inside a lock
        // I take over the responsibilty of sending...
        send_till_empty();
      }
      else {
        // sender may be still sending.
        // lets wait till it is done
        sched_yield();
      }
    }
    comm->send2(target, 
                reinterpret_cast<char*>(&hdr), sizeof(packet_hdr),
                data,  len);
    sendbuf.end_critical_section();
  }
  
  

}

void dc_buffered_stream_send::send_till_empty() {
  while(1) {
    char* c = NULL;
    std::streamsize readlen = sendbuf.introspective_read(c, 65536);
    if (readlen == 0) break;
    comm->send(target, c,  readlen);
    sendbuf.advance_head(readlen);
  }
}

void dc_buffered_stream_send::send_loop() {
  while (!done) {
    char *c = NULL;
    std::streamsize readlen = sendbuf.blocking_introspective_read(c, 65536);
    if (readlen == 0 && sendbuf.is_done()) break;
    comm->send(target, c,  readlen);
    sendbuf.advance_head(readlen);
  }
}

void dc_buffered_stream_send::shutdown() {
  done = true;
  sendbuf.stop_reader();
  thr.join();
}

} // namespace dc_impl
} // namespace graphlab


