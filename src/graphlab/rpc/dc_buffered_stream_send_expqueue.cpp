#include <iostream>
#include <boost/iostreams/stream.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_buffered_stream_send_expqueue.hpp>

namespace graphlab {
namespace dc_impl {

void dc_buffered_stream_send_expqueue::send_data(procid_t target_, 
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

void dc_buffered_stream_send_expqueue::send_data(procid_t target, 
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
  hdr.packet_type_mask = packet_type_mask;
  
  std::streamsize numbytes_needed = sizeof(packet_hdr) + len;
  expqueue_entry eentry;
  eentry.len = numbytes_needed;
  eentry.c = (char*)malloc(numbytes_needed);
  memcpy(eentry.c, &hdr, sizeof(packet_hdr));
  memcpy(eentry.c + sizeof(packet_hdr), data, len);
  sendqueue.enqueue(eentry);
}


void dc_buffered_stream_send_expqueue::write_combining_send(expqueue_entry e) {
  char sendcombining[combine_upper_threshold];
  size_t length = e.len;
  memcpy(sendcombining, e.c, e.len);
  free(e.c);
  
  while(1) {
    std::pair<expqueue_entry, bool> data = sendqueue.try_dequeue();
    if (data.second) {
      // has data
      if (length + data.first.len <= combine_upper_threshold) {
        // if resultant length is still below threshold, write into the combiner
        memcpy(sendcombining + length, data.first.c , data.first.len);
        length += data.first.len;
        free(data.first.c);
      }
      else {
        // send the combined data
        comm->send(target, sendcombining, length);
        // and send what we just read
        comm->send(target, data.first.c, data.first.len);
        free(data.first.c);
        break;
      }
    }
    else {
      // no more data. just send what we have
      comm->send(target, sendcombining, length);
      break;
    }
  }
}

void dc_buffered_stream_send_expqueue::send_loop() {
  
  while (1) {
    std::pair<expqueue_entry, bool> data = sendqueue.dequeue();
    if (data.second == false) break;
    
    // if the length is small (below the combining threshold
    // and if it fits in the combining buffer. Write it into the buffer
    if (data.first.len <= combine_lower_threshold) {
      write_combining_send(data.first);        
    }
    else {
      comm->send(target, data.first.c, data.first.len);
      free(data.first.c);
    }
  }
}

void dc_buffered_stream_send_expqueue::shutdown() {
  sendqueue.stop_blocking();
  thr.join();
}

} // namespace dc_impl
} // namespace graphlab


