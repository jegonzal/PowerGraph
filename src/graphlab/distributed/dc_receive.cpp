#include <sys/types.h>
#include <sys/socket.h>
#include <aio.h>


#include <poll.h>



#include <iostream>
#include <string>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>


#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/dc_internal.hpp>

namespace graphlab {
// background messaging thread

void distributed_control::message_dispatch_thread::run(distributed_control *dc_) {
  distributed_control &dc = *dc_;
  logger(LOG_INFO, "message thread started");
  while(1) {
    std::pair<dispatch_req_data, bool> ret = dc.dispatch_requests.dequeue();
    if (ret.second == false) break;

    uint16_t packtype = *(uint16_t*)(ret.first.buf);
    switch(packtype){
      case REMOTECALL_ID:
        dc.receive_call_message(ret.first.buf, ret.first.len);
        break;
      case REMOTECALLX_ID:
        dc.receive_callx_message(ret.first.buf, ret.first.len);
        break;
      case REMOTECALLXS_ID:
        dc.receive_callxs_message(ret.first.buf, ret.first.len);
        break;
      case REMOTECALL_CONTROL_ID:
        dc.receive_callx_message(ret.first.buf, ret.first.len);
        break;      
      default:
        ASSERT_MSG(false, "Invalid Packet Type %d", packtype);
    }
    free(ret.first.buf);
    if (packtype != REMOTECALL_CONTROL_ID) dc.msgprocessed.inc();
  }
}

  
  
// The background receiving thread

distributed_control::messageproc_thread::messageproc_thread() {
}


distributed_control::messageproc_thread::~messageproc_thread() {

}

void distributed_control::messageproc_thread::run() {
  // build the selects

  logger(LOG_INFO, "Message Processing Thread %d/%d started",
                    thread::thread_id() + 1, dc->procs.size());
  size_t j = 0;
  size_t sockid = 0;
  for (size_t i = 0;i < dc->socks.size(); ++i) {
    // if it is not to my self
    if (i != dc->procid()) {
      if (j == thread::thread_id()) {
        sockid = i;
        break;
      }
      ++j;
    }
  }

  while(!(*done)) {
    dc->receive_handler(dc->socks[sockid], dc->buffer[sockid]);
  }

  logger(LOG_INFO, "Message Processing Thread stopped");
}


void distributed_control::receive_handler(sockfd_t sock,recv_buffer &buffer) {
  while(1) {
    if (done) return;
    // receive the message
    ASSERT_GT((int)(buffer.buflen) - (int)(buffer.buftail), 0);
    int msglen = recv(sock, buffer.buffer + buffer.buftail, buffer.buflen - buffer.buftail, 0);
    // check for various failure conditions
    if (msglen < 0) {
      if (errno == EINTR) {
        continue;
      }
      else if (errno == EAGAIN) {
        continue;
      }
      else {
        logstream(LOG_FATAL) << "Receive Error: " << strerror(errno) << std::endl;
        ASSERT_TRUE(false);
      }
    }
    if (done) return;
    else if (msglen == 0) {
      // nothing left to do
      // in fact this means connection closed
      ASSERT_MSG(false, "Unexpected connection close");
      
      return;
    }
    //logger(LOG_INFO, "Receiving %d on sock %d", msglen, sock);
    // push the new location of buftail
    buffer.buftail += msglen;
    ASSERT_LE(buffer.buftail, buffer.buflen);
    // this the desired packet length
    if (buffer.buftail < 4) continue;
    // Update the buffer counters
    pheaderlen_t packetlen = *reinterpret_cast<pheaderlen_t*>(buffer.buffer);

    buffer.num_recvcalls++;
    
    
    // if my buffer will not fit. resize the buffer
    // remember to factor in the header size
    size_t eff_packetsize = packetlen + sizeof(pheaderlen_t);
    if (buffer.buflen < eff_packetsize) {
      buffer.buffer = (char*)realloc(buffer.buffer, 2 * eff_packetsize);
      buffer.buflen = 2 * eff_packetsize;
    }

    ASSERT_GE(buffer.buflen, buffer.minimum_buflen);
    // while I have complete packets in my buffer
    size_t bufhead = 0;
    while(1) {
      // if there is not even room to read the header, break
      if (bufhead + sizeof(pheaderlen_t) > buffer.buftail) break;
      
      pheaderlen_t packetlen = *reinterpret_cast<pheaderlen_t*>(buffer.buffer + bufhead);
      size_t eff_packetsize = packetlen + sizeof(pheaderlen_t);
      //logger(LOG_INFO, "Receiving packet of size %d on sock %d", packetlen, sock);
      // do we have the entire packet?
      if (buffer.buftail - bufhead >= eff_packetsize) {
        //logger(LOG_INFO, "Dispatch");
        // yep we have room
        // shift the head
        bufhead += sizeof(pheaderlen_t);
        
        dispatch_req_data dispatchreq;
        dispatchreq.buf = (char*)malloc(packetlen);
        dispatchreq.len = packetlen;
        
        ASSERT_MSG(dispatchreq.buf != NULL, "Could not allocate memory in dc_receive! Errno=%d", errno);
        
        memcpy(dispatchreq.buf, buffer.buffer + bufhead, packetlen);
        dispatch_requests.enqueue(dispatchreq);

        // shift the bufhead
        bufhead += packetlen;
      }
      else {
        break;
      }
    }
    // roll up the buffer
    memmove(buffer.buffer, buffer.buffer + bufhead, buffer.buftail - bufhead);
    buffer.buftail = buffer.buftail - bufhead;
  }
}

                             
void distributed_control::receive_call_message(char* msg, size_t len) {
  remotecall_packdata* packd = (remotecall_packdata*)msg;
  // check the size
  size_t estpacklength = sizeof(remotecall_packdata) + 
                        sizeof(handlerarg_t) * packd->numargs + packd->len;
  ASSERT_EQ(estpacklength, len);
  ASSERT_LE(packd->numargs, 8);
  
  void* blobptr = NULL;
  if (packd->len > 0) blobptr = msg + len - packd->len;
  handler_type fn = (handler_type)(packd->fnptr);
  //logstream(LOG_INFO) << "Calling " << &fn << " from " << packd->srcnodeid << "\n";
  //logstream(LOG_INFO) << "bloblength " << packd->len <<". numargs " << packd->numargs << "\n";
  switch(packd->numargs) {
   case 0:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len); 
    break;
   case 1:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0]); 
    break;
   case 2:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0], packd->args[1]); 
    break;
   case 3:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0], packd->args[1],packd->args[2]); 
    break;
   case 4:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0], packd->args[1],packd->args[2],
          packd->args[3]); 
    break;
   case 5:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0], packd->args[1],packd->args[2],
          packd->args[3], packd->args[4]); 
    break;
   case 6:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0], packd->args[1],packd->args[2],
          packd->args[3], packd->args[4], packd->args[5]); 
    break;
   case 7:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0], packd->args[1],packd->args[2],
          packd->args[3], packd->args[4], packd->args[5],
          packd->args[6]); 
    break;
   case 8:
    (fn)(*this, packd->srcnodeid, blobptr, packd->len, 
          packd->args[0], packd->args[1],packd->args[2],
          packd->args[3], packd->args[4], packd->args[5],
          packd->args[6], packd->args[7]); 
    break;
   default:
    break;
  }
}

void distributed_control::receive_callx_message(char* msg, size_t len) {
  remotecallx_packdata* packd = (remotecallx_packdata*)msg;
  // check the size
  size_t estpacklength = sizeof(remotecallx_packdata) + 
                        packd->len + packd->stacklen;
  ASSERT_EQ(estpacklength, len);
  
  void* blobptr = NULL;
  void* stackptr = NULL;
  if (packd->len > 0) blobptr = msg + sizeof(remotecallx_packdata);
  stackptr = msg + sizeof(remotecallx_packdata) + packd->len;

  xdispatcher xdispatchptr = (xdispatcher)(packd->fnptr);
  (xdispatchptr)(*this, packd->srcnodeid, blobptr, packd->len, stackptr); 
}

void distributed_control::receive_callxs_message(char* msg, size_t len) {
  remotecallxs_packdata* packd = (remotecallxs_packdata*)msg;
  // check the size
  size_t estpacklength = sizeof(remotecallxs_packdata) + 
                        packd->len + packd->stacklen;
  ASSERT_EQ(estpacklength, len);
  
  void* blobptr = NULL;
  void* stackptr = NULL;
  if (packd->len > 0) blobptr = msg + sizeof(remotecallxs_packdata);
  stackptr = msg + sizeof(remotecallxs_packdata) + packd->len;
  boost::iostreams::stream<boost::iostreams::array_source>  streamin((char*)stackptr, packd->stacklen);
  iarchive arc(streamin);
  xsdispatcher xsdispatchptr = (xsdispatcher)(packd->fnptr);
  (xsdispatchptr)(*this, packd->srcnodeid, blobptr, packd->len, arc); 
}

  
}
