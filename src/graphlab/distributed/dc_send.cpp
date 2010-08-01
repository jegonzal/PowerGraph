#include <sys/types.h>
#include <sys/socket.h>
//#include <aio.h>


#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/dc_internal.hpp>

namespace graphlab {
  

void distributed_control::send_to_sock(procid_t target, char* buf, size_t len) {
  pheaderlen_t *p = (pheaderlen_t*)buf;
  (*p) = len;
  //sendlocks[target].lock();
  // logger(LOG_INFO, "Sending %d", len);
  int numsent = 0;
  int totallen = len + sizeof(pheaderlen_t);
  while (numsent < totallen) {
    int ret = send(socks[target], buf + numsent, totallen - numsent, 0);
    if (ret < 0) {
      logstream(LOG_FATAL) << "send error: " << strerror(errno) << std::endl;
      ASSERT_TRUE(false);
    }
    numsent += ret;
  }
  //sendlocks[target].unlock();
  free(buf);
}

void distributed_control::send_call_message(procid_t target, 
                                         void* remote_function, 
                                         void* ptr, size_t len,
                                         size_t numargs, handlerarg_t *arr) {
  msgsent.inc();
  DCHECK_NE(target, id);
  DCHECK_LE(len, std::numeric_limits<uint32_t>::max());
  DCHECK_LE(numargs, 8);

  size_t packlength = sizeof(remotecall_packdata) + 
                      sizeof(handlerarg_t) * numargs + len;

  char * sendbuffer = (char*)malloc(packlength + sizeof(pheaderlen_t));
  char * realbufpos = sendbuffer + sizeof(pheaderlen_t);
  ASSERT_NE(sendbuffer, NULL);

  // pack the data
  remotecall_packdata* p = (remotecall_packdata*)(realbufpos);
  memset(p, 0, sizeof(remotecall_packdata));
  p->packtype = REMOTECALL_ID;
  p->srcnodeid = id;
  p->fnptr = remote_function;
  p->len = uint32_t(len);
  p->numargs = numargs;
  memcpy(p->args, arr, sizeof(handlerarg_t) * numargs);
  
  // now to pack in the rest of the blob
  size_t blobbegin = packlength - len;
  memcpy(realbufpos + blobbegin, ptr, len);
  
  // send! print output on failure and quit
  //logger(LOG_INFO, "Total Sent %d", ctr);
  //send_to_sock(target, sendbuffer, packlength);
  send_requests.enqueue(send_req_data(target, sendbuffer, packlength));
}


void distributed_control::send_callx_message(procid_t target, void* remote_function, 
                        void* ptr, size_t len,void* stackbegin, size_t stacklen) {
  msgsent.inc();
  DCHECK_NE(target, id);
  DCHECK_LE(len, std::numeric_limits<uint32_t>::max());


  // compute the message size 
  size_t packlength = sizeof(remotecallx_packdata) +len +stacklen;
  // allocate    
  char* sendbuffer = (char*)malloc(packlength + sizeof(pheaderlen_t));
  char* realbufpos = sendbuffer + sizeof(pheaderlen_t);
  ASSERT_NE(sendbuffer, NULL);

  // pack the data 
  remotecallx_packdata* p =  (remotecallx_packdata*)(realbufpos);
  memset(p, 0, sizeof(remotecallx_packdata));
  p->packtype = REMOTECALLX_ID;
  p->srcnodeid = id;
  p->fnptr = remote_function;
  p->len = uint32_t(len);
  p->stacklen = uint32_t(stacklen);
  memcpy(realbufpos + sizeof(remotecallx_packdata), ptr, len);
  memcpy(realbufpos + sizeof(remotecallx_packdata) +len ,
         stackbegin, stacklen);
         

  send_requests.enqueue(send_req_data(target, sendbuffer, packlength));
  
}

void distributed_control::send_callxs_message(procid_t target, void* remote_function, 
                        void* ptr, size_t len,const std::string &stack) {
  msgsent.inc();
  DCHECK_NE(target, id);
  DCHECK_LE(len, std::numeric_limits<uint32_t>::max());

  // compute the message size 
  size_t packlength = sizeof(remotecallxs_packdata) + len + stack.length();
  // allocate    
  char* sendbuffer = (char*)malloc(packlength + sizeof(pheaderlen_t));
  char* realbufpos = sendbuffer + sizeof(pheaderlen_t);
  ASSERT_NE(sendbuffer, NULL);

  // pack the data 
  remotecallxs_packdata* p =  (remotecallxs_packdata*)(realbufpos);
  memset(p, 0, sizeof(remotecallxs_packdata));
  p->packtype = REMOTECALLXS_ID;
  p->srcnodeid = id;
  p->fnptr = remote_function;
  p->len = uint32_t(len);
  p->stacklen = uint32_t(stack.length());
  memcpy(realbufpos + sizeof(remotecallx_packdata), ptr, len);
  memcpy(realbufpos + sizeof(remotecallx_packdata) +len ,
         stack.c_str(), stack.length());

  send_requests.enqueue(send_req_data(target, sendbuffer, packlength));
  
}




void distributed_control::send_call_control_message(procid_t target, void* remote_function,
                        void* ptr, size_t len,void* stackbegin, size_t stacklen) {
  DCHECK_NE(target, id);
  DCHECK_LE(len, std::numeric_limits<uint32_t>::max());


  // compute the message size
  size_t packlength = sizeof(remotecallx_packdata) +len +stacklen;
  // allocate
  char* sendbuffer = (char*)malloc(packlength + sizeof(pheaderlen_t));
  char* realbufpos = sendbuffer + sizeof(pheaderlen_t);
  ASSERT_NE(sendbuffer, NULL);

  // pack the data
  remotecallx_packdata* p =  (remotecallx_packdata*)(realbufpos);
  memset(p, 0, sizeof(remotecallx_packdata));
  p->packtype = REMOTECALL_CONTROL_ID;
  p->srcnodeid = id;
  p->fnptr = remote_function;
  p->len = uint32_t(len);
  p->stacklen = uint32_t(stacklen);
  memcpy(realbufpos + sizeof(remotecallx_packdata), ptr, len);
  memcpy(realbufpos + sizeof(remotecallx_packdata) +len ,
         stackbegin, stacklen);


  send_requests.enqueue(send_req_data(target, sendbuffer, packlength));

}


}
