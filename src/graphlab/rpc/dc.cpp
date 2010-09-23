#include <boost/unordered_map.hpp>
#include <boost/bind.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc_sctp_comm.hpp>

#include <graphlab/rpc/dc_stream_send.hpp>
#include <graphlab/rpc/dc_stream_receive.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>

namespace graphlab {

 /**
Callback function. This function is called whenever data is received
*/
void dc_recv_callback(void* tag, procid_t src, const char* buf, size_t len) {
  distributed_control *dc = (distributed_control*)(tag);
  dc->receivers[src]->incoming_data(src, buf, len);
}

distributed_control::~distributed_control() {
  logstream(LOG_INFO) << "Shutting down distributed control " << std::endl;
  comm->close();
  for (size_t i = 0;i < senders.size(); ++i) delete senders[i];
  for (size_t i = 0;i < receivers.size(); ++i) delete receivers[i];
  senders.clear();
  receivers.clear();
  // shutdown function call handlers
  fcallqueue.stop_blocking();
  fcallhandlers.join();
  delete comm;
}
  
void distributed_control::exec_function_call(procid_t source, std::istream &istrm) {
  // extract the dispatch function
  iarchive arc(istrm);
  size_t f; 
  arc >> f;
  // a regular funcion call
  if (f != 0) {
    dc_impl::dispatch_type dispatch = (dc_impl::dispatch_type)f;
    dispatch(*this, source, istrm);
  }
  else {
    // f is NULL!. This is a portable call. deserialize the function name
    std::string s;
    arc >> s;
    char isrequest;
    arc >> isrequest;
    if (isrequest == 0) {
      std::cout << "portable call to " << s << std::endl;
      // look for the registration
      dc_impl::dispatch_map_type::const_iterator iter = portable_dispatch_call_map.find(s);
      if (iter == portable_dispatch_call_map.end()) {
        logstream(LOG_ERROR) << "Unable to locate dispatcher for function " << s << std::endl;
        return;
      }
      // dispatch
      iter->second(*this, source, istrm);

    }
    else {
     std::cout << "portable request to " << s << std::endl;
     dc_impl::dispatch_map_type::const_iterator iter = portable_dispatch_request_map.find(s);
      if (iter == portable_dispatch_request_map.end()) {
        logstream(LOG_ERROR) << "Unable to locate dispatcher for function " << s << std::endl;
        return;
      }
      // dispatch
      iter->second(*this, source, istrm);

    }
  }
} 
 
 
void distributed_control::deferred_function_call(procid_t source, char* buf, size_t len) {
  fcallqueue.enqueue(function_call_block(source, buf, len));
}

void distributed_control::fcallhandler_loop() {
  // pop an element off the queue
  while(1) {
    std::pair<function_call_block, bool> entry;
    entry = fcallqueue.dequeue();
    // if the queue is empty and we should quit
    if (entry.second == false) return;
    
    //create a stream containing all the data
    boost::iostreams::stream<boost::iostreams::array_source> 
                                istrm(entry.first.data, entry.first.len);
    exec_function_call(entry.first.source, istrm);
    receivers[entry.first.source]->function_call_completed();
  }
}


void distributed_control::init(const std::vector<std::string> &machines,
            const std::string &initstring,
            procid_t curmachineid,
            size_t numhandlerthreads,
            dc_comm_type commtype) {   
  REGISTER_RPC((*this), reply_increment_counter);
  if (commtype == TCP_COMM) {
    comm = new dc_impl::dc_tcp_comm();
  }
  else if (commtype == SCTP_COMM) {
    #ifdef DHAS_SCTP
    comm = new dc_impl::dc_sctp_comm();
    #else
    logger(LOG_FATAL, "SCTP support was not compiled");
    #endif
  }
  else {
    ASSERT_MSG(false, "Unexpected value for comm type");
  }
  // create the receiving objects
  if (comm->capabilities() && dc_impl::COMM_STREAM) {
    for (size_t i = 0; i < machines.size(); ++i) {
      receivers.push_back(new dc_impl::dc_stream_receive(this));
      senders.push_back(new dc_impl::dc_stream_send(this, comm));
    }
  }
  else {
    // TODO
    logstream(LOG_FATAL) << "Datagram handlers not implemented yet" << std::endl;
  }

  // create the handler threads
  // store the threads in the threadgroup
  for (size_t i = 0;i < numhandlerthreads; ++i) {
    launch_in_new_thread(fcallhandlers,
                          boost::bind(&distributed_control::fcallhandler_loop, 
                                      this));
  }

  // start the machines
  comm->init(machines, initstring, curmachineid, 
            dc_recv_callback, this); 
  
  }

}
