#include <unistd.h> 
#include <sys/param.h>
#include <stdlib.h>
#include <sys/types.h>
#include <ifaddrs.h>
#include <netinet/in.h>

#include <map>
#include <sstream>

#include <boost/unordered_map.hpp>
#include <boost/bind.hpp>

#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/net_util.hpp>
#include <graphlab/metrics/metrics.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc_sctp_comm.hpp>

#include <graphlab/rpc/dc_stream_send.hpp>
#include <graphlab/rpc/dc_stream_receive.hpp>
#include <graphlab/rpc/dc_buffered_stream_send.hpp>
#include <graphlab/rpc/dc_buffered_stream_send_expqueue.hpp>
#include <graphlab/rpc/dc_buffered_stream_receive.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/dc_services.hpp>

namespace graphlab {

namespace dc_impl {
bool thrlocal_resizing_array_key_initialized = false;
pthread_key_t thrlocal_resizing_array_key;

struct dc_tls_data{
  resizing_array_sink ras;
  boost::iostreams::stream<resizing_array_sink_ref> strm;    
  
  dc_tls_data(): ras(128),strm(ras){ };
  ~dc_tls_data() {
    strm.flush();
    free(ras.c_str());
  }
  
};

boost::iostreams::stream<resizing_array_sink_ref>& 
get_thread_local_stream() {
  dc_tls_data* curptr = reinterpret_cast<dc_tls_data*>(
                        pthread_getspecific(thrlocal_resizing_array_key));
  if (curptr != NULL) {
    curptr->ras.clear();
    return curptr->strm;
  }
  else {
    dc_tls_data* ras = new dc_tls_data;
    int err = pthread_setspecific(thrlocal_resizing_array_key, ras);
    ASSERT_EQ(err, 0);
    return ras->strm;
  }
}

void thrlocal_destructor(void* v){ 
  dc_tls_data* s = reinterpret_cast<dc_tls_data*>(v);
  if (s != NULL) {
    delete s;
  }
  pthread_setspecific(thrlocal_resizing_array_key, NULL);
}
}

static std::string get_working_dir() {
#ifdef _GNU_SOURCE
  char* path = get_current_dir_name();
  assert(path != NULL);
  std::string ret = path;
  free(path);
#else
  char path[MAXPATHLEN];
  assert(getcwd(path, MAXPATHLEN) != NULL);
  std::string ret = path;
#endif
  return ret;
}

 /**
Callback function. This function is called whenever data is received
*/
void dc_recv_callback(void* tag, procid_t src, const char* buf, size_t len) {
  distributed_control *dc = (distributed_control*)(tag);
  dc->receivers[src]->incoming_data(src, buf, len);
}

distributed_control::~distributed_control() {
  distributed_services->barrier();
  logstream(LOG_INFO) << "Shutting down distributed control " << std::endl;
  size_t bytessent = bytes_sent();
  for (size_t i = 0;i < senders.size(); ++i) {
    senders[i]->shutdown();
    delete senders[i];
  }
  comm->close();
  size_t bytesreceived = bytes_received();
  for (size_t i = 0;i < receivers.size(); ++i) {
    receivers[i]->shutdown();
    delete receivers[i];
  }
  senders.clear();
  receivers.clear();
  // shutdown function call handlers
  fcallqueue.stop_blocking();
  fcallhandlers.join();
  delete comm;
  logstream(LOG_INFO) << "Bytes Sent: " << bytessent << std::endl;
  logstream(LOG_INFO) << "Calls Sent: " << calls_sent() << std::endl;
  logstream(LOG_INFO) << "Bytes Received: " << bytesreceived << std::endl;
  logstream(LOG_INFO) << "Calls Received: " << calls_received() << std::endl;

}
  
void distributed_control::exec_function_call(procid_t source, 
                                            unsigned char packet_type_mask, 
                                            std::istream &istrm) {
  // extract the dispatch function
  iarchive arc(istrm);
  size_t f; 
  arc >> f;
  // a regular funcion call
  if (f != 0) {
    dc_impl::dispatch_type dispatch = (dc_impl::dispatch_type)f;
    dispatch(*this, source, packet_type_mask, istrm);
  }
  else {
    // f is NULL!. This is a portable call. deserialize the function name
    std::string s;
    arc >> s;
    char isrequest;
    arc >> isrequest;
    if (isrequest == 0) {
      // std::cout << "portable call to " << s << std::endl;
      // look for the registration
      dc_impl::dispatch_map_type::const_iterator iter = portable_dispatch_call_map.find(s);
      if (iter == portable_dispatch_call_map.end()) {
        logstream(LOG_ERROR) << "Unable to locate dispatcher for function " << s << std::endl;
        return;
      }
      // dispatch
      iter->second(*this, source, packet_type_mask, istrm);

    }
    else {
     // std::cout << "portable request to " << s << std::endl;
     dc_impl::dispatch_map_type::const_iterator iter = portable_dispatch_request_map.find(s);
      if (iter == portable_dispatch_request_map.end()) {
        logstream(LOG_ERROR) << "Unable to locate dispatcher for function " << s << std::endl;
        return;
      }
      // dispatch
      iter->second(*this, source, packet_type_mask, istrm);

    }
  }
  if ((packet_type_mask & CONTROL_PACKET) == 0) inc_calls_received(source);
} 
 
 
void distributed_control::deferred_function_call(procid_t source, unsigned char packet_type_mask,
                                                char* buf, size_t len) {
  fcallqueue.enqueue(function_call_block(source, packet_type_mask, buf, len));
}

void distributed_control::fcallhandler_loop(size_t id) {
  // pop an element off the queue
  while(1) {
    std::pair<function_call_block, bool> entry;
    entry = fcallqueue.dequeue(id);
    // if the queue is empty and we should quit
    if (entry.second == false) return;
    
    //create a stream containing all the data
    boost::iostreams::stream<boost::iostreams::array_source> 
                                istrm(entry.first.data, entry.first.len);
    exec_function_call(entry.first.source, entry.first.packet_type_mask, istrm);
    receivers[entry.first.source]->function_call_completed(entry.first.packet_type_mask);
  }
  std::cerr << "Handler " << id << " died." << std::endl;
}


std::map<std::string, std::string> distributed_control::parse_options(std::string initstring) {
  std::map<std::string, std::string> options;
  std::replace(initstring.begin(), initstring.end(), ',', ' ');
  std::replace(initstring.begin(), initstring.end(), ';', ' ');        
  std::string opt, value;
  // read till the equal
  std::stringstream s(initstring);
  while(s.good()) {
    getline(s, opt, '=');
    if (s.bad() || s.eof()) break;
    getline(s, value, ' ');
    if (s.bad()) break;
    options[trim(opt)] = trim(value);
  }
  return options;
}

void distributed_control::init(const std::vector<std::string> &machines,
            const std::string &initstring,
            procid_t curmachineid,
            size_t numhandlerthreads,
            dc_comm_type commtype) {   
  // initialize thread local storage
  if (dc_impl::thrlocal_resizing_array_key_initialized == false) {
    dc_impl::thrlocal_resizing_array_key_initialized = true;
    int err = pthread_key_create(&dc_impl::thrlocal_resizing_array_key, 
                        dc_impl::thrlocal_destructor);
    ASSERT_EQ(err, 0);
  }
  //-------- Initialize the full barrier ---------
  full_barrier_in_effect = false;
  procs_complete.resize(machines.size());
  //-----------------------------------------------
  
  REGISTER_RPC((*this), reply_increment_counter);
  // parse the initstring
  std::map<std::string,std::string> options = parse_options(initstring);
  bool buffered_send = false;
  bool buffered_recv = false;
  bool buffered_queued_send = false;
  if (options["buffered_send"] == "true" || 
    options["buffered_send"] == "1" ||
    options["buffered_send"] == "yes") {
    buffered_send = true;
    std::cerr << "Buffered Send Option is ON." << std::endl;
  }

  if (options["buffered_queued_send"] == "true" || 
    options["buffered_queued_send"] == "1" ||
    options["buffered_queued_send"] == "yes") {
    buffered_queued_send = true;
    if (buffered_send == true) {
      std::cerr << "buffered_queued_send and buffered_send cannot be on simultaneously" << std::endl;
      exit(1);
    }
    std::cerr << "Buffered Queued Send Option is ON." << std::endl;
  }

  if (options["buffered_recv"] == "true" ||
    options["buffered_recv"] == "1" ||
    options["buffered_recv"] == "yes") {
    buffered_recv = true;
    std::cerr << "Buffered Recv Option is ON." << std::endl;
  }
  
  if (commtype == TCP_COMM) {
    comm = new dc_impl::dc_tcp_comm();
    std::cerr << "TCP Communication layer constructed." << std::endl;
  }
  else if (commtype == SCTP_COMM) {
    #ifdef HAS_SCTP
    comm = new dc_impl::dc_sctp_comm();
    std::cerr << "SCTP Communication layer constructed." << std::endl;
    #else
    logger(LOG_FATAL, "SCTP support was not compiled");
    #endif
  }
  else {
    ASSERT_MSG(false, "Unexpected value for comm type");
  }
  global_calls_sent.resize(machines.size());
  global_calls_received.resize(machines.size());
  fcallqueue.init(numhandlerthreads);
  // create the receiving objects
  if (comm->capabilities() && dc_impl::COMM_STREAM) {
    for (size_t i = 0; i < machines.size(); ++i) {
      if (buffered_recv) {
        receivers.push_back(new dc_impl::dc_buffered_stream_receive(this));
      }
      else {
        receivers.push_back(new dc_impl::dc_stream_receive(this));
      }
      
      if (buffered_send) {
        senders.push_back(new dc_impl::dc_buffered_stream_send(this, comm, i));
      }
      else if (buffered_queued_send) {
        senders.push_back(new dc_impl::dc_buffered_stream_send_expqueue(this, comm, i));
      }
      else{
        senders.push_back(new dc_impl::dc_stream_send(this, comm, i));
      }
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
                                      this, i));
  }

  // start the machines
  comm->init(machines, options, curmachineid, 
              receivers); 
  
  // set the local proc values
  localprocid = comm->procid();
  localnumprocs = comm->numprocs();
  
  // construct the services
  distributed_services = new dc_services(*this);
  compute_master_ranks();
}
  
dc_services& distributed_control::services() {
  return *distributed_services;
}


void distributed_control::comm_barrier(procid_t targetmachine) {
  ASSERT_LT(targetmachine, numprocs());
  if (targetmachine != procid() && senders[targetmachine]->channel_active(targetmachine)) {
    std::stringstream strm;
    senders[targetmachine]->send_data(targetmachine, BARRIER | CONTROL_PACKET, strm, 0);
  }
}

void distributed_control::comm_barrier() {
  std::stringstream strm;
  for (size_t i = 0;i < senders.size(); ++i) {
    if (i != procid() && senders[i]->channel_active(i)) {
      senders[i]->send_data(i, BARRIER | CONTROL_PACKET, strm, 0);
    }
  }
}

void distributed_control::barrier() {
  distributed_services->barrier();
}



void distributed_control::compute_master_ranks() {
  uint32_t local_ip = get_local_ip();
  std::string localipandpath = tostr(local_ip) + get_working_dir();
  std::vector<std::string> ipandpath(numprocs());
  ipandpath[procid()] = localipandpath;
  all_gather(ipandpath);
  for(size_t i = 0; i < ipandpath.size(); ++i) {
    if(ipandpath[i] == localipandpath) {
      masterid = i;
      break;
    }
  }
}

 /*****************************************************************************
                      Implementation of Full Barrier
*****************************************************************************/
/* It is unfortunate but this is copy paste code from dc_dist_object.hpp
  I thought for a long time how to implement this without copy pasting and 
  I can't think of a simple enough solution.
  
  Part of the issue is that the "context" concept was not built into to the
  RPC system to begin with and is currently folded in through the dc_dist_object system.
  As a result, the global context becomes very hard to define properly.
  Including a dc_dist_object as a member only resolves the high level contexts
  such as barrier, broadcast, etc which do not require intrusive access into
  deeper information about the context. The full barrier however, requires deep
  information about the context which cannot be resolved easily.
*/

/**
This barrier ensures globally across all machines that
all calls issued prior to this barrier are completed before
returning. This function could return prematurely if
other threads are still issuing function calls since we
cannot differentiate between calls issued before the barrier
and calls issued while the barrier is being evaluated.
*/
void distributed_control::full_barrier() {
  // gather a sum of all the calls issued to machine 0
  std::vector<size_t> calls_sent_to_target(numprocs(), 0);
  for (size_t i = 0;i < numprocs(); ++i) {
    calls_sent_to_target[i] = global_calls_sent[i].value;
  }
  
  // tell node 0 how many calls there are
  std::vector<std::vector<size_t> > all_calls_sent(numprocs());
  all_calls_sent[procid()] = calls_sent_to_target;
  all_gather(all_calls_sent, true);
  
  // get the number of calls I am supposed to receive from each machine
  calls_to_receive.clear(); calls_to_receive.resize(numprocs(), 0);
  for (size_t i = 0;i < numprocs(); ++i) {
    calls_to_receive[i] += all_calls_sent[i][procid()];
//    std::cout << "Expecting " << calls_to_receive[i] << " calls from " << i << std::endl;
  }
  // clear the counters
  num_proc_recvs_incomplete.value = numprocs();
  procs_complete.clear();
  // activate the full barrier
  full_barrier_in_effect = true;
  // begin one pass to set all which are already completed
  for (size_t i = 0;i < numprocs(); ++i) {
    if (global_calls_received[i].value >= calls_to_receive[i]) {
      if (procs_complete.set_bit(i) == false) {
        num_proc_recvs_incomplete.dec();
      }
    }
  }
  
  full_barrier_lock.lock();
  while (num_proc_recvs_incomplete.value > 0) full_barrier_cond.wait(full_barrier_lock);
  full_barrier_lock.unlock();
  full_barrier_in_effect = false;
//   for (size_t i = 0; i < numprocs(); ++i) {
//     std::cout << "Received " << global_calls_received[i].value << " from " << i << std::endl;
//   }
}

std::map<std::string, size_t> distributed_control::gather_statistics(){
    std::map<std::string, size_t> ret;

    std::vector<collected_statistics> stats(numprocs());
    stats[procid()].callssent = calls_sent();
    stats[procid()].bytessent = bytes_sent();
    
    gather(stats, 0, true);
    if (procid() == 0) {
      collected_statistics cs;
      for (size_t i = 0;i < numprocs(); ++i) {
        cs.callssent += stats[i].callssent;
        cs.bytessent += stats[i].bytessent;
      }
      ret["total_calls_sent"] = cs.callssent;
      ret["total_bytes_sent"] = cs.bytessent;
    }
    return ret; 
}

void distributed_control::fill_metrics() {
  std::map<std::string, size_t> ret = gather_statistics();
  if (procid() == 0) {
    metrics& engine_metrics = metrics::create_metrics_instance("RPC", true);
    engine_metrics.set("nodes", numprocs(), INTEGER);
    engine_metrics.set("total_calls_sent", ret["total_calls_sent"], INTEGER);
    engine_metrics.set("total_bytes_sent", ret["total_bytes_sent"], INTEGER);
  }
}

} //namespace graphlab
