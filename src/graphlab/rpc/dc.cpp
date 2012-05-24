/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#include <unistd.h> 
#include <sys/param.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <ifaddrs.h>
#include <netinet/in.h>

#include <map>
#include <sstream>

#include <boost/unordered_map.hpp>
#include <boost/bind.hpp>
//#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/net_util.hpp>
#include <graphlab/util/mpi_tools.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
//#include <graphlab/rpc/dc_sctp_comm.hpp>
#include <graphlab/rpc/dc_buffered_stream_send2.hpp>
#include <graphlab/rpc/dc_stream_receive.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/dc_services.hpp>

#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>

namespace graphlab {

namespace dc_impl {


bool thrlocal_resizing_array_key_initialized = false;
pthread_key_t thrlocal_resizing_array_key;

bool thrlocal_sequentialization_key_initialized = false;
pthread_key_t thrlocal_sequentialization_key;

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
    if (curptr->ras.buffer_size > 65536) curptr->ras.clear(1024);
    else curptr->ras.clear();
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

unsigned char distributed_control::set_sequentialization_key(unsigned char newkey) {
  size_t oldval = reinterpret_cast<size_t>(pthread_getspecific(dc_impl::thrlocal_sequentialization_key));
  size_t newval = newkey;
  pthread_setspecific(dc_impl::thrlocal_sequentialization_key, reinterpret_cast<void*>(newval));
  assert(oldval < 256);
  return (unsigned char)oldval;
}

unsigned char distributed_control::new_sequentialization_key() {
  size_t oldval = reinterpret_cast<size_t>(pthread_getspecific(dc_impl::thrlocal_sequentialization_key));
  size_t newval = (oldval + 1) % 256;
  pthread_setspecific(dc_impl::thrlocal_sequentialization_key, reinterpret_cast<void*>(newval));
  assert(oldval < 256);
  return (unsigned char)oldval;
}

unsigned char distributed_control::get_sequentialization_key() {
  size_t oldval = reinterpret_cast<size_t>(pthread_getspecific(dc_impl::thrlocal_sequentialization_key));
  assert(oldval < 256);
  return (unsigned char)oldval;
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


distributed_control::distributed_control() {
  dc_init_param initparam;
  if (init_param_from_env(initparam)) {
    logstream(LOG_INFO) << "Distributed Control Initialized from Environment" << std::endl;
  } else if (mpi_tools::initialized() && init_param_from_mpi(initparam)) {
      logstream(LOG_INFO) << "Distributed Control Initialized from MPI" << std::endl;
  }
  else {
    logstream(LOG_INFO) << "Shared Memory Execution" << std::endl;
    // get a port and socket
    std::pair<size_t, int> port_and_sock = get_free_tcp_port();
    size_t port = port_and_sock.first;
    int sock = port_and_sock.second;

    initparam.machines.push_back(std::string("localhost:") + tostr(port));
    initparam.curmachineid = 0;
    initparam.initstring = std::string(" __sockhandle__=") + tostr(sock) + " ";
    initparam.numhandlerthreads = RPC_DEFAULT_NUMHANDLERTHREADS;
    initparam.commtype = RPC_DEFAULT_COMMTYPE;
  }
  init(initparam.machines, 
        initparam.initstring, 
        initparam.curmachineid, 
        initparam.numhandlerthreads,
        initparam.commtype);
  INITIALIZE_TRACER(dc_receive_queuing, "dc: time spent on enqueue");
  INITIALIZE_TRACER(dc_receive_multiplexing, "dc: time spent exploding a chunk");
  INITIALIZE_TRACER(dc_call_dispatch, "dc: time spent issuing RPC calls");
}

distributed_control::distributed_control(dc_init_param initparam) {
  init(initparam.machines, 
        initparam.initstring, 
        initparam.curmachineid, 
        initparam.numhandlerthreads,
        initparam.commtype);
  INITIALIZE_TRACER(dc_receive_queuing, "dc: time spent on enqueue");
  INITIALIZE_TRACER(dc_receive_multiplexing, "dc: time spent exploding a chunk");
  INITIALIZE_TRACER(dc_call_dispatch, "dc: time spent issuing RPC calls");
}


distributed_control::~distributed_control() {
  PERMANENT_DESTROY_DIST_EVENT_LOG(eventlog);
  distributed_services->full_barrier();
  logstream(LOG_INFO) << "Shutting down distributed control " << std::endl;
  
  size_t bytessent = bytes_sent();
  for (size_t i = 0;i < senders.size(); ++i) {
    senders[i]->flush();
  }
  
  comm->close();
  
  for (size_t i = 0;i < senders.size(); ++i) {
    delete senders[i];
  }
  size_t bytesreceived = bytes_received();
  for (size_t i = 0;i < receivers.size(); ++i) {
    receivers[i]->shutdown();
    delete receivers[i];
  }
  senders.clear();
  receivers.clear();
  // shutdown function call handlers
  for (size_t i = 0;i < fcallqueue.size(); ++i) fcallqueue[i].stop_blocking();
  fcallhandlers.join();
  logstream(LOG_INFO) << "Bytes Sent: " << bytessent << std::endl;
  logstream(LOG_INFO) << "Calls Sent: " << calls_sent() << std::endl;
  logstream(LOG_INFO) << "Network Sent: " << network_bytes_sent() << std::endl;
  logstream(LOG_INFO) << "Bytes Received: " << bytesreceived << std::endl;
  logstream(LOG_INFO) << "Calls Received: " << calls_received() << std::endl;
  
  delete comm;

}
  

void distributed_control::exec_function_call(procid_t source,
                                            unsigned char packet_type_mask,
                                            const char* data,
                                            const size_t len) {
  BEGIN_TRACEPOINT(dc_call_dispatch);
  // not a POD call
  if ((packet_type_mask & POD_CALL) == 0) {
    // extract the dispatch function
    boost::iostreams::stream<boost::iostreams::array_source> strm(data, len);
    iarchive arc(strm);
    size_t f;
    arc >> f;
    // a regular funcion call
    dc_impl::dispatch_type dispatch = (dc_impl::dispatch_type)f;
    dispatch(*this, source, packet_type_mask, strm);
  }
  else {
    dc_impl::dispatch_type2 dispatch2 = *reinterpret_cast<const dc_impl::dispatch_type2*>(data);
    dispatch2(*this, source, packet_type_mask, data, len);
  }
  if ((packet_type_mask & CONTROL_PACKET) == 0) inc_calls_received(source);
  END_TRACEPOINT(dc_call_dispatch);
}

void distributed_control::deferred_function_call_chunk(char* buf, size_t len, procid_t src) {
  BEGIN_TRACEPOINT(dc_receive_queuing);
  fcallqueue_entry* fc = new fcallqueue_entry;
  fc->chunk_src = buf;
  fc->chunk_len = len;
  fc->chunk_ref_counter = NULL;
  fc->is_chunk = true;
  fc->source = src;
  fcallqueue_length.inc();
  fcallqueue[src % fcallqueue.size()].enqueue(fc);
  END_TRACEPOINT(dc_receive_queuing);
}


void distributed_control::process_fcall_block(fcallqueue_entry &fcallblock) {
  if (fcallblock.is_chunk == false) {
    for (size_t i = 0;i < fcallblock.calls.size(); ++i) {
      fcallqueue_length.dec();
      exec_function_call(fcallblock.source, fcallblock.calls[i].packet_mask,
                        fcallblock.calls[i].data, fcallblock.calls[i].len);
    }
    if (fcallblock.chunk_ref_counter != NULL) {
      if (fcallblock.chunk_ref_counter->dec(fcallblock.calls.size()) == 0) {
        delete fcallblock.chunk_ref_counter;
        free(fcallblock.chunk_src);
      }
    }
  }
  else {
    fcallqueue_length.dec();
    BEGIN_TRACEPOINT(dc_receive_multiplexing);
    fcallqueue_entry* queuebufs[fcallqueue.size()];
    atomic<size_t>* refctr = new atomic<size_t>(0);
    
    for (size_t i = 0;i < fcallqueue.size(); ++i) {
      queuebufs[i] = new fcallqueue_entry;
      queuebufs[i]->chunk_src = fcallblock.chunk_src;
      queuebufs[i]->chunk_ref_counter = refctr;
      queuebufs[i]->chunk_len = 0;
      queuebufs[i]->source = fcallblock.source;
      queuebufs[i]->is_chunk = false;
    }
    
    //parse the data in fcallblock.data
    char* data = fcallblock.chunk_src;
    size_t remaininglen = fcallblock.chunk_len;
    PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, BYTES_EVENT, remaininglen);
    size_t stripe = 0;
    while(remaininglen > 0) {
      ASSERT_GE(remaininglen, sizeof(dc_impl::packet_hdr));
      dc_impl::packet_hdr hdr = *reinterpret_cast<dc_impl::packet_hdr*>(data);
      ASSERT_LE(hdr.len, remaininglen);
      
      if ((hdr.packet_type_mask & CONTROL_PACKET) == 0) {
        global_bytes_received[hdr.src].inc(hdr.len);
      }
      refctr->value++;
      if (hdr.sequentialization_key == 0) {
        queuebufs[stripe]->calls.push_back(function_call_block(
                                            data + sizeof(dc_impl::packet_hdr), 
                                            hdr.len,
                                            hdr.packet_type_mask));
        ++stripe;
        if (stripe == (fcallblock.source % fcallqueue.size())) ++stripe;
        if (stripe >= fcallqueue.size()) stripe -= fcallqueue.size();
      }
      else {
        size_t idx = (hdr.sequentialization_key % (fcallqueue.size()));
        queuebufs[idx]->calls.push_back(function_call_block(
                                            data + sizeof(dc_impl::packet_hdr), 
                                            hdr.len,
                                            hdr.packet_type_mask));
      }
      data += sizeof(dc_impl::packet_hdr) + hdr.len;
      remaininglen -= sizeof(dc_impl::packet_hdr) + hdr.len;
    }
    END_TRACEPOINT(dc_receive_multiplexing);
    BEGIN_TRACEPOINT(dc_receive_queuing);
    for (size_t i = 0;i < fcallqueue.size(); ++i) { 
      if (queuebufs[i]->calls.size() > 0) {
        fcallqueue_length.inc(queuebufs[i]->calls.size());
        fcallqueue[i].enqueue(queuebufs[i]);
      }
      else {
        delete queuebufs[i];
      }
    }
    END_TRACEPOINT(dc_receive_queuing);
  }
}

void distributed_control::stop_handler_threads(size_t threadid,
                                                size_t total_threadid) {
  for (size_t i = threadid;i < fcallqueue.size(); i += total_threadid) {
    fcallqueue[i].stop_blocking();
    while (fcall_handler_active[i]) usleep(1);
  }
}

void distributed_control::stop_handler_threads_no_wait(size_t threadid,
                                                       size_t total_threadid) {
  for (size_t i = threadid;i < fcallqueue.size(); i += total_threadid) {
    fcallqueue[i].stop_blocking();
  }
}


void distributed_control::start_handler_threads(size_t threadid,
                                                size_t total_threadid) {
  for (size_t i = threadid;i < fcallqueue.size(); i += total_threadid) fcallqueue[i].start_blocking();
  for (size_t i = threadid;i < fcallqueue.size(); i += total_threadid) {
    fcallhandlers.launch(boost::bind(&distributed_control::fcallhandler_loop,
                                      this, i));
  }
}

void distributed_control::handle_incoming_calls(size_t threadid,
                                                size_t total_threadid) {
  for (size_t i = threadid;i < fcallqueue.size(); i += total_threadid) {
    if (fcallqueue[i].empty_unsafe() == false) {
      std::deque<fcallqueue_entry*> q;
      fcallqueue[i].swap(q);
      while (!q.empty()) {
        fcallqueue_entry* entry;
        entry = q.front();
        q.pop_front();

        process_fcall_block(*entry);
        delete entry;
      }
    }
  }
}

void distributed_control::fcallhandler_loop(size_t id) {
  // pop an element off the queue
//  float t = lowres_time_seconds();
  fcall_handler_active[id].inc();
  while(1) {
    fcallqueue[id].wait_for_data();
    if (fcallqueue[id].is_alive() == false) break;
    std::deque<fcallqueue_entry*> q;
    fcallqueue[id].swap(q);
    while (!q.empty()) {
      fcallqueue_entry* entry;
      entry = q.front();
      q.pop_front();
      
      process_fcall_block(*entry);
      delete entry;
    }
    //  std::cerr << "Handler " << id << " died." << std::endl;
  }
  fcall_handler_active[id].dec();
}
  

std::map<std::string, std::string> 
  distributed_control::parse_options(std::string initstring) {
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
  ASSERT_MSG(machines.size() <= RPC_MAX_N_PROCS, 
             "Number of processes exceeded hard limit of %d", RPC_MAX_N_PROCS);
    
  // initialize thread local storage
  if (dc_impl::thrlocal_resizing_array_key_initialized == false) {
    dc_impl::thrlocal_resizing_array_key_initialized = true;
    int err = pthread_key_create(&dc_impl::thrlocal_resizing_array_key, 
                        dc_impl::thrlocal_destructor);
    ASSERT_EQ(err, 0);
  }
  
  if (dc_impl::thrlocal_sequentialization_key_initialized == false) {
    dc_impl::thrlocal_sequentialization_key_initialized = true;
    int err = pthread_key_create(&dc_impl::thrlocal_sequentialization_key, NULL);
    ASSERT_EQ(err, 0);
  }

  //-------- Initialize the full barrier ---------
  full_barrier_in_effect = false;
  procs_complete.resize(machines.size());
  //-----------------------------------------------

  // initialize the counters
  
  global_calls_sent.resize(machines.size());
  global_calls_received.resize(machines.size());
  global_bytes_received.resize(machines.size());
  fcallqueue.resize(numhandlerthreads);

  
  // parse the initstring
  std::map<std::string,std::string> options = parse_options(initstring);

  if (commtype == TCP_COMM) {
    comm = new dc_impl::dc_tcp_comm();
  }
/*  else if (commtype == SCTP_COMM) {
    #ifdef HAS_SCTP
    comm = new dc_impl::dc_sctp_comm();
    std::cerr << "SCTP Communication layer constructed." << std::endl;
    #else
    logger(LOG_FATAL, "SCTP support was not compiled");
    #endif
  }*/
  else {
    ASSERT_MSG(false, "Unexpected value for comm type");
  }
  // create the receiving objects
  if (comm->capabilities() && dc_impl::COMM_STREAM) {
    for (procid_t i = 0; i < machines.size(); ++i) {
      receivers.push_back(new dc_impl::dc_stream_receive(this, i));
      senders.push_back(new dc_impl::dc_buffered_stream_send2(this, comm, i));
    }
  }
  // create the handler threads
  // store the threads in the threadgroup
  fcall_handler_active.resize(numhandlerthreads);
  
  fcallhandlers.resize(numhandlerthreads);
  for (size_t i = 0;i < numhandlerthreads; ++i) {
    fcallhandlers.launch(boost::bind(&distributed_control::fcallhandler_loop, 
                                      this, i));
  }

  
  // set the local proc values
  localprocid = curmachineid;
  localnumprocs = machines.size();
  
  // construct the services
  distributed_services = new dc_services(*this);
  // start the machines
  comm->init(machines, options, curmachineid, 
              receivers, senders);
  std::cerr << "TCP Communication layer constructed." << std::endl;
  compute_master_ranks();
  
#ifdef USE_EVENT_LOG
    PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, *this, std::cout, 3000, dist_event_log::RATE_BAR);
#else
    PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, *this, std::cout, 3000, dist_event_log::LOG_FILE);
#endif
    PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, CALLS_EVENT, "Total RPC Calls");
    PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, BYTES_EVENT, "Total Bytes Communicated");

}

size_t distributed_control::set_sender_option(std::string opt, size_t value) {
  size_t oldval = 0;
  // we assume that all senders are identical
  for (size_t i = 0;i < senders.size(); ++i) {
    oldval = senders[i]->set_option(opt, value);
  }
  return oldval;
}

dc_services& distributed_control::services() {
  return *distributed_services;
}


void distributed_control::barrier() {
  distributed_services->barrier();
}

void distributed_control::flush() {
  for (procid_t i = 0;i < senders.size(); ++i) {
    senders[i]->flush();
  }
}


void distributed_control::compute_master_ranks() {
  uint32_t local_ip = get_local_ip();
  std::string localipandpath = tostr(local_ip) + get_working_dir();
  std::vector<std::string> ipandpath(numprocs());
  ipandpath[procid()] = localipandpath;
  all_gather(ipandpath);
  for(procid_t i = 0; i < ipandpath.size(); ++i) {
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
  __asm("mfence");
  // begin one pass to set all which are already completed
  for (procid_t i = 0;i < numprocs(); ++i) {
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
  barrier();
//   for (size_t i = 0; i < numprocs(); ++i) {
//     std::cout << "Received " << global_calls_received[i].value << " from " << i << std::endl;
//   }
}


} //namespace graphlab

