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


#ifndef DC_INTERNAL_TYPES_HPP
#define DC_INTERNAL_TYPES_HPP
#include <boost/function.hpp>
#include <boost/unordered_map.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/util/resizing_array_sink.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
namespace graphlab {
class distributed_control;

namespace dc_impl {  

/** 
 * \ingroup rpc_internal
 * 
 * The type of the callback function used by the communications
classes when data is received*/
typedef void (*comm_recv_callback_type)(void* tag, procid_t src, 
                                        const char* buf, size_t len);

/** 
 * \ingroup rpc_internal
 * 
 * The type of the local function call dispatcher */
typedef void (*dispatch_type)(distributed_control& dc, procid_t, unsigned char, std::istream&);

typedef boost::unordered_map<std::string, dispatch_type> dispatch_map_type;

// commm capabilities
const size_t COMM_STREAM = 1;
const size_t COMM_DATAGRAM = 0;

/** 
 * \ingroup rpc_internal
 * 
 * The header form of each packet */
struct packet_hdr {
  uint64_t len; /// length of the packet
  procid_t src; /// source machine
  unsigned char packet_type_mask; /// the types are in dc_packet_mask.hpp
  unsigned char sequentialization_key;
};

/** 
 * \ingroup rpc_internal
 * 
 * special handling for the only pointer datatype 
we natively support serialization for. Basically,
we must delete it. if charstring_free is called on a
char*, it will be deleted. Otherwise it will not do anything*/
template <typename T> 
inline void charstring_free(T& t) { }

template <>
inline void charstring_free<char*>(char* &c){
    delete [] c;
};

// if 1, uses semaphore
// if 0 uses spin wait
#define REQUEST_WAIT_METHOD 1


/** 
 * \ingroup rpc_internal
 * 
 * 
 * The data needed to receive the matched send / recvs */
struct recv_from_struct {
  inline recv_from_struct():tag(0), hasdata(false) { }
  
  std::string data;
  size_t tag;
  mutex lock;
  conditional cond;
  bool hasdata;
  
};

/**
 * \ingroup rpc_internal
 * 
 * Used for termination detection
 */
struct terminator_token {
  terminator_token():calls_sent(0),calls_recv(0),terminate(false) { }
  terminator_token(size_t sent, size_t recv):calls_sent(sent),
                          calls_recv(recv),terminate(false) { }
  size_t calls_sent;
  size_t calls_recv;
  bool terminate;
};

extern boost::iostreams::stream<resizing_array_sink_ref>& get_thread_local_stream();

}
}

SERIALIZABLE_POD(graphlab::dc_impl::terminator_token);
#endif

