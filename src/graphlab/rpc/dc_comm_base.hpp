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


#ifndef DC_COMM_BASE_HPP
#define DC_COMM_BASE_HPP
#include <vector>
#include <string>
#include <map>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_receive.hpp>
namespace graphlab {
namespace dc_impl {  

  
/**
The base class of all comms implementations
*/
class dc_comm_base {
 public:
   
  dc_comm_base();
  
  virtual size_t capabilities() const  = 0;
  /**
   Parses initialization parameters. Most of these parameters are
   user provided, or provided on a higher level initialization system.
   It is entirely up to the comm implementation how these parameters to be treated.
   The descriptions here are largely prescriptive.
   All machines are called with the same initialization parameters (of course with the 
   exception of curmachineid)

   The expected behavior is that 
   this fuction should pause until all communication has been set up
   and returns the number of systems in the network.
   After which, all other remaining public functions (numprocs(), send(), etc)
   should operate normally. Every received message should immediate trigger the 
   attached receiver
   
   machines: a vector of string over machine IDs. This is typically provided by the user
             or through some other initialization mechanism
   initstring: Additional parameters passed by the user
   curmachineid: The ID of the current machine. Will be size_t(-1) if this is not available.
                 (Some comm protocols will negotiate this itself.)
   
   receiver: the receiving object
  */
  virtual void init(const std::vector<std::string> &machines,
            const std::map<std::string,std::string> &initopts,
            procid_t curmachineid,
            std::vector<dc_receive*> receiver) = 0;

  /// Must close all connections when this function is called
  virtual void close() = 0;
  
  virtual ~dc_comm_base() {}
  virtual procid_t numprocs() const = 0;
  
  virtual procid_t procid() const = 0;
  
  virtual size_t network_bytes_sent() const = 0;
  virtual size_t network_bytes_received() const = 0;

  
  /** returns true if the channel to the target
  machine is truly open. The dc_comm_base specification allows
  for lazy channels which are not created until it is used.
  For such implementations, this function should return true
  if the channel has been created, and false otherwise. Non-lazy
  implementations should return true all the time.
  The invariant to ensure is that this function must return true
  for a target machine ID if a packet has been sent from this machine
  to the target before this call.
  */
  virtual bool channel_active(size_t target) const = 0;
  
  /**
   Sends the string of length len to the target machine dest.
   Only valid after call to init();
   Establishes a connection if necessary
  */
  virtual void send(size_t target, const char* buf, size_t len) = 0;
 
  virtual void send2(size_t target, 
             const char* buf1, const size_t len1,
             const char* buf2, const size_t len2) = 0; 

  // not required and not used
  virtual void flush(size_t target) = 0;
};

} // namespace dc_impl
} // namespace graphlab
#endif

