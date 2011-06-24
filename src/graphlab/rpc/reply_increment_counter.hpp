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


#ifndef REPLY_INCREMENT_COUNTER_HPP
#define REPLY_INCREMENT_COUNTER_HPP
#include <string>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

namespace graphlab {

class distributed_control;

namespace dc_impl {
/**
\ingroup rpc
A wrapper around a char array. This structure 
is incapable of freeing itself and must be managed externally
*/
struct blob {
  /// Constructs a blob containing a pointer to a character array with length len
  blob(char* c, size_t len):c(c),len(len) { };
  blob():c(NULL), len(0){ };
  
  char *c;  ///< stored pointer 
  size_t len; ///< stored length
  
  
  /// serialize the char array
  void save(oarchive& oarc) const {
    oarc << len;
    if (len > 0) serialize(oarc, c, len);
  }
  
  /// deserializes a char array. If there is already a char array here, it will be freed
 void load(iarchive& iarc) {
    if (c) ::free(c);
    c = NULL;
    iarc >> len;
    if (len > 0) {
      c = (char*) malloc(len);
      deserialize(iarc, c, len);
    }
  }
  
  /// Free the stored char array.
  void free() {
    if (c) {
      ::free(c);
      c = NULL;
      len = 0;
    }
  }
};

/**
Defines a really useful function that performs an atomic
increment of a flag when called. This is useful for waiting
for a reply to a request
\note: usemutex = false probably does not work and should be deprecated.
\ingroup rpc
\see reply_increment_counter
*/
struct reply_ret_type{
  atomic<size_t> flag;
  blob val;
  bool usemutex;
  mutex mut;
  conditional cond;
  /**
   * Constructs a reply object which waits for 'retcount' replies.
   * usemutex should always be true
   */
  reply_ret_type(bool usemutex, size_t retcount = 1):flag(retcount), 
                                                     usemutex(true) { 
  }
  
  ~reply_ret_type() { }

  /**
   * Waits for all replies to complete. It is up to the 
   * reply implementation to decrement the counter.
   */
  inline void wait() {
    if (usemutex) {
      mut.lock();
      while(flag.value != 0) cond.wait(mut);
      mut.unlock();
    }
    else {
      while(flag.value != 0) sched_yield();
    }
  }
};



/**
 * Like reply_ret_type but can store a blob for each reply. 
 * \ingroup rpc
 * \see stored_increment_counter
 */
struct stored_ret_type{
  atomic<size_t> flag;
  std::map<procid_t, blob> val;
  mutex mut;
  conditional cond;
  /**
   * Constructs a reply object which waits for 'retcount' replies.
   * usemutex should always be true
   */
  stored_ret_type(size_t retcount = 1):flag(retcount) { 
  }
  
  ~stored_ret_type() { }

  /**
   * Waits for all replies to complete. It is up to the 
   * reply implementation to decrement the counter.
   */
  inline void wait() {
      mut.lock();
      while(flag.value != 0) cond.wait(mut);
      mut.unlock();
  }
};

}


/**
 * \ingroup rpc
 * A simple RPC call which converts ptr to a pointer to a reply_ret_type,
 * stores the blob in it, and decrements its reply counter.
 * \see reply_ret_type
 */
void reply_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, dc_impl::blob ret);

/**
 * \ingroup rpc
 * A simple RPC call which converts ptr to a pointer to a stored_ret_type,
 * stores the blob in it, and decrements its reply counter.
 * \see stored_ret_type
 */
void stored_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, dc_impl::blob ret);

}

#endif

