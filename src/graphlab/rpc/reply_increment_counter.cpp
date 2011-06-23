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


#include <string>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>

namespace graphlab {

void reply_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, dc_impl::blob ret) {
  dc_impl::reply_ret_type *a = reinterpret_cast<dc_impl::reply_ret_type*>(ptr);
  a->mut.lock();
  a->val = ret;
  size_t retval = a->flag.dec();  
  if (retval == 0 && a->usemutex) {
    a->cond.signal();
  }
  a->mut.unlock();
}

void stored_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, dc_impl::blob ret) {
  dc_impl::stored_ret_type *a = reinterpret_cast<dc_impl::stored_ret_type*>(ptr);
  a->mut.lock();
  a->val[src] = ret;
  size_t retval = a->flag.dec();  
  if (retval == 0) {
    a->cond.signal();
  }
  a->mut.unlock();
}

}

