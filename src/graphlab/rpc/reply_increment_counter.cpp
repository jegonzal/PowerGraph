/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>

namespace graphlab {

void reply_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, dc_impl::blob ret) {
  dc_impl::reply_ret_type *a = reinterpret_cast<dc_impl::reply_ret_type*>(ptr);
  a->mut.lock();
  a->val=ret;
  size_t retval = a->flag.dec();  
  if (retval == 0 && a->usemutex) {
    a->cond.signal();
  }
  a->mut.unlock();
}

}
