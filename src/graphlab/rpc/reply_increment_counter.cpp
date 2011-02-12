#include <string>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>

namespace graphlab {

void reply_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, dc_impl::blob ret) {
  dc_impl::reply_ret_type *a = reinterpret_cast<dc_impl::reply_ret_type*>(ptr);
  a->val=ret;
  size_t retval = a->flag.dec();  
  if (retval == 0 && a->usemutex) a->cond.signal();
}

}
