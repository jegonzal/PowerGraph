#include <string>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>

namespace graphlab {

void reply_increment_counter(distributed_control &dc, procid_t src, 
                             size_t ptr, blob ret) {
  reply_ret_type *a = reinterpret_cast<reply_ret_type*>(ptr);
  a->val=ret;
  a->flag.inc();  
  if (a->usesem) a->sem.post();
}

}
