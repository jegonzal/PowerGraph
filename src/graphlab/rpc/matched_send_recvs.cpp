#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/matched_send_recvs.hpp>

namespace graphlab {
namespace dc_impl {




void block_and_wait_for_recv(distributed_control &dc,
                             procid_t src,
                             std::string& str,
                             size_t tag) {
  dc.recv_froms[src].lock.lock();
  dc.recv_froms[src].data = str;
  dc.recv_froms[src].tag = tag;
  dc.recv_froms[src].hasdata = true;
  dc.recv_froms[src].cond.signal();
  dc.recv_froms[src].lock.unlock();
}



}
}
