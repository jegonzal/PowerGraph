#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/matched_send_recvs.hpp>
#include <graphlab/rpc/dc_dist_object_base.hpp>

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



void obj_block_and_wait_for_recv(distributed_control &dc,
                                 procid_t src,
                                 size_t objid,
                                 std::string& str,
                                 size_t tag) {
  dc_dist_object_base* rmi = dc.get_rmi_instance(objid);
  rmi->recv_froms[src].lock.lock();
  rmi->recv_froms[src].data = str;
  rmi->recv_froms[src].tag = tag;
  rmi->recv_froms[src].hasdata = true;
  rmi->recv_froms[src].cond.signal();
  rmi->recv_froms[src].lock.unlock();
}



}
}
