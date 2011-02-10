#ifndef GRAPHLAB_RPC_MATCHED_SEND_RECV_HPP
#define GRAPHLAB_RPC_MATCHED_SEND_RECV_HPP
#include <string>
#include <graphlab/rpc/dc_types.hpp>

namespace graphlab {

class distributed_control;

namespace dc_impl {

void block_and_wait_for_recv(distributed_control &dc,
                             procid_t src,
                             std::string& str,
                             size_t tag);


}

};
#endif
