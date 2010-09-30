#ifndef DISTRIBUTED_CONTROL_TYPES_HPP
#define DISTRIBUTED_CONTROL_TYPES_HPP
#include <inttypes.h>
namespace graphlab {
  typedef uint16_t procid_t;
  
  enum dc_comm_type {
    TCP_COMM, SCTP_COMM
  };
};
#include <graphlab/rpc/dc_packet_mask.hpp>
#endif
