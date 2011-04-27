#ifndef DISTRIBUTED_CONTROL_TYPES_HPP
#define DISTRIBUTED_CONTROL_TYPES_HPP
#include <inttypes.h>
namespace graphlab {
  /// The type used for numbering processors \ingroup rpc
  typedef uint16_t procid_t;
  
  /**
   * \ingroup rpc
   * The underlying communication protocol
   */
  enum dc_comm_type {
    TCP_COMM,   ///< TCP/IP
    SCTP_COMM   ///< SCTP (limited support)
  };
};
#include <graphlab/rpc/dc_packet_mask.hpp>
#endif
