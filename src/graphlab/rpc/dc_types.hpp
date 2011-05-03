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
