#ifndef GET_LAST_DC_PROCID_HPP
#define GET_LAST_DC_PROCID_HPP

#include <graphlab/rpc/dc_types.hpp>

namespace graphlab {
namespace dc_impl {
  /**
   * \brief Returns the procid of the current process as set by the latest
   * constructed distributed_control object
   */
  procid_t get_last_dc_procid();
} // dc_impl
} // graphlab
 
#endif // GET_LAST_DC_PROCID_HPP
