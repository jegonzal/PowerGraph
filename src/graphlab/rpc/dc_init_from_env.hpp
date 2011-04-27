#ifndef GRAPHLAB_DC_INIT_FROM_ENV_HPP
#define GRAPHLAB_DC_INIT_FROM_ENV_HPP
#include <graphlab/rpc/dc.hpp>
namespace graphlab {
  /** 
   * \ingroup rpc
   * initializes parameters from environment. Returns true on success */
  bool init_param_from_env(dc_init_param& param);
}

#endif // GRAPHLAB_DC_INIT_FROM_ENV_HPP



