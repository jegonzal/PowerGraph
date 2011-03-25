#ifndef GRAPHLAB_DC_INIT_FROM_MPI_HPP
#define GRAPHLAB_DC_INIT_FROM_MPI_HPP
#include <graphlab/rpc/dc.hpp>
namespace graphlab {
  /** initializes parameters from MPI. Returns true on success
      MPI must be initialized before calling this funciton */
  bool init_param_from_mpi(dc_init_param& param, dc_comm_type commtype = TCP_COMM);
}

#endif // GRAPHLAB_DC_INIT_FROM_MPI_HPP