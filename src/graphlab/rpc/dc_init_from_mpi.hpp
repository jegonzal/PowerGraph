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

#ifndef GRAPHLAB_DC_INIT_FROM_MPI_HPP
#define GRAPHLAB_DC_INIT_FROM_MPI_HPP
#include <graphlab/rpc/dc.hpp>
namespace graphlab {
  /**
   * \ingroup rpc 
   * initializes parameters from MPI. Returns true on success
      MPI must be initialized before calling this function */
  bool init_param_from_mpi(dc_init_param& param, dc_comm_type commtype = TCP_COMM);
}

#endif // GRAPHLAB_DC_INIT_FROM_MPI_HPP


