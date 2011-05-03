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

#include <cstdio>
#include <cstdlib>
#include <string>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/logger/logger.hpp>
namespace graphlab {

bool init_param_from_env(dc_init_param& param) {
  char* nodeid = getenv("SPAWNID");
  if (nodeid == NULL) {
    return false;
  }
  param.curmachineid = atoi(nodeid);

  char* nodes = getenv("SPAWNNODES");
  std::string nodesstr = nodes;
  if (nodes == NULL) {
    return false;
  }

  param.machines = strsplit(nodesstr, ",");
  for (size_t i = 0;i < param.machines.size(); ++i) {
    param.machines[i] = param.machines[i] + ":" + tostr(10000 + i);
  }
  // set defaults
  param.numhandlerthreads = DEFAULT_NUMHANDLERTHREADS;
  param.commtype = DEFAULT_COMMTYPE;
  return true;
}

} // namespace graphlab

