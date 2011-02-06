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