#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/net_util.hpp>
#include <graphlab/logger/logger.hpp>

#ifdef HASMPI
#include <graphlab/util/mpi_tools.hpp>
#endif
namespace graphlab {

bool init_param_from_mpi(dc_init_param& param,dc_comm_type commtype) {
#ifdef HASMPI
  ASSERT_MSG(commtype == TCP_COMM, "MPI initialization only supports TCP at the moment");
  // Look for a free port to use. 
  std::pair<size_t, int> port_and_sock = get_free_tcp_port();
  size_t port = port_and_sock.first;
  int sock = port_and_sock.second;
  
  std::string ipaddr = get_local_ip_as_str();
  ipaddr = ipaddr + ":" + tostr(port);
  // now do an allgather
  logstream(LOG_INFO) << "Will Listen on: " << ipaddr << std::endl;
  std::vector<std::string> machines;
  mpi_tools::all_gather(ipaddr, param.machines);
  // set defaults
  param.curmachineid = mpi_tools::rank();

  param.numhandlerthreads = DEFAULT_NUMHANDLERTHREADS;
  param.commtype = commtype;
  param.initstring = param.initstring + std::string(" __sockhandle__=") + tostr(sock) + " ";
  return true;
#else
  std::cerr << "MPI Support not compiled!" << std::endl;
  exit(0);
#endif
}

} // namespace graphlab


