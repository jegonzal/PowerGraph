#include <iostream>
#include <string>
#include <map>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/async_consensus.hpp>
using namespace graphlab;

async_consensus* cons;

void test(distributed_control &dc, procid_t procid, size_t i) {
  if (i > 0) {
    cons->cancel();
    size_t targetproc = (dc.procid() + 1) % dc.numprocs();
    dc.remote_call(targetproc, test, i - 1);
  }
  else {
    std::cout << "Done!" << std::endl;
  }
}



int main(int argc, char ** argv) {
  /** Initialization */
  mpi_tools::init(argc, argv);
  global_logger().set_log_level(LOG_DEBUG);

  dc_init_param param;
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  distributed_control dc(param);

  async_consensus consensus(dc);
  cons = &consensus;
  dc.barrier();
  test(dc, dc.procid(), 5000);
  while(!consensus.done());
  dc.barrier();
  mpi_tools::finalize();
}
