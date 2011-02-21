#include <iostream>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/util/timer.hpp>
using namespace graphlab;


int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);
  global_logger().set_log_level(LOG_INFO);

  dc_init_param param;
  if (init_param_from_mpi(param) == false) {
    return 0;
  }

  global_logger().set_log_level(LOG_DEBUG);
  distributed_control dc(param);
  timer ti;
  ti.start();
  for (size_t i = 0;i < 100000; ++i) {
    dc.full_barrier();
  }
  if (dc.procid() == 0) {
    std::cout << "100K full barriers in: " << ti.current_time() << std::endl;
  }

  ti.start();
  for (size_t i = 0;i < 100000; ++i) {
    dc.barrier();
  }
  if (dc.procid() == 0) {
    std::cout << "100K barriers in: " << ti.current_time() << std::endl;
  }



  ti.start();
  for (size_t i = 0;i < 100000; ++i) {
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (dc.procid() == 0) {
    std::cout << "100K mpi barriers in: " << ti.current_time() << std::endl;
  }
  mpi_tools::finalize();
}
