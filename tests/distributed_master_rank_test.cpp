#include <string>
#include <algorithm>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>

using namespace graphlab;



int main(int argc, char** argv) {
  dc_init_param param;
  assert(init_param_from_env(param));
  global_logger().set_log_level(LOG_DEBUG);

  distributed_control dc(param);
  if (dc.is_master_rank()) {
    std::cout << "Is Master!\n";
  }
  else {
    std::cout << "Is not Master!\n";
  }
  std::cout << "Master node is " << dc.master_rank() << "\n";
  dc.barrier();
}

