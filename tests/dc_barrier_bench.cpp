#include <iostream>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/util/timer.hpp>
using namespace graphlab;


int main(int argc, char ** argv) {
  dc_init_param param;

  // if not running in DC environment, make atoms
  if (init_param_from_env(param) == false) {
    std::cout << "Must run in RPC environment\n";
    return 0;
  }

  global_logger().set_log_level(LOG_DEBUG);
  distributed_control dc(param);
  timer ti;
  ti.start();
  std::vector<std::vector<size_t> > all_calls_sent(dc.numprocs());


  for (size_t i = 0;i < 100000; ++i) {
    all_calls_sent[dc.procid()].resize(dc.numprocs(), rand());
    dc.all_gather(all_calls_sent, true);    
  }
  if (dc.procid() == 0) {
    std::cout << "100K barriers in: " << ti.current_time() << std::endl;
  }
}
