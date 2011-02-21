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

  std::vector<std::string> data(dc.numprocs());
  
  std::string tempstring;
  tempstring.resize(10*1024*1024,'a');
  data[dc.procid()] = tempstring;

  // test gather
  std::vector<std::string> test1 = data;
  dc.gather(test1, 0);
  if (dc.procid() == 0) {
    for (size_t i = 0;i < test1.size(); ++i) {
      assert(data[dc.procid()] == tempstring);
    }
  }

  // test all gather
  test1 = data;
  dc.all_gather(test1);
  for (size_t i = 0;i < test1.size(); ++i) {
    assert(data[dc.procid()] == tempstring);
  }

  std::vector<size_t> numbers(dc.numprocs(), dc.procid());
  dc.all_gather(numbers);
  for (size_t i = 0;i < dc.numprocs(); ++i) {
    ASSERT_EQ(numbers[i], i);
  }

  dc.barrier();
}
