#include <stdint.h>
#include <string>

#include <graphlab/distributed/graph/cloned_graph.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
#include <graphlab.hpp>


using namespace graphlab;

void apply_fun(size_t index,
               const ishared_data<blob_cloned_graph>& sdm,
               any& current_data,
               const any& acc) {
  current_data.as<size_t>() += acc.as<size_t>();
}
  
int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);

  distributed_control dc(&argc, &argv);
  
  dc.init_message_processing(2);
  distributed_shared_data<blob_cloned_graph> dsdm(dc);
  std::cout << std::endl;
  dc.barrier();
  if (dc.procid() == 0){
    dsdm.set_constant(1, size_t(1));
    dsdm.set_constant(2, size_t(2));
    dsdm.set_constant(3, size_t(3));
    dsdm.atomic_set(4, size_t(4));
  }
  dc.barrier();
  ASSERT_EQ(dsdm.get_constant(1).as<size_t>(), 1);
  ASSERT_EQ(dsdm.get_constant(2).as<size_t>(), 2);
  ASSERT_EQ(dsdm.get_constant(3).as<size_t>(), 3);
  ASSERT_EQ(dsdm.get(4).as<size_t>(), 4);
  dc.barrier();
  dsdm.atomic_apply(4, apply_fun, size_t(1));
  dc.barrier();
  ASSERT_EQ(dsdm.get(4).as<size_t>(), 4 + dc.numprocs());
  dc.barrier();
}




