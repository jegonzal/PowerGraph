#include <string>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/distributed2/distributed_glshared.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>

using namespace graphlab;

distributed_glshared<size_t> anumber;
distributed_glshared<std::string> astring;

int main(int argc, char** argv) {
  dc_init_param param;
  global_logger().set_log_level(LOG_DEBUG);

  //distributed_control dc(machines,"buffered_send=yes,buffered_recv=yes", machineid, 8, SCTP_COMM);
  distributed_control dc(param);
  distributed_glshared_manager glmanager(dc);
  
}
