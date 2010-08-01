#include <graphlab/distributed/distributed_control.hpp>
using namespace graphlab;

int main(int argc, char **argv) {
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(&argc, &argv);
}
