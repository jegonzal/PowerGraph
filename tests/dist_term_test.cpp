#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_terminator.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/util/timer.hpp>
using namespace graphlab;

atomic<size_t> numrecvs;

void pong_handler(distributed_control& dc, size_t source, void* ptr, size_t len) {
  std::cout << "pong" << std::endl;
  numrecvs.inc();
}


void testterm(distributed_control &dc,
                distributed_terminator &term) {
  size_t numtransmitted = 0;
  for (size_t i = 0;i < 10; ++i) {
    dc.remote_call(rand() % dc.numprocs(), pong_handler, NULL, 0);
  }
  numtransmitted = 10;
  while (!term.done(numtransmitted, numrecvs.value)) {
    sched_yield();
  }
}




int main(int argc, char **argv) {
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(1);
  distributed_terminator term(dc);
  numrecvs.value = 0;
  dc.barrier();

  testterm(dc, term);
  dc.barrier();
}
