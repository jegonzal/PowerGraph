#include <graphlab.hpp>
int main(int argc, char** argv) {
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  std::cout << "Hello from " << dc.procid() << "!\n";
  dc.cout() << "Hello World!\n";

  graphlab::mpi_tools::finalize();
}

