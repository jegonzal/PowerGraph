

#include <iostream>



#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>


#include <distributed_graphlab.hpp>

using namespace graphlab;

typedef double vertex_data;
typedef double edge_data;

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

size_t index(const size_t nrows, const size_t i, const size_t j) {
  return i*nrows + j;
}


int main(int argc, char** argv) {
  std::cout << "Testing graph ingress." << std::endl;

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
 
  std::cout << dc.procid() << ": Starting." << std::endl;
  graph_type graph(dc);
  const size_t nrows = 1000;
  size_t eid = 0;
  for(size_t i = 0; i < nrows; ++i) {
    for(size_t j = 0; j < nrows; ++j) {
      if(eid++ % dc.numprocs() == dc.procid()) {
        const size_t source = index(nrows, i, j);
        if(i+1 < nrows) {
          const size_t target = index(nrows, i+1, j);
          graph.add_edge(source, target);
          graph.add_edge(target, source);
        }
        if(j+1 < nrows) {
          const size_t target = index(nrows, i, j+1);
          graph.add_edge(source, target);
          graph.add_edge(target, source);
        }      
      } // end of if(eid ...
    }
  } // end of for(i ...

  std::cout << dc.procid() << ": Enter Finalize" << std::endl;
  graph.finalize();
  std::cout << dc.procid() << ": Finished" << std::endl;
  std::cout << "Finished!" << std::endl;
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main
