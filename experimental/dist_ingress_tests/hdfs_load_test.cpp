#include <iostream>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>

#define DEBUG_GRAPH

#include <distributed_graphlab.hpp>


#include "hdfs_util.hpp"





#include <graphlab/macros_def.hpp>


typedef double vertex_data;
typedef double edge_data;

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

int main(int argc, char** argv) {
  ///! Initialize control plain using mpi
   graphlab::mpi_tools::init(argc, argv);
   graphlab::dc_init_param rpc_parameters;
   graphlab::init_param_from_mpi(rpc_parameters);
   graphlab::distributed_control dc(rpc_parameters);

  if (argc < 2)
  {
    std::cout << "1 load domain graph\n 2 load web graph" << std::endl;
    exit(0);
  }
  std::string graph_to_load;
  if (atoi(argv[1]) == 0) {
    graph_to_load = "toy domain graph";
  } 
  else if (atoi(argv[1]) == 1) {
    graph_to_load = "domain graph";
  } else if (atoi(argv[1]) == 2) {
    graph_to_load = "web graph";
  } else {
    std::cout << "Invalid argument" << std::endl;
    std::cout << "1 load domain graph\n 2 load web graph" << std::endl;
    exit(0);
  }

  graphlab::timer mytimer; mytimer.start();
  std::cout << "Testing load " << graph_to_load << std::endl;
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);

  std::string edata_path; 

  if (graph_to_load == "toy domain graph")  {
    edata_path = "/users/haijieg/webgraph/domain_graph/edata_splits0/";
  } else if (graph_to_load == "domain graph")  {
    edata_path = "/users/haijieg/webgraph/domain_graph/edata_splits/";
  } else {
    edata_path = "/users/haijieg/webgraph/full_graph/edata_splits/";
  }
  std::string format = "adj";

  std::cout << dc.procid() << ": Starting." << std::endl;
  graph_type graph(dc);

  // Hadoop hdfs init 
  graphlab::hdfs hdfs;
  std::vector<std::string> file_list = hdfs.list_files(edata_path);
  std::sort(file_list.begin(), file_list.end());
  for (size_t i = 0; i < file_list.size(); ++i) {
    const std::string fname = file_list[i];
    if (i % dc.numprocs() == dc.procid()) {
      std::cout << "Proc #" << dc.procid() << " open file " << fname << std::endl;
      const std::string format = "adj";
      const bool success = 
        graphlab::graph_ops::load_structure(fname, format, graph);
      ASSERT_TRUE(success);
    }
  }
  

  std::cout << dc.procid() << ": Enter Finalize" << std::endl;
  graph.finalize();

  std::cout << " ===============================================================" << std::endl;
  std::cout << dc.procid() << ": Finished in " << mytimer.current_time() << std::endl;

  std::cout << "========== Graph statistics on proc " << dc.procid() << " ==============="
    << "\n Num vertices: " << graph.num_vertices()
    << "\n Num edges: " << graph.num_edges()
    << "\n Num replica: " << graph.num_replicas()
    << "\n Replica to vertex ratio: " << (float)graph.num_replicas()/graph.num_vertices()
    << "\n --------------------------------------------" 
    << "\n Num local own vertices: " << graph.num_local_own_vertices()
    << "\n Num local vertices: " << graph.num_local_vertices()
    << "\n Replica to own ratio: " << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
    << "\n Num local edges: " << graph.num_local_edges()
    << "\n Begin edge id: " << graph.global_eid(0)
    << "\n Edge balance ratio: " << (float)graph.num_local_edges()/graph.num_edges()
    << std::endl;



  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main

#include <graphlab/macros_undef.hpp>
