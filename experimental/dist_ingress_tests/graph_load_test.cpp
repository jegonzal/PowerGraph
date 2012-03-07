#include <iostream>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/../../toolkits/shared/io.hpp>

#define DEBUG_GRAPH

#include <distributed_graphlab.hpp>
#include <graphlab/graph/graph_ops.hpp>
#include <google/malloc_extension.h>
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

  graphlab::command_line_options clopts("Distributed graph load test.");
  clopts.use_distributed_options();

  size_t suffix_len = 3;
  std::string graphpath; 
  bool savebin = false;
  bool loadbin = false;
  std::string binpath= "./"; 
  std::string binprefix= "subgraphx"; 
  std::string format = "adj";

  clopts.attach_option("graph", &graphpath, graphpath,
                       "The graph path \n");

  clopts.attach_option("savebin", &savebin, savebin,
                       "Option to save the graph as binary\n");

  clopts.attach_option("loadbin", &loadbin, loadbin,
                       "Option to load the graph from binary\n");

  clopts.attach_option("binpath", &binpath, binpath,
                       "The path for load/save binary file\n");

  clopts.attach_option("binprefix", &binprefix, binprefix,
                       "The prefix for load/save binary file\n");

  clopts.attach_option("format", &format, format,
                       "format of the graph: {adj, snap}\n");

  clopts.attach_option("suffix_len", &suffix_len, suffix_len,
                       "length of the id suffix of the graph file: {3,4,5...}\n");


  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  graphlab::timer mytimer; mytimer.start();
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graph_type graph(dc, clopts);

  if (loadbin){
    logstream(LOG_INFO) << "Load graph from binary." << std::endl;
    graph.load(binpath, binprefix);
    dc.barrier();
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
  }

  std::cout << dc.procid() << ": Starting." << std::endl;
  std::cout << dc.procid() << ": Load from " << graphpath << std::endl;

  std::vector<std::string> files;
  if(boost::starts_with(graphpath, "hdfs://")) {
      graphlab::hdfs hdfs;
      files = hdfs.list_files(graphpath);
  } else {
      graphlab::fs_util::list_files_with_prefix(graphpath, "x", files);
      for(size_t i = 0; i < files.size(); ++i)
        files[i] = graphpath + files[i];
  }

  std::sort(files.begin(), files.end());
  for(size_t i = 0; i < files.size(); ++i) {
    if (i % dc.numprocs() == dc.procid()) {
       std::cout << "Loading graph from structure file: " << files[i] << std::endl;
       bool success = graphlab::graph_ops::load_structure(files[i], format, graph);
       ASSERT_TRUE(success);
    }
  }

  std::cout << dc.procid() << ": Enter Finalize" << std::endl;
  graph.finalize();

  if (savebin) {
    graph.save(binpath, binprefix);
  }

  std::cout 
    << " ===============================================================" 
    << std::endl;
  std::cout << dc.procid() << ": Finished in " << mytimer.current_time() << std::endl;

  std::cout 
    << "========== Graph statistics on proc " << dc.procid() << " ==============="
    << "\n Num vertices: " << graph.num_vertices()
    << "\n Num edges: " << graph.num_edges()
    << "\n Num replica: " << graph.num_replicas()
    << "\n Replica to vertex ratio: " 
    << (float)graph.num_replicas()/graph.num_vertices()
    << "\n --------------------------------------------" 
    << "\n Num local own vertices: " << graph.num_local_own_vertices()
    << "\n Num local vertices: " << graph.num_local_vertices()
    << "\n Replica to own ratio: " 
    << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
    << "\n Num local edges: " << graph.num_local_edges()
    << "\n Begin edge id: " << graph.global_eid(0)
    << "\n Edge balance ratio: " << (float)graph.num_local_edges()/graph.num_edges()
    << std::endl;

   size_t value;
   MallocExtension::instance()->GetNumericProperty("generic.heap_size", &value);
   std::cout << "Heap Size: " << (double)value/(1024*1024) << "MB" << "\n";
   MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &value);
   std::cout << "Allocated Size: " << (double)value/(1024*1024) << "MB" << "\n";

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main

#include <graphlab/macros_undef.hpp>
