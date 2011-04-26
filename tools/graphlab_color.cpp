#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>


#include <graphlab.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/schedulers/multiqueue_fifo_scheduler.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/distributed_locking_engine.hpp>
#include <graphlab/distributed2/distributed_chromatic_engine.hpp>
#include <graphlab/distributed2/distributed_glshared.hpp>
using namespace graphlab;


#include <graphlab/macros_def.hpp>

/// GLOBAL CONSTANTS

enum color_state {
  VALID_BUT_NOT_MINIMAL,
  MINIMAL
};

color_state donelevel;

///////////////////////////////////////////////////////////////////////////////
///////////////////////// Types ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
typedef vertex_color_type vertex_data_type;

struct edge_data_type { };
SERIALIZABLE_POD(edge_data_type);


typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;
typedef distributed_locking_engine< graph_type, multiqueue_fifo_scheduler<graph_type> > engine_type;
//typedef distributed_chromatic_engine< graph_type > engine_type;
typedef engine_type::iscope_type iscope_type;
typedef engine_type::icallback_type icallback_type;
typedef engine_type::icallback_type icallback_type;
typedef engine_type::update_task_type update_task_type;















void color_update_function(iscope_type& scope,
                          icallback_type& callback) {
  vertex_color_type curcolor = scope.const_vertex_data();
  dense_bitset db(4096);
  db.fill();
  foreach(edge_id_t ine, scope.in_edge_ids()) {
    db.clear_bit_unsync(scope.const_neighbor_vertex_data(scope.source(ine)));
  }
  foreach(edge_id_t oute, scope.out_edge_ids()) {
    db.clear_bit_unsync(scope.const_neighbor_vertex_data(scope.target(oute)));
  }
  uint32_t sc;
  assert(db.first_bit(sc));
  vertex_color_type smallestcolor = (vertex_color_type)(sc);
  // current color is not the smallest color!
  if (smallestcolor != curcolor) {
    // then if we want a minimal coloring, make this the smallest color
    // or if it is not a valid coloring, make it so
    if (donelevel == MINIMAL || db.get(curcolor) == 0) {
      scope.vertex_data() = smallestcolor;
      
      // schedule neighbors
      foreach(edge_id_t ine, scope.in_edge_ids()) {
        callback.add_task(scope.source(ine), color_update_function);
      }
      foreach(edge_id_t oute, scope.out_edge_ids()) {
        callback.add_task(scope.target(oute), color_update_function);
      }
    }
  }
  
  
}






















int main(int argc, char** argv) {
  // set the global logger
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  // Initialize the mpi tools
  graphlab::mpi_tools::init(argc, argv);

  // Parse the command lines
  std::string path;
  std::string donestring;
  std::string outputfile = "colors.txt";
  graphlab::command_line_options clopts("Graph Colorizer");
  clopts.attach_option("atomindex", &path, path, "atom index file");
  clopts.attach_option("outputfile", &outputfile, outputfile, "color output file");
  clopts.attach_option("donelevel", &donestring, donestring, "\'minimal\' or \'valid\'");
  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }
  donelevel = MINIMAL;

  dc_init_param param;
  if( ! init_param_from_mpi(param) ) 
    logstream(LOG_FATAL) 
      << "Failed MPI laucher!" << std::endl; 
  param.initstring = "buffered_queued_send=yes";
  param.numhandlerthreads = 4;
  distributed_control dc(param);
  
  const bool NO_LOAD_DATA(true);
  graph_type  graph(dc, path, NO_LOAD_DATA);
  

  engine_type engine(dc, graph, clopts.get_ncpus());

  engine.add_task_to_all(color_update_function, 1.0);
  
    // Scheduling tasks
  std::cout << "Running Colorizer." << std::endl;
  engine.start();
  
  std::cout << "Finished" << std::endl;
  
  
  typedef std::vector<std::pair<vertex_id_t, vertex_color_type> > vector_of_pairs;
  std::vector<vector_of_pairs> vertex_and_colors(dc.numprocs());
  foreach(vertex_id_t ownedv, graph.owned_vertices()) {
    vertex_and_colors[dc.procid()].push_back(std::make_pair(ownedv, graph.get_vertex_data(ownedv)));
  }
  dc.gather(vertex_and_colors, 0);

  if (dc.procid() == 0) {
    vector_of_pairs result;
    for (size_t i = 0;i < dc.numprocs(); ++i) {
      std::copy(vertex_and_colors[i].begin(), vertex_and_colors[i].end(),
                std::back_inserter(result));
      vector_of_pairs().swap(vertex_and_colors[i]);
      
    }
    std::sort(result.begin(), result.end()); 
    std::ofstream fout;
    fout.open(outputfile.c_str());
    for (size_t i = 0;i < result.size(); ++i) {
      fout << (unsigned int)(result[i].second) << "\n";
    }
  }
  
  // Gather metrics
  dc.fill_metrics();
  graph.fill_metrics();

  
  if(dc.procid() == 0) {
    basic_reporter reporter;
    dc.report_metrics(reporter);
    graph.report_metrics(reporter);
  } 
  
  logstream(LOG_INFO) << "Finished " << dc.procid() << std::endl;
  dc.full_barrier();
  graphlab::mpi_tools::finalize();  
  return EXIT_SUCCESS;
}

