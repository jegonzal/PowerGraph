#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <fstream>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/util/command_line_options.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>
using namespace graphlab;
  
typedef distributed_graph<vertex_color_type, char> graph_type;

enum color_state {
  NOT_DONE,
  VALID_BUT_NOT_MINIMAL,
  MINIMAL
};

dense_bitset db(4096);

color_state color_owned_subgraph(graph_type &g) {
  color_state c = MINIMAL;
  
  foreach(vertex_id_t ownedv, g.owned_vertices()) {
    db.fill();
    foreach(vertex_id_t v, g.in_vertices(ownedv)) {
      db.clear_bit_unsync(g.get_vertex_data(v));
    }
    foreach(vertex_id_t v, g.out_vertices(ownedv)) {
      db.clear_bit_unsync(g.get_vertex_data(v));
    }
    vertex_color_type curcolor = g.get_vertex_data(ownedv);
    uint32_t sc;
    assert(db.first_bit(sc));
    vertex_color_type smallestcolor = (vertex_color_type)(sc);
    
    if (smallestcolor != curcolor) {
      if (c == MINIMAL && db.get(curcolor)) {
        c = VALID_BUT_NOT_MINIMAL;
      }
      else {
        c = NOT_DONE;
      }
      g.set_vertex_data_async(ownedv, smallestcolor);      
    }
  }
  return c;
}


bool check_is_done(std::vector<size_t> vcal, color_state donelevel) {
  for (size_t i = 0;i < vcal.size(); ++i) {
    color_state c = (color_state)(vcal[i]);
    if (c == NOT_DONE) return false;
    else if (c == VALID_BUT_NOT_MINIMAL && donelevel == MINIMAL) return false;
  }
  return true;
}

int main(int argc, char** argv) {

  
  global_logger().set_log_level(LOG_DEBUG);
  
  // Parse the command lines
  std::string path;
  std::string donestring;
  std::string outputfile = "colors.txt";
  graphlab::command_line_options clopts("Graph Colorizer", true);
  size_t roundspersync = 1; 
  clopts.attach_option("atomindex", &path, path, "atom index file");
  clopts.attach_option("outputfile", &outputfile, outputfile, "color output file");
  clopts.attach_option("donelevel", &donestring, donestring, "\'minimal\' or \'valid\'");
  clopts.attach_option("roundspersync", &roundspersync, roundspersync, "rounds per synchronize");
  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }
  color_state donelevel = MINIMAL;

  if (donestring == "minimal") {
    donelevel = MINIMAL;
  }
  else if (donestring == "valid") {
    donelevel = VALID_BUT_NOT_MINIMAL;
  }
  else {
    logstream(LOG_FATAL) << "Invalid done level: " << donestring << std::endl;
    return 1;
  }

  mpi_tools::init(argc, argv);
  dc_init_param param;
  param.initstring="buffered_send=yes";
  // if not running in DC environment, make atoms
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  distributed_control dc(param);


  graph_type dg(dc, path, true);
  size_t iter = 0;
  while(1) {
    iter++;
    for (size_t j = 0; j < roundspersync - 1; ++j) {
      color_state c = color_owned_subgraph(dg);
    }
    dc.full_barrier();
    color_state c = color_owned_subgraph(dg);
    size_t cval = (size_t)c;
    std::vector<size_t> vcal(dc.numprocs(), 0);
    vcal[dc.procid()] = cval;
    dc.all_gather(vcal);
    // check if done
    bool isdone = check_is_done(vcal, donelevel);
    if (isdone) break;
    if (dc.procid() == 0) {
      logstream(LOG_INFO) << "Iteration " << iter << std::endl;
    }
    //dg.push_all_owned_vertices_to_replicas(true);
    dg.rmi.full_barrier();
    if (dc.procid() == 0) {
      logstream(LOG_INFO) << "Synchronized " << iter << std::endl;
    }
    dc.barrier();
  }
  typedef std::vector<std::pair<vertex_id_t, vertex_color_type> > vector_of_pairs;
  std::vector<vector_of_pairs> vertex_and_colors(dc.numprocs());
  foreach(vertex_id_t ownedv, dg.owned_vertices()) {
    vertex_and_colors[dc.procid()].push_back(std::make_pair(ownedv, dg.get_vertex_data(ownedv)));
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
  dc.barrier();
  mpi_tools::finalize();
}
