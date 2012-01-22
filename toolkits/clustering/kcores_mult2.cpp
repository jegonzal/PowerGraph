/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */



#include <cmath>
#include <cstdio>
#include <limits>
#include <iostream>
#include "graphlab.hpp"
#include "graphlab/graph/graph3.hpp"
#include "graphlab/graph/multigraph.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include "../shared/stats.hpp"

#include <graphlab/macros_def.hpp>
using namespace graphlab;

bool debug = false;
int nodes = 0;
uint * edge_count;
unsigned long long total_edges = 0;

struct vertex_data {
}; 

struct edge_data {
};

typedef graphlab::multigraph<vertex_data, edge_data> multigraph_type;
typedef graphlab::graph3<vertex_data, edge_data> graph_type;
graph_type * reference_graph = NULL;

struct kcore_update :
  public graphlab::iupdate_functor<graph_type, kcore_update> {
  void operator()(icontext_type& context) {
  } 
};

class accumulator :
  public graphlab::iaccumulator<graph_type, kcore_update, accumulator> {
public:
  accumulator(){ }

  void operator()(icontext_type& context) {
   

    vertex_data & vdata = context.vertex_data();
    if (debug)
      logstream(LOG_INFO)<<"Entering node: " << context.vertex_id() << std::endl;

    edge_list_type outedgeid = context.out_edges();
    for(size_t i = 0; i < outedgeid.size(); i++) {
      const edge_type& edge = reference_graph->find(context.vertex_id(), outedgeid[i].target());
      if (edge.offset() != (uint)-1)
          edge_count[edge.offset()]++;
      total_edges++;
    }
  };

  void operator+=(const accumulator& other) { 
  }

  void finalize(iglobal_context_type& context) {
   }
}; // end of  accumulator





int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile;
  std::string format = "matrixmarket";
  int unittest = 0;
  int lineformat = MATRIX_MARKET_4;
  bool gzip = true;
  bool stats = false;
  std::string filter = "";
  int reference = 0;
  int max_graph = 1000;
  std::string list_dir = "/usr2/bickson/daily.sorted/";
  std::string dir_path = "/usr2/bickson/bin.graphs/";
  std::string out_dir = "/usr0/bickson/"; 
 
  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("lineformat", &lineformat, lineformat, "matrix line format");
  
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=TBD");
  clopts.attach_option("nodes", &nodes, nodes, "number of nodes"); 
  clopts.attach_option("gzip", &gzip, gzip, "gzipped input file?");
  clopts.attach_option("stats", &stats, stats, "calculate graph stats and exit");
  clopts.attach_option("filter", & filter, filter, "Filter - parse files starting with prefix");
  clopts.attach_option("references", &reference, reference, "reference - why day to compare to?"); 
  clopts.attach_option("max_graph", &max_graph, max_graph, "maximum number of graphs parsed");
  clopts.attach_option("list_dir", &list_dir, list_dir, "directory with a list of file names to parse");
  clopts.attach_option("dir_path", &dir_path, dir_path, "actual directory where files are found");
  clopts.attach_option("outdir", &out_dir, out_dir, "output dir");
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }


  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab Linear solver library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Gaussian Belief Propagation, "
    << "Jacobi method, Conjugate Gradient" << std::endl;



  // Create a core
  graphlab::core<graph_type, kcore_update> core;
  core.set_options(clopts); // Set the engine options

  nodes = 121408373;
  core.set_scope_type("vertex");

    multigraph_type multigraph;
    multigraph.load(list_dir, dir_path, filter, true);
  logstream(LOG_INFO)<<"Going to load reference graph: " << reference << std::endl;
  graphlab::timer mytimer; mytimer.start();
  multigraph.doload(reference);
  reference_graph = multigraph.graph(0);
  edge_count = new uint[multigraph.graph(0)->num_edges()];

   for (int i=0; i< std::min(multigraph.num_graphs(),max_graph); i++){
       if (i != reference){
       accumulator acum;
       multigraph.doload(i);
       core.graph() = *multigraph.graph(1);
       assert(multigraph.get_node_vdata()->size() == (uint)nodes);
       core.add_sync("sync", acum, 1000);
       core.add_global("NUM_ACTIVE", int(0));
       core.sync_now("sync");
       logstream(LOG_INFO)<<mytimer.current_time()<<") Finished giong over graph number " << i << std::endl;
       multigraph.unload(1);
       } 
     }

  
  std::cout << "finished in " << mytimer.current_time() << std::endl;
  for (int i=0; i< 1000; i++)
     std::cout<<i<<": "<<edge_count[i]<<std::endl;

    write_output_vector_binary(out_dir + boost::lexical_cast<std::string>(reference) + "edge_count.bin", edge_count, reference_graph->num_edges());
    return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
