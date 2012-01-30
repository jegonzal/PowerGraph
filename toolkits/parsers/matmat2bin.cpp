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

#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include <graphlab/graph/graph2.hpp>
#include <graphlab/graph/graph3.hpp>
#include <graphlab/graph/graph.hpp>

using namespace graphlab;
using namespace std;

bool debug;

struct vertex_data{ 
 void add_self_edge(double val){};

};
struct edge_data{ 
 double val;
 edge_data(double _val){ val = _val; };
};

typedef graphlab::graph2<vertex_data, edge_data> graph_type;
int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile, outdir;
  std::string format = "matrixmarket";
  size_t sync_interval = 10000;
  int unittest = 0;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix input file");
  clopts.add_positional("data");
  clopts.attach_option("outdir", &outdir, outdir, "Output directory");
 
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab V2 matrix factorization library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Lanczos" << std::endl;


  graph_type graph2;
  std::cout << "Load matrix " << datafile << std::endl;
  bipartite_graph_descriptor info;
  load_graph(datafile, format, info, graph2);
  graph2.finalize();
  save_to_bin(outdir + datafile, graph2, true);
}
