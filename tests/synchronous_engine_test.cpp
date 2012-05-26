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

#include <vector>
#include <algorithm>
#include <iostream>


// #include <cxxtest/TestSuite.h>

#include <graphlab.hpp>



class count_in_neighbors : 
  public graphlab::ivertex_program<int, int, int>,
  public graphlab::IS_POD_TYPE {
public:
  edge_dir_type 
  gather_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::IN_EDGES;
  }
  gather_type 
  gather(icontext_type& context, const vertex_type& vertex, 
         edge_type& edge) const {
    return 1;
  }
  void apply(icontext_type& context, vertex_type& vertex, 
             const gather_type& total) {
    ASSERT_EQ( total, int(vertex.num_in_edges()) );
  }
  edge_dir_type 
  scatter_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
}; // end of count neighbors




int main(int argc, char** argv) {

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);


  graphlab::command_line_options clopts("Test code.");
  
  std::cout << "Creating a powerlaw graph" << std::endl;
  typedef count_in_neighbors::graph_type graph_type;
  graph_type graph(dc, clopts);
  graphlab::graph_ops::load_synthetic_powerlaw(graph, 100);

  std::cout << "Constructing a syncrhonous engine" << std::endl;
  typedef graphlab::synchronous_engine<count_in_neighbors> engine_type;
  engine_type engine(dc, graph, clopts);
  engine.initialize();
  std::cout << "Scheduling all vertices to count their neighbors" << std::endl;
  engine.signal_all();

  std::cout << "Running!" << std::endl;
  engine.start();

  std::cout << "Finished" << std::endl;




  graphlab::mpi_tools::finalize();
} // end of main





