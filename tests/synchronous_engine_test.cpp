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

typedef graphlab::distributed_graph<int,int> graph_type;


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
    context.signal(vertex);
  }
  edge_dir_type 
  scatter_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
}; // end of count neighbors


void test_in_neighbors(graphlab::distributed_control& dc,
                       graphlab::command_line_options& clopts,
                       graph_type& graph) {
  std::cout << "Constructing a syncrhonous engine for in neighbors" << std::endl;
  typedef graphlab::synchronous_engine<count_in_neighbors> engine_type;
  engine_type engine(dc, graph, clopts);
  engine.initialize();
  std::cout << "Scheduling all vertices to count their neighbors" << std::endl;
  engine.signal_all();
  std::cout << "Running!" << std::endl;
  engine.start();
  std::cout << "Finished" << std::endl;
}


class count_out_neighbors : 
  public graphlab::ivertex_program<int, int, int>,
  public graphlab::IS_POD_TYPE {
public:
  edge_dir_type 
  gather_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::OUT_EDGES;
  }
  gather_type 
  gather(icontext_type& context, const vertex_type& vertex, 
         edge_type& edge) const {
    return 1;
  }
  void apply(icontext_type& context, vertex_type& vertex, 
             const gather_type& total) {
    ASSERT_EQ( total, int(vertex.num_out_edges()) );
    context.signal(vertex);
  }
  edge_dir_type 
  scatter_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
}; // end of count neighbors

void test_out_neighbors(graphlab::distributed_control& dc,
                        graphlab::command_line_options& clopts,
                        graph_type& graph) {
  std::cout << "Constructing a syncrhonous engine for out neighbors" << std::endl;
  typedef graphlab::synchronous_engine<count_out_neighbors> engine_type;
  engine_type engine(dc, graph, clopts);
  engine.initialize();
  std::cout << "Scheduling all vertices to count their neighbors" << std::endl;
  engine.signal_all();
  std::cout << "Running!" << std::endl;
  engine.start();
  std::cout << "Finished" << std::endl;
}


class count_all_neighbors : 
  public graphlab::ivertex_program<int, int, int>,
  public graphlab::IS_POD_TYPE {
public:
  edge_dir_type 
  gather_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }
  gather_type 
  gather(icontext_type& context, const vertex_type& vertex, 
         edge_type& edge) const {
    return 1;
  }
  void apply(icontext_type& context, vertex_type& vertex, 
             const gather_type& total) {
    ASSERT_EQ( total, int(vertex.num_in_edges() + vertex.num_out_edges() ) );
    context.signal(vertex);
  }
  edge_dir_type 
  scatter_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
}; // end of count neighbors

void test_all_neighbors(graphlab::distributed_control& dc,
                        graphlab::command_line_options& clopts,
                        graph_type& graph) {
  std::cout << "Constructing a syncrhonous engine for all neighbors" << std::endl;
  typedef graphlab::synchronous_engine<count_all_neighbors> engine_type;
  engine_type engine(dc, graph, clopts);
  engine.initialize();
  std::cout << "Scheduling all vertices to count their neighbors" << std::endl;
  engine.signal_all();
  std::cout << "Running!" << std::endl;
  engine.start();
  std::cout << "Finished" << std::endl;
}




class basic_messages : 
  public graphlab::ivertex_program<int, int, int, int>,
  public graphlab::IS_POD_TYPE {
  int message_value;
public:

  void recv_message(icontext_type& context, const vertex_type& vertex,
                    const message_type& msg) {
    message_value = msg;
  } 

  edge_dir_type 
  gather_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::IN_EDGES;
  }
 
  gather_type gather(icontext_type& context, const vertex_type& vertex, 
         edge_type& edge) const {
    return 1;
  }

  void apply(icontext_type& context, vertex_type& vertex, 
             const gather_type& total) {
    context.signal(vertex, 0);
    if(message_value < 0) {
      // first iteration has wrong messages
      return;
    }
    ASSERT_EQ(total, message_value);

  }

  edge_dir_type 
  scatter_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::OUT_EDGES;
  }

  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {
    context.signal(edge.target(), 1);
  }

}; // end of test_messages

void test_messages(graphlab::distributed_control& dc,
                   graphlab::command_line_options& clopts,
                   graph_type& graph) {
  std::cout << "Testing messages" << std::endl;
  typedef graphlab::synchronous_engine<basic_messages> engine_type;
  engine_type engine(dc, graph, clopts);
  engine.initialize();
  std::cout << "Scheduling all vertices to test messages" << std::endl;
  engine.signal_all(-1);
  std::cout << "Running!" << std::endl;
  engine.start();
  std::cout << "Finished" << std::endl;
}







int main(int argc, char** argv) {

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);

  graphlab::command_line_options clopts("Test code.");
  clopts.engine_args.add_option("max_iterations", 10);
  std::cout << "Creating a powerlaw graph" << std::endl;
  graph_type graph(dc, clopts);
  graphlab::graph_ops::load_synthetic_powerlaw(graph, 100);

  test_in_neighbors(dc, clopts, graph);
  test_out_neighbors(dc, clopts, graph);
  test_all_neighbors(dc, clopts, graph);
  test_messages(dc, clopts, graph);

  graphlab::mpi_tools::finalize();
} // end of main





