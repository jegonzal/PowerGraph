/*  
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


#include <graphlab.hpp>
#include <math.h>
#include <cstdlib>
#include <ctime>

double infection_chance;
double recovery_chance;

enum Status {INFECTED, SUSCEPTIBLE, RECOVERED};

// The vertex data is its status (S, I, or R) 
typedef Status vertex_data_type;

// infected_status counts the number of infected neighbors, since
// the number of infected neighbors determines how likely a susceptible node is
// to be infected
struct infected_status: public graphlab::IS_POD_TYPE {
  int value;
  vertex_data_type status;

  infected_status() {
    value = 0;
    status = INFECTED;
  }

  infected_status& operator+=(const infected_status& other) {
    if (other.status == INFECTED) {
      this->value++;
    }

    return *this;
  }
};


typedef infected_status gather_type;
 
// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, graphlab::empty> graph_type;

bool line_parser(graph_type& graph, const std::string& filename, const std::string& textline) {
  std::stringstream strm(textline);
  graphlab::vertex_id_type vid;

  char label;
  // first entry in the line is a vertex ID
  strm >> vid;

  // next entry is their status (S, I, or R)
  strm >> label;

  vertex_data_type statusLabel;
  if (label == 'S') {
    statusLabel = SUSCEPTIBLE;
  } else if (label == 'I') {
    statusLabel = INFECTED;
  } else {
    statusLabel = RECOVERED;
  }

  
  // insert this vertex with its label 
  graph.add_vertex(vid, statusLabel);

  // while there are elements in the line, continue to read until we fail
  while(1) {
    graphlab::vertex_id_type other_vid;
    strm >> other_vid;
    if (strm.fail()) {
      break;
    }
    graph.add_edge(vid, other_vid);
  }

  return true;
}

class cascades:
  public graphlab::ivertex_program<graph_type, gather_type>,
  public graphlab::IS_POD_TYPE {

  public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::ALL_EDGES;
    }

    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      // figure out which data to get from the edge.
      bool isEdgeSource = (vertex.id() == edge.source().id());
      vertex_data_type neighbor_status = isEdgeSource ? edge.target().data() : edge.source().data();

      // create infected_status and add neighbor's status to it. 

      infected_status status;
      status.status = neighbor_status;
      status.value = 1;
      
      return status;
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
      vertex_data_type old_data = vertex.data();

      vertex_data_type result = old_data;
      double random_value;

      // if vertex.data == RECOVERED, don't do anything
      // if vertex.data == INFECTED, roll on recovery_chance to see if recovery
      // occurs 
      // if vertex.data == SUSCEPTIBLE, then do (total) dice rolls (each time comparing
      // the result to infection_chance). if any of them show up as positive,
      // set vertex.data to INFECTED.  

      if (old_data != RECOVERED) {
        if (old_data == INFECTED) {
          random_value = ((double)rand())/RAND_MAX;
          if (random_value <= recovery_chance) {
            result = RECOVERED;
          }
        } else if (old_data == SUSCEPTIBLE) {
          for (int i = 0; i < total.value; i++) {
            random_value = ((double)rand())/RAND_MAX;
            if (random_value <= infection_chance) {
              result = INFECTED;
              break;
            }
          }
        }
      }

      vertex.data() = result;

      if (result == INFECTED) {
        context.signal(vertex);
      }
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
        return graphlab::NO_EDGES;
    }
  };

struct cascades_writer{
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;

    Status status = v.data();
    char vertex_data;
    
    // Convert the status back into a char
    if (status == INFECTED) {
      vertex_data = 'I';
    } else if (status == SUSCEPTIBLE) {
      vertex_data = 'S';
    } else {
      vertex_data = 'R';
    }
    strm << v.id() << "\t" << vertex_data << "\n";
    return strm.str();
  }
  std::string save_edge (graph_type::edge_type e) { return ""; }
};


int main(int argc, char** argv) {
  srand((unsigned) time(0));
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);
  
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Label Propagation algorithm.");
  std::string graph_dir;
  std::string execution_type = "synchronous";
  double recovery = -1;
  double infection = -1;
  size_t iterations = -1;

  clopts.attach_option("graph", graph_dir, "The graph file. Required ");
  clopts.add_positional("graph");
  clopts.attach_option("execution", execution_type, "Execution type (synchronous or asynchronous)");

  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  clopts.attach_option("recovery chance", recovery, "Chance of recovery for an infected individual at each step. Required.");
  clopts.attach_option("infection chance", infection, "Chance of infection for a susceptible individual per person at each step. Required.");

  clopts.attach_option("iterations", iterations, "If set, will force the use of synchronous engine overriding any engine option set by the --engine parameter. Runs cascades for a fixed number of iterations. Also overrides the max_iterations option in the engine.");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (graph_dir == "") {
    dc.cout() << "Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }

  if (recovery == -1) {
    dc.cout() << "Recovery chance not specified. Cannot continue";
    return EXIT_FAILURE;
  }

  if (infection == -1) {
    dc.cout() << "Infection chance not specified. Cannot continue";
    return EXIT_FAILURE;
  }

  infection_chance = infection;
  recovery_chance = recovery;
  
  if (iterations != -1) {
    // make sure this is the synchronous engine
    dc.cout() << "--iterations set. Forcing Synchronous engine, and running for "
    << iterations  << " iterations." << std::endl;
    clopts.get_engine_args().set_option("type", "synchronous");
    clopts.get_engine_args().set_option("max_iterations", iterations);
  }
 
  // Build the graph ----------------------------------------------------------
  graph_type graph(dc);
  dc.cout() << "Loading graph using line parser" << std::endl;
  graph.load(graph_dir, line_parser);
  // must call finalize before querying the graph
  graph.finalize();

  dc.cout() << "#vertices: " << graph.num_vertices() << " #edges:" << graph.num_edges() << std::endl;

  graphlab::omni_engine<cascades> engine(dc, graph, execution_type, clopts);

  engine.signal_all();
  engine.start();

  const float runtime = engine.elapsed_seconds();
  dc.cout() << "Finished Running engine in " << runtime << " seconds." << std::endl;

  if (saveprefix != "") {
    graph.save(saveprefix, cascades_writer(),
       false,  // do not gzip
       true,   //save vertices
       false); // do not save edges 
  }
  

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
}
