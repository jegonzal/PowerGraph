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


/**
 * \file
 *
 * Written by Danny Bickson, CMU
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab.hpp>
#include <graphlab/engine/gl3engine.hpp>
#include <graphlab/macros_def.hpp>

#define BOND_PERCULATION_MAP_REDUCE 0
#define BOND_PERCULATION_TRANSFORM 1
static bool debug;
int max_iter = 100000;
std::string predictions;

struct vertex_data : public graphlab::IS_POD_TYPE{
  unsigned int comp_id;
}; 

std::size_t hash_value(vertex_data const& b) {
  return b.comp_id;
}


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data sgdo stores the most recent error estimate.
 */
struct edge_data : public graphlab::IS_POD_TYPE {
   unsigned int id;
   unsigned int comp_id;
   edge_data(unsigned int id) : id(id) { comp_id = id; };
   edge_data(){ id = comp_id = -1; }

}; // end of edge data

std::size_t hash_value(edge_data const& b) {
  return b.comp_id;
}


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::gl3engine<graph_type> engine_type;



/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph,
                         const std::string& filename,
                         const std::string& line) {
  ASSERT_FALSE(line.empty());
  // Parse the line
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  unsigned int edge_id(-1);
  strm >> source_id >> target_id >> edge_id;

  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(edge_id));
  return true; // successful load
} // end of graph_loader

size_t count_component(const graph_type::edge_type & edge) {
  int diff = (edge.source().data().comp_id != edge.target().data().comp_id);
  return diff;
} 

unsigned int bond_perculation_map(const graph_type::vertex_type& center,
                         graph_type::edge_type& edge,
                         const graph_type::vertex_type& other) {
   return edge.data().id;
}

void store_min_component(const graph_type::vertex_type& center,
                         graph_type::edge_type& edge,
                         const graph_type::vertex_type& other) {
   edge.data().comp_id = center.data().comp_id;
}

//find min component of two edges
void bond_perculation_combine(unsigned int& v1, const unsigned int& v2) {
    v1 = std::min(v1, v2);
}

//the main update function
void bond_perculation_function(engine_type::context_type& context,
                  graph_type::vertex_type& vertex) {
       
     vertex.data().comp_id =  context.map_reduce<unsigned int>(BOND_PERCULATION_MAP_REDUCE, graphlab::ALL_EDGES);
     context.edge_transform(BOND_PERCULATION_TRANSFORM, graphlab::ALL_EDGES);

}


struct model_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  
  std::string save_vertex(const vertex_type& vertex) const {
    return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return boost::lexical_cast<std::string>(edge.data().id) + " " + boost::lexical_cast<std::string>(edge.data().comp_id);
  }
}; 




int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description =
      "Compute connected component - by edges";
  graphlab::command_line_options clopts(description);
  std::string input_dir;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("predictions", predictions,
                       "The prefix (folder and filename) to save predictions.");
  clopts.attach_option("max_iter", max_iter,
                       "number of iterations");
  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer;
  graph_type graph(dc, clopts);
  graph.load(input_dir, graph_loader);
  dc.cout() << "Loading graph. Finished in "
            << timer.current_time() << std::endl;

  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in "
            << timer.current_time() << std::endl;


  dc.cout()
      << "========== Graph statistics on proc " << dc.procid()
      << " ==============="
      << "\n Num vertices: " << graph.num_vertices()
      << "\n Num edges: " << graph.num_edges()
      << "\n Num replica: " << graph.num_replicas()
      << "\n Replica to vertex ratio: "
      << float(graph.num_replicas())/graph.num_vertices()
      << "\n --------------------------------------------"
      << "\n Num local own vertices: " << graph.num_local_own_vertices()
      << "\n Num local vertices: " << graph.num_local_vertices()
      << "\n Replica to own ratio: "
      << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
      << "\n Num local edges: " << graph.num_local_edges()
      //<< "\n Begin edge id: " << graph.global_eid(0)
        << "\n Edge balance ratio: "
        << float(graph.num_local_edges())/graph.num_edges()
        << std::endl;

  dc.cout() << "Creating engine" << std::endl;

  engine_type engine(dc, graph, clopts);
  engine.register_map_reduce(BOND_PERCULATION_MAP_REDUCE, bond_perculation_map, bond_perculation_combine);
  engine.register_edge_transform(BOND_PERCULATION_TRANSFORM, store_min_component);
  for (int i=0; i< max_iter; i++){
     engine.parfor_all_local_vertices(bond_perculation_function);
     engine.wait();
     size_t diff = graph.map_reduce_edges<size_t>(count_component);
     dc.cout() << "iter = " << i << " diff= " << diff << std::endl;
     if (diff == 0)
       break;
  }

  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
            << std::endl
            << "Final Runtime (seconds):   " << runtime;

  // Compute the final training error -----------------------------------------
  dc.cout() << "Final error: " << std::endl;
  // Make predictions ---------------------------------------------------------
  if(!predictions.empty()) {
    std::cout << "Saving predictions" << std::endl;
    const bool gzip_output = false;
    const bool save_vertices = false;
    const bool save_edges = true;
    const size_t threads_per_machine = 2;

    //save the output
    graph.save(predictions, model_saver(),
		gzip_output, save_edges, save_vertices, threads_per_machine);
  }
 

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



