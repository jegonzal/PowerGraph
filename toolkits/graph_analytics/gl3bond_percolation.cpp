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
 * See algorithm description in Wikipedia: http://en.wikipedia.org/wiki/Percolation_theory
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab.hpp>
#include <graphlab/engine/gl3engine.hpp>
#include <graphlab/macros_def.hpp>
#include <map>

#define BOND_PERCOLATION_MAP_REDUCE 0
#define SYNC
bool debug;
int max_iter = 100000;
std::string output_file;
int n = 0;
float p = 0;

struct vertex_data : public graphlab::IS_POD_TYPE{
  unsigned int comp_id;
  vertex_data(): comp_id(-1) {}
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
  if (edge_id == (uint)-1)
    logstream(LOG_FATAL)<<"Input file must contains line with the following format: [from] [ to] [edge_id]\n " << std::endl << " where edge_id is a consecutive integer " << std::endl;

  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(edge_id));
  return true; // successful load
} // end of graph_loader

size_t count_component(const graph_type::edge_type & edge) {
  int diff = (edge.source().data().comp_id != edge.target().data().comp_id);
  if (debug && diff)
    std::cout<<"Adding diff between node: " << edge.source().id() << " to: " << edge.target().id()<< " compA: " << edge.source().data().comp_id << " compB: " <<
      edge.target().data().comp_id << std::endl;
  return diff;
}

unsigned int bond_percolation_map(const graph_type::vertex_type& center,
                         graph_type::edge_type& edge,
                         const graph_type::vertex_type& other) {
   if (debug)
     std::cout<<"Comparing: " << center.data().comp_id << " : " << edge.data().id << " : " << other.data().comp_id << std::endl;
   edge.data().comp_id =  std::min(std::min(center.data().comp_id, edge.data().id), other.data().comp_id);
   if (debug)
     std::cout<<"Setting edge id to: " << edge.data().comp_id << " from: " << center.id() << " to: " << other.id() << std::endl;
   return edge.data().comp_id;
}


//find min component of two edges
void bond_percolation_combine(unsigned int& v1, const unsigned int& v2) {
    v1 = std::min(v1, v2);
    if (debug)
      std::cout<<"Comparing two edge ids: " << v1 << " : " << v2 << std::endl;
}

//the main update function

void bond_percolation_function(engine_type::context_type& context,
                  graph_type::vertex_type& vertex
#ifndef SYNC
                  , const engine_type::message_type& unused) {
#else
   ){
#endif

     uint comp_id = vertex.data().comp_id;
     vertex.data().comp_id =  context.map_reduce<unsigned int>(BOND_PERCOLATION_MAP_REDUCE, graphlab::ALL_EDGES);
     if (debug && comp_id != vertex.data().comp_id)
       std::cout<<"node: " << vertex.id() << " min edge component found: " << vertex.data().comp_id << std::endl;
#ifndef SYNC
     if (comp_id != vertex.data().comp_id)
       context.broadcast_signal(graphlab::ALL_EDGES);
#endif
}


struct model_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;

  std::string save_vertex(const vertex_type& vertex) const {
    return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return boost::lexical_cast<std::string>(edge.data().id) + " " + boost::lexical_cast<std::string>(edge.data().comp_id) + std::string("\n");
  }
};

//reduce sizes of connected components
struct label_counter {
  std::map<uint, uint> counts;

  label_counter() { }

  explicit label_counter(size_t labelId) {
    counts[labelId] = 1;
  }

  label_counter& operator+=(const label_counter& other) {
    std::map<uint, uint>::const_iterator iter = other.counts.begin();
    while(iter != other.counts.end()) {
      counts[iter->first] += iter->second;
      ++iter;
    }
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << counts;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> counts;
  }
};


label_counter get_comp_ide(const graph_type::edge_type& edge) {
  return label_counter(edge.data().comp_id);
}
label_counter get_comp_idv(const graph_type::vertex_type& vertex) {
  return label_counter(vertex.data().comp_id);
}





int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description =
      "Compute connected component - by edges";
  graphlab::command_line_options clopts(description);
  std::string input_dir;
  clopts.attach_option("graph", input_dir,
                       "The directory containing the graph file");
  clopts.add_positional("graph");
  clopts.attach_option("output_file", output_file,
                       "The prefix (folder and filename) to save output_file.");
  clopts.attach_option("max_iter", max_iter, "max number of iterations");
  clopts.attach_option("debug", debug, "debug (verbose) mode");
  clopts.attach_option("p", p, "percentage for active node (optional)");
  clopts.attach_option("n", n, "total number of nodes (optional)");

  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  //dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer;
  graph_type graph(dc, clopts);
  graph.load(input_dir, graph_loader);
  dc.cout() << "Loading graph. Finished in "
            << timer.current_time() << std::endl;

  //dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  //dc.cout() << "Finalizing graph. Finished in "
  //          << timer.current_time() << std::endl;


  /*dc.cout()
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
  */

  if (debug)
    //omp_set_num_threads(1);

  /* THE MAIN LOOP */
  dc.cout() << "Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts);
#ifndef SYNC
  engine.set_vertex_program(bond_percolation_function);
#endif
  engine.register_map_reduce(BOND_PERCOLATION_MAP_REDUCE, bond_percolation_map, bond_percolation_combine);

  engine.signal_all();
#ifndef SYNC
  engine.wait();
#endif

#ifdef SYNC
  /* FOR EACH ITERATION */
  for (int i=0; i< max_iter; i++){
     /* PERFORM UPDATE FUNCTION */
     engine.parfor_all_local_vertices(bond_percolation_function);
     /* WAIT UNTIL COMPLETION */
     engine.wait();
     /* CHECK FOR CONVERGENCE */
     size_t diff = graph.map_reduce_edges<size_t>(count_component);
     dc.cout() << "iter = " << i << " diff= " << diff << std::endl;
     if (diff == 0)
       break;
  }
#endif

  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
            << std::endl
            << "Final Runtime (seconds):   " << runtime;

  dc.cout() << " Num updates: " << engine.num_updates() << std::endl;
  // Make output_file ---------------------------------------------------------
  if(!output_file.empty()) {
    std::cout << "Saving output_file" << std::endl;
    const bool gzip_output = false;
    const bool save_vertices = false;
    const bool save_edges = true;
    const size_t threads_per_machine = 1;

    //save the output
    graph.save(output_file, model_saver(), gzip_output, save_vertices, save_edges, threads_per_machine);
  }

  //take statistics
  label_counter state = graph.map_reduce_edges<label_counter>(get_comp_ide);
  label_counter statv = graph.map_reduce_vertices<label_counter>(get_comp_idv);
  if (p > 0)
    dc.cout() << "site fraction p= " << p << std::endl;
  if (n > 0){
    dc.cout() << "n=" << n*p << std::endl;
    dc.cout() << "isolated sites: " << p*n-graph.num_vertices() << std::endl;
  }
  dc.cout() << "Number of sites: " << graph.num_vertices() << std::endl;
  dc.cout() << "Number of bonds: " << graph.num_edges() << std::endl;
  if (n){
  dc.cout() << "Percentage of sites: " << (double)graph.num_vertices() / n << std::endl;
  dc.cout() << "Percentage of bonds: " << (double)graph.num_edges() / (2.0*n) << std::endl;
  }
  dc.cout() << "SITES RESULT:\nsize\tcount\n";
  std::map<uint,uint> final_countsv;
  uint total_sites = 0;
  for (std::map<uint, uint>::const_iterator iter = statv.counts.begin();
      iter != statv.counts.end(); iter++) {
      //dc.cout() << iter->first << "\t" << iter->second << "\n";
      final_countsv[iter->second] += 1;
      total_sites += iter->second;
  }
  for (std::map<uint, uint>::const_iterator iter = final_countsv.begin();
      iter != final_countsv.end(); iter++) {
      dc.cout() << iter->first << "\t" << iter->second << "\n";
  }

  dc.cout() << "BONDS RESULT:\nsize\tcount\n";
  std::map<uint,uint> final_countse;
  uint total_bonds = 0;
  for (std::map<uint, uint>::const_iterator iter = state.counts.begin();
      iter != state.counts.end(); iter++) {
      //dc.cout() << iter->first << "\t" << iter->second << "\n";
      final_countse[iter->second] += 1;
      total_bonds += iter->second;
  }
  for (std::map<uint, uint>::const_iterator iter = final_countse.begin();
      iter != final_countse.end(); iter++) {
      dc.cout() << iter->first << "\t" << iter->second << "\n";
  }
  assert(total_sites == graph.num_vertices());
  assert(total_bonds == graph.num_edges());

  //shutdown MPI
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



