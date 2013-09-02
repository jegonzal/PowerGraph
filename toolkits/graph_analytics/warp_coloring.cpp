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

#include <vector>
#include <string>
#include <fstream>

#include <graphlab.hpp>
#include <graphlab/engine/gl3engine.hpp>
#include <graphlab/macros_def.hpp>
using namespace graphlab;

#define UNIQUE_COLOR_MAP_REDUCE 0
#define SIGNAL_IF_CHANGE 1

// The vertex data is the color of the vertex
typedef graphlab::vertex_id_type color_type;



struct edge_data_type : public graphlab::IS_POD_TYPE {
  bool dirty;
  bool locked;
  char owned_by; // 0 is source, 1 is other
  bool requested; // 1 is requested by other party;
};

std::ostream& operator<<(std::ostream& o, const edge_data_type& e) {
  o << e.dirty << " " << e.locked << " " << (int)e.owned_by << " " << e.requested;
  return o;
}

// The graph type is determined by the vertex and edge data types
typedef distributed_graph<color_type, edge_data_type> graph_type;

typedef gl3engine<graph_type> engine_type;

engine_type* eng;

bool EDGE_CONSISTENT = false;



/*
 * This is the gathering type which accumulates an (unordered) set of
 * all neighboring colors
 * It is a simple wrapper around a boost::unordered_set with
 * an operator+= which simply performs a set union.
 *
 * This struct can be significantly accelerated for small sets.
 * Small collections of vertex IDs should not require the overhead
 * of the unordered_set.
 */
struct set_union_gather {
  boost::unordered_set<color_type> colors;

  /*
   * Combining with another collection of vertices.
   * Union it into the current set.
   */
  set_union_gather& operator+=(const set_union_gather& other) {
    foreach(graphlab::vertex_id_type othervid, other.colors) {
      colors.insert(othervid);
    }
    return *this;
  }

  // serialize
  void save(graphlab::oarchive& oarc) const {
    oarc << colors;
  }

  // deserialize
  void load(graphlab::iarchive& iarc) {
    iarc >> colors;
  }
};

/**************************************************************************/
/*                                                                        */
/*                           Coloring Functions                           */
/*                                                                        */
/**************************************************************************/

set_union_gather unique_color_map(const graph_type::vertex_type& center,
                                  graph_type::edge_type& e,
                                  const graph_type::vertex_type& other) {
  set_union_gather gather;
  color_type other_color = other.data();
  gather.colors.insert(other_color);
  return gather;
}

void unique_color_combine(set_union_gather& v1, const set_union_gather& v2) {
  v1 += v2;
}


void schedule_neighbors_if_change(const graph_type::vertex_type& center,
                                  graph_type::edge_type& e,
                                  const graph_type::vertex_type& other) {
  if (center.data() == other.data()) {
    eng->get_context().signal(other);
  }
}

/**************************************************************************/
/*                                                                        */
/*                         Chandy Misra Functions                         */
/*                                                                        */
/**************************************************************************/
void initialize_chandy_misra(graph_type::edge_type& edge) {
  edge.data().dirty = true;
  edge.data().locked = false;
  edge.data().owned_by = edge.source().id() < edge.target().id() ? 0 : 1;
  edge.data().requested = false;
//   std::cout << edge.source().id() << "->" << edge.target().id() << " " << edge.data() << "\n";
}


#define LOCK_IF_OWNED 3
size_t lock_if_owned_map(const graph_type::vertex_type& center,
                         graph_type::edge_type& e,
                         const graph_type::vertex_type& other) {
  char m =  (e.source().id() == center.id()) ? 0 : 1;
  if (e.data().owned_by == m) e.data().locked = true;
  else if (e.data().dirty && e.data().locked == false) {
    e.data().owned_by = m;
    e.data().dirty = false;
    e.data().requested = false;
    e.data().locked = true;
  } else {
    e.data().requested = true;
  }
//   std::cout << "Lock If Owned: " << center.id() << ": Fork = "<< e.source().id() << "->" << e.target().id() << " "  << e.data() << "\n";
  return (e.data().owned_by == m && e.data().locked);
}

void lock_if_owned_combine(size_t& v1, const size_t& v2) {
  v1 += v2;
}


#define STOP_EATING 4
void stop_eating(const graph_type::vertex_type& center,
                  graph_type::edge_type& e,
                  const graph_type::vertex_type& other) {
  e.data().locked = false;
  e.data().dirty = true;
  if (e.data().requested) {
    // switch owner
    e.data().owned_by = !e.data().owned_by;
    e.data().dirty = false;
    e.data().requested = false;
  }

//   std::cout << "Stop Eating: " << center.id() << ": Fork = "<< e.source().id() << "->" << e.target().id() << " "  << e.data() << "\n";
}


#define UNLOCK_FORKS_MAINTAIN_REQUEST 5
void unlock_forks_maintain_request(const graph_type::vertex_type& center,
                                   graph_type::edge_type& e,
                                   const graph_type::vertex_type& other) {
  char m =  (e.source().id() == center.id()) ? 0 : 1;
  if (e.data().owned_by == m && e.data().locked) {
    e.data().locked = false;
  }

  if (e.data().owned_by == m && e.data().dirty && e.data().requested) {
    e.data().owned_by = !e.data().owned_by;
    e.data().dirty = false;
    e.data().requested = true;
  }

//   std::cout << "Release: " << center.id() << ": Fork = " << e.source().id() << "->" << e.target().id() << " " << e.data() << "\n";
}





void update_function(engine_type::context_type& context,
                     graph_type::vertex_type& vertex,
                     const engine_type::message_type& unused) {

  // acquire locks
  if (EDGE_CONSISTENT) {
    size_t expected_num_locks = vertex.num_in_edges() + vertex.num_out_edges();
    while(1) {
      size_t numnbr = context.map_reduce<size_t>(LOCK_IF_OWNED, ALL_EDGES);
      if (numnbr == expected_num_locks) {
        break;
      }
      else {
        context.edge_transform(UNLOCK_FORKS_MAINTAIN_REQUEST, ALL_EDGES);
      }
      graphlab::fiber_group::yield();
    }
  }
  set_union_gather neighborhood =
      context.map_reduce<set_union_gather>(UNIQUE_COLOR_MAP_REDUCE,
                                           ALL_EDGES);

  bool color_changed = false;

  size_t neighborhoodsize = neighborhood.colors.size();
  for (color_type curcolor = 0; curcolor < neighborhoodsize + 1; ++curcolor) {
    if (neighborhood.colors.count(curcolor) == 0) {
      vertex.data() = curcolor;
      break;
    }
  }
  if (EDGE_CONSISTENT) context.edge_transform(STOP_EATING, ALL_EDGES);

  context.edge_transform(SIGNAL_IF_CHANGE, ALL_EDGES, false);
}


/**************************************************************************/
/*                                                                        */
/*                         Validation   Functions                         */
/*                                                                        */
/**************************************************************************/
size_t validate_conflict(graph_type::edge_type& edge) {
  return edge.source().data() == edge.target().data();
}



int main(int argc, char** argv) {
  // Initialize control plain using mpi
  mpi_tools::init(argc, argv);
  distributed_control dc;

  // Parse command line options -----------------------------------------------
  command_line_options clopts("PageRank algorithm.");
  std::string graph_dir;
  std::string format = "adj";
  clopts.set_scheduler_type("fifo");
  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "The graph file format");
  clopts.attach_option("edgescope", EDGE_CONSISTENT,
                       "Set to 1 if edge consistency is to be used");
  size_t powerlaw = 0;
  clopts.attach_option("powerlaw", powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  graphlab::launch_metric_server();
  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  if(powerlaw > 0) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, 2.1, 100000000);
  }
  else if (graph_dir.length() > 0) { // Load the graph from a file
    dc.cout() << "Loading graph in format: "<< format << std::endl;
    graph.load_format(graph_dir, format);
  }
  else {
    dc.cout() << "graph or powerlaw option must be specified" << std::endl;
    clopts.print_description();
    return 0;
  }
  // must call finalize before querying the graph
  graph.finalize();

  graph.transform_edges(initialize_chandy_misra);

  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;

  // Running The Engine -------------------------------------------------------
  engine_type engine(dc, graph, clopts);
  eng = &engine;

  engine.register_map_reduce(UNIQUE_COLOR_MAP_REDUCE,
                             unique_color_map,
                             unique_color_combine);

  engine.register_edge_transform(SIGNAL_IF_CHANGE,
                                 schedule_neighbors_if_change);


  engine.register_map_reduce(LOCK_IF_OWNED,
                             lock_if_owned_map,
                             lock_if_owned_combine);

  engine.register_edge_transform(STOP_EATING,
                                 stop_eating);

  engine.register_edge_transform(UNLOCK_FORKS_MAINTAIN_REQUEST,
                                 unlock_forks_maintain_request);

  engine.set_vertex_program(update_function);
  timer ti; ti.start();
  engine.signal_all();

  engine.wait();
  dc.cout() << "Finished Running engine in " << ti.current_time()
            << " seconds." << std::endl;
  dc.cout() << engine.num_updates()
            << " updates." << std::endl;

  size_t conflict_count = graph.map_reduce_edges<size_t>(validate_conflict);
  dc.cout() << "Num conflicts = " << conflict_count << "\n";


  graphlab::stop_metric_server();

  mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


// We render this entire program in the documentation


