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


/*
 * Graph coloring algorithm, such that vertex programs are scheduled in 
 * order of vertex degree. Includes trade-off featre for determining the
 * fraction of the graph to conduct ordered execution over random execution,
 * allowing user to specify run time - colour quality trade-off
 */

#include <boost/unordered_set.hpp>
#include <graphlab.hpp>
#include <graphlab/ui/metrics_server.hpp>
#include <graphlab/macros_def.hpp>
#include <cmath>  /* for std::abs(double) */

typedef graphlab::vertex_id_type color_type;

/*
 * Vertex data: color and degree of node
 */
typedef struct {
  int color;
  int degree;

   // serialize
  void save(graphlab::oarchive& oarc) const {
    oarc << color << degree;
  }

  // deserialize
  void load(graphlab::iarchive& iarc) {
    iarc >> color >> degree;
  }

} vertex_data_type;

#define UNCOLORED -1
/*
 * no edge data
 */
typedef graphlab::empty edge_data_type;
bool EDGE_CONSISTENT = false;
bool TRADE = false;

size_t graph_size = 0;
size_t fraction = 0;
int max_degree = 0;
int low_degree = INT_MAX;
signed int current_degree;

size_t already_signalled = 0;
std::set<int> used_colors;
std::set<int> degrees;
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
/*
 * Define the type of the graph
 */
typedef graphlab::distributed_graph<vertex_data_type,
                                    edge_data_type> graph_type;


/*
 * On gather, we accumulate a set of all adjacent colors.
 */
class graph_coloring:
      public graphlab::ivertex_program<graph_type,
                                      set_union_gather>,
      /* I have no data. Just force it to POD */
      public graphlab::IS_POD_TYPE  {
public:
  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  } 

  /*
   * For each edge, figure out the ID of the "other" vertex
   * and accumulate a set of the neighborhood vertex IDs.
   */
  gather_type gather(icontext_type& context,
                     const vertex_type& vertex,
                     edge_type& edge) const {
    set_union_gather gather;
    color_type other_color = edge.source().id() == vertex.id() ?
                                 edge.target().data().color: edge.source().data().color;

    gather.colors.insert(other_color);
    return gather;
  }

  /*
   * the gather result now contains the colors in the neighborhood.
   * pick a different color and store it 
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& neighborhood) {
    // find the smallest color not described in the neighborhood
    size_t neighborhoodsize = neighborhood.colors.size();
    for (color_type curcolor = 0; curcolor < neighborhoodsize + 1; ++curcolor) {
      if (neighborhood.colors.count(curcolor) == 0) {
        used_colors.insert(curcolor);
        vertex.data().color = curcolor;
        break;
      }
    }
  }

  edge_dir_type scatter_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    if (EDGE_CONSISTENT) return graphlab::NO_EDGES;
    else return graphlab::ALL_EDGES;
  } 


  /*
   * For each edge, count the intersection of the neighborhood of the
   * adjacent vertices. This is the number of triangles this edge is involved
   * in.
   */
  void scatter(icontext_type& context,
              const vertex_type& vertex,
              edge_type& edge) const {
    // both points have different colors!
    if (edge.source().data().color == edge.target().data().color) {
      context.signal(edge.source().id() == vertex.id() ? 
                      edge.target() : edge.source());
    }
  }
};

void initialize_vertex_values(graph_type::vertex_type& v) {
  v.data().degree = v.num_out_edges();
  degrees.insert(v.data().degree);
  v.data().color = UNCOLORED;
  if (v.data().degree > max_degree)
    max_degree = v.data().degree;
  if (v.data().degree < low_degree)
    low_degree = v.data().degree;
}


/*
 * A saver which saves a file where each line is a vid / color pair
 */
struct save_colors{
  std::string save_vertex(graph_type::vertex_type v) { 
    return graphlab::tostr(v.id()) + "\t" +
           graphlab::tostr(v.data().color) + "\n";
  }
  std::string save_edge(graph_type::edge_type e) {
    return "";
  }
};

typedef graphlab::async_consistent_engine<graph_coloring> engine_type;

graphlab::empty signal_vertices_at_degree(engine_type::icontext_type& ctx,
                                     const graph_type::vertex_type& vertex) {
  if (vertex.data().degree == current_degree) {
    already_signalled++;
    ctx.signal(vertex);
  }
  return graphlab::empty();
}

graphlab::empty signal_uncolored(engine_type::icontext_type& ctx,
                                     const graph_type::vertex_type& vertex) {
  if (vertex.data().color == UNCOLORED) {
    ctx.signal(vertex);
  }
  return graphlab::empty();
}

struct max_deg_vertex_reducer: public graphlab::IS_POD_TYPE {
  size_t degree;
  graphlab::vertex_id_type vid;
  max_deg_vertex_reducer& operator+=(const max_deg_vertex_reducer& other) {
    if (degree < other.degree) {
      (*this) = other;
    }
    return (*this);
  }
};

max_deg_vertex_reducer find_max_deg_vertex(const graph_type::vertex_type vtx) {
  max_deg_vertex_reducer red;
  red.degree = vtx.num_out_edges();
  red.vid = vtx.id();
  return red;
}

/**************************************************************************/
/*                                                                        */
/*                         Validation   Functions                         */
/*                                                                        */
/**************************************************************************/
size_t validate_conflict(graph_type::edge_type& edge) {
  return edge.source().data().color == edge.target().data().color;
}

inline bool isEqual(double x, double y)
{
  const double epsilon = 1e-5;
  return std::abs(x - y) <= epsilon * std::abs(x);
}

int main(int argc, char** argv) {

  //global_logger().set_log_level(LOG_INFO);

  // Initialize control plane using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;


  dc.cout() << "This program computes a simple graph coloring of a"
            "provided graph.\n\n";

  graphlab::command_line_options clopts("Graph coloring. "
    "Given a graph, this program computes a graph coloring of the graph."
    "The Asynchronous engine is used.");
  std::string prefix, format;
  std::string output;
  float alpha = 2.1;
  size_t powerlaw = 0;
  double trade = 0;
  clopts.attach_option("graph", prefix,
                       "Graph input. reads all graphs matching prefix*");
  clopts.attach_option("format", format,
                       "The graph format");
   clopts.attach_option("output", output,
                       "A prefix to save the output.");
   clopts.attach_option("powerlaw", powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
   clopts.attach_option("alpha", alpha,
                       "Alpha in powerlaw distrubution");
  clopts.attach_option("trade", trade,
                       "Execute tradeoff version. Probability of degree execution for node (0.0 to 1.0)");
  clopts.attach_option("edgescope", EDGE_CONSISTENT,
                       "Use Locking. ");

  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (prefix.length() == 0 && powerlaw == 0) {
    clopts.print_description();
    return EXIT_FAILURE;
  }
  if (output == "") {
    dc.cout() << "Warning! Output will not be saved\n";
  }

  graphlab::launch_metric_server();
  // load graph
  graph_type graph(dc, clopts);

  if (!isEqual(0.0, trade)) {
    TRADE = true;
  }
  
  if(powerlaw > 0) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, alpha, 100000000);
  } else { // Load the graph from a file
    if (prefix == "") {
      dc.cout() << "--graph is not optional\n";
      return EXIT_FAILURE;
    }
    else if (format == "") {
      dc.cout() << "--format is not optional\n";
      return EXIT_FAILURE;
    }
    graph.load_format(prefix, format);
  }
  graph.finalize();

  dc.cout() << "Number of vertices: " << graph.num_vertices() << std::endl
    << "Number of edges:    " << graph.num_edges() << std::endl;

  graphlab::timer ti;

  
  dc.cout() << "Initialising vertex data..." <<std::endl;
  graph.transform_vertices(initialize_vertex_values);
 
  dc.cout() << "Degrees range from "<< low_degree << " to " << max_degree << std::endl;

  // create engine to count the number of triangles
  dc.cout() << "Coloring..." << std::endl;
  if (EDGE_CONSISTENT) {
    clopts.get_engine_args().set_option("factorized", false);
  } else {
    clopts.get_engine_args().set_option("factorized", true);
  } 
  graphlab::async_consistent_engine<graph_coloring> engine(dc, graph, clopts);

  //Tradeoff between ordered and random vertex execution
  if (TRADE) {
    graph_size = graph.num_vertices();
    fraction = (int) graph_size * trade;
    dc.cout() << "Degree ordered coloring for " << fraction << " in " << graph_size << " vertices." << std::endl;
  }
  for (int x = max_degree; x >= low_degree; x--){
    if (degrees.find(x) != degrees.end()) {
      current_degree = x;
      engine.map_reduce_vertices<graphlab::empty>(signal_vertices_at_degree);  
      if (TRADE) {
        //Already signalled vertices for degree ordered execution
        if(already_signalled >= fraction) {
          engine.start();
          //Signal remaining vertices randomly
          engine.map_reduce_vertices<graphlab::empty>(signal_uncolored);  
          engine.start();
          break;
        }
      }
    }
  }

  if (!TRADE) {
    engine.start();
  }

  size_t conflict_count = graph.map_reduce_edges<size_t>(validate_conflict);
  if (conflict_count > 0) {
    dc.cout() << "Still uncolored, finalising..." << std::endl;
    engine.map_reduce_vertices<graphlab::empty>(signal_uncolored);
    engine.start();
    conflict_count = graph.map_reduce_edges<size_t>(validate_conflict);
  }

  dc.cout() << "Colored in " << ti.current_time() << " seconds" << std::endl;
  dc.cout() << "Colored using " << used_colors.size() << " colors" << std::endl;

  dc.cout() << "Num conflicts = " << conflict_count << "\n";
  if (output != "") {
    graph.save(output,
              save_colors(),
              false, /* no compression */
              true, /* save vertex */
              false, /* do not save edge */
              1); /* one file per machine */
  }
  
  graphlab::stop_metric_server();

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main