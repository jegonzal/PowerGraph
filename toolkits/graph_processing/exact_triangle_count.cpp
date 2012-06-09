#include <boost/unordered_set.hpp>
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>
/**
 *  
 * In this program we implement the "hash-table" version of the
 * "edge-iterator" algorithm described in
 * 
 *    T. Schank. Algorithmic Aspects of Triangle-Based Network Analysis.
 *    Phd in computer science, University Karlsruhe, 2007.
 *
 * The procedure is quite straightforward:
 *   - each vertex maintains a list of all of its neighbors in a hash table.
 *   - For each edge (u,v) in the graph, count the number of intersections
 *     of the neighbor set on u and the neighbor set on v. (of course,
 *     excluding the vertex v in the neighbor set of u, and vice versa).
 *   - We store the size of the intersection on the edge.
 * 
 * This will count every triangle exactly 3 times. Summing across all the
 * edges and dividing by 3 gives the desired result.
 *
 * The preprocessing stage take O(|E|) time, and it has been shown that this
 * algorithm takes $O(|E|^(3/2))$ time.
 *
 * We then introduce a optimization that is similar to the "forward" algorithm
 * described in thesis above. Instead of maintaining a complete list of all
 * neighbors, each vertex only maintains a list of all neighbors with
 * ID greater than itself. This implicitly generates a topological sort
 * of the graph.
 *
 * Then you can see that each triangle
 *
 * \verbatim
  
     A----->C
     |     ^
     |   /
     v /
     B
   
 * \endverbatim
 * Must be counted only once. (Only when processing edge AB, can one
 * observe that A and B have intersecting out-neighbor sets).
 *
 *
 * \note The implementation here is built to be easy to understand
 * and not necessarily optimal. In particular the unordered_set is slow
 * for small number of entries. A union of a small set which does not rely
 * on malloc, and an unordered_set is probably much more efficient.
 */

/*
 * Each vertex maintains a list of all its neighbors.
 */
typedef boost::unordered_set<graphlab::vertex_id_type> vertex_data_type;

/*
 * Each edge is simply a counter of triangles
 */
typedef size_t edge_data_type;

// To collect the set of neighbors, we need a message type which is
// basically a set of vertex IDs

struct set_union_gather {
  boost::unordered_set<graphlab::vertex_id_type> vid_set;
  
  set_union_gather& operator+=(const set_union_gather& other) {
    foreach(graphlab::vertex_id_type othervid, other.vid_set) {
      vid_set.insert(othervid);
    }
    return *this;
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << vid_set;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> vid_set;
  }
};

/*
 * Define the type of the graph
 */
typedef graphlab::distributed_graph<vertex_data_type,
                                    edge_data_type> graph_type;


class triangle_count :
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

  
  gather_type gather(icontext_type& context,
                     const vertex_type& vertex,
                     edge_type& edge) const {
    set_union_gather gather;
    // Insert the opposite end of the edge IF the opposite end has
    // ID greater than the current vertex
    vertex_id_type otherid = edge.source().id() == vertex.id() ?
                             edge.target().id() : edge.source().id();
    if (otherid > vertex.id()) gather.vid_set.insert(otherid);
    return gather;
  }

  /*
   * Simply store the gather data on the vertex 
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& neighborhood) {
    vertex.data() = neighborhood.vid_set;
  } // end of apply

  /*
   * Scatter over all edges to compute the intersection.
   * I only need to touch each edge once, so if I scatter just on the
   * out edges, that is sufficient.
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::OUT_EDGES;
  }

  static size_t count_set_intersect(
               const boost::unordered_set<vertex_id_type>& smaller_set,
               const boost::unordered_set<vertex_id_type>& larger_set) {
    size_t count = 0;
    foreach(vertex_id_type vid, smaller_set) {
      count += larger_set.count(vid);
    }
    return count;
  }
  
  void scatter(icontext_type& context,
              const vertex_type& vertex,
              edge_type& edge) const {
    const vertex_data_type& srclist = edge.source().data();
    const vertex_data_type& targetlist = edge.target().data();
    if (srclist.size() >= targetlist.size()) {
      edge.data() = count_set_intersect(targetlist, srclist);
    }
    else {
      edge.data() = count_set_intersect(srclist, targetlist);
    }
  }
};


size_t get_edge_data(const graph_type::edge_type& e) {
  return e.data();
}


int main(int argc, char** argv) {
  std::cout << "This program counts the exact number of triangles in the "
            "provided graph.\n";

  graphlab::command_line_options clopts("Exact Triangle Counting");
  std::string prefix, format;
  clopts.attach_option("graph",
                       &prefix, prefix,
                       "Graph Prefix");
  clopts.attach_option("format",
                       &format, format,
                       "The graph format");
  
  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (prefix == "") {
    std::cout << "--graph is not optional\n";
    return EXIT_FAILURE;
  }
  else if (format == "") {
    std::cout << "--format is not optional\n";
    return EXIT_FAILURE;
  }
  
  // Initialize control plane using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  // load graph
  graph_type graph(dc, clopts);
  graph.load_format(prefix, format);
  graph.finalize();
  dc.cout() << "Number of vertices: " << graph.num_vertices() << std::endl
            << "Number of edges:    " << graph.num_edges() << std::endl;

  graphlab::timer ti;
  
  // create engine to generate neighborhood
  dc.cout() << "Counting Triangles..." << std::endl;
  graphlab::synchronous_engine<triangle_count> engine(dc, graph, clopts);
  engine.signal_all();
  engine.start();

  size_t count =
    graph.map_reduce_edges<size_t>(get_edge_data);

  dc.cout() << count << " Triangles"  << std::endl;
  dc.cout() << "Counted in " << ti.current_time() << " seconds" << std::endl;
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main

