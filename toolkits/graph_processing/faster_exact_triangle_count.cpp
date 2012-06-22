#include <boost/unordered_set.hpp>
#include <graphlab/util/cuckoo_set_pow2.hpp>
#include <graphlab.hpp>
#include <graphlab/ui/metrics_server.hpp>
#include <graphlab/util/cuckoo_set_pow2.hpp>
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
 *     of the neighbor set on u and the neighbor set on v.
 *   - We store the size of the intersection on the edge.
 * 
 * This will count every triangle exactly 3 times. Summing across all the
 * edges and dividing by 3 gives the desired result.
 *
 * The preprocessing stage take O(|E|) time, and it has been shown that this
 * algorithm takes $O(|E|^(3/2))$ time.
 *
 * If we only require total counts, we can introduce a optimization that is
 * similar to the "forward" algorithm
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
typedef graphlab::cuckoo_set_pow2<graphlab::vertex_id_type, 3> hash_set;
 
/*
 * Each vertex maintains a list of all its neighbors.
 * and a final count for the number of triangles it is involved in
 */
struct vertex_data_type {
  vertex_data_type():vid_set(-1,1,1),
                     num_triangles(0),has_large_neighbors(true){ }
  // A list of all its neighbors
  hash_set vid_set;
  // The number of triangles this vertex is involved it.
  // only used if "per vertex counting" is used
  uint32_t num_triangles;
  bool has_large_neighbors;
  void save(graphlab::oarchive &oarc) const {
    oarc << vid_set << num_triangles << has_large_neighbors;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> vid_set >> num_triangles >> has_large_neighbors;
  }
};


/*
 * Each edge is simply a counter of triangles
 */
typedef uint32_t edge_data_type;

// To collect the set of neighbors, we need a message type which is
// basically a set of vertex IDs

bool PER_VERTEX_COUNT = false;

// if phases are used, only vertices with id % CUR_PHASE will be considered
// in the round
unsigned short NUM_PHASES = 1;
unsigned short CUR_PHASE = 0;

// if the neighborhood size is <= this number, it will ignore the phase number
size_t MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION = 32;


/*
 * This is the gathering type which accumulates an array of
 * all neighboring vertices.
 * It is a simple wrapper around a vector with
 * an operator+= which simply performs a  +=
 */
struct set_union_gather {
  bool large_neighbors;
  std::vector<graphlab::vertex_id_type> vid_vec;

  set_union_gather():large_neighbors(false) {
  }
  /*
   * Combining with another collection of vertices.
   * Union it into the current set.
   */
  set_union_gather& operator+=(const set_union_gather& other) {
    if (other.vid_vec.size() > 0) {
      vid_vec.reserve(vid_vec.size() + other.vid_vec.size());
      foreach(graphlab::vertex_id_type othervid, other.vid_vec) {
        vid_vec.push_back(othervid);
      }
    }
    large_neighbors |= other.large_neighbors;
    return *this;
  }
  
  // serialize
  void save(graphlab::oarchive& oarc) const {
    oarc << large_neighbors << vid_vec;
  }

  // deserialize
  void load(graphlab::iarchive& iarc) {
    iarc >> large_neighbors >> vid_vec;
  }
};

/*
 * Define the type of the graph
 */
typedef graphlab::distributed_graph<vertex_data_type,
                                    edge_data_type> graph_type;


/*
 * This class implements the triangle counting algorithm as described in
 * the header. On gather, we accumulate a set of all adjacent vertices.
 * If per_vertex output is not necessary, we can use the optimization
 * where each vertex only accumulates neighbors with greater vertex IDs.
 */
class triangle_count :
      public graphlab::ivertex_program<graph_type,
                                      set_union_gather>,
      /* I have no data. Just force it to POD */
      public graphlab::IS_POD_TYPE  {
public:
  bool do_not_scatter;

  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    if (CUR_PHASE == 0  ||
        vertex.num_in_edges() + vertex.num_out_edges() 
                            >  MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION) {
      return graphlab::ALL_EDGES;
    }
    else {
      return graphlab::NO_EDGES;
    }
  } 

  /*
   * For each edge, figure out the ID of the "other" vertex
   * and accumulate a set of the neighborhood vertex IDs.
   */
  gather_type gather(icontext_type& context,
                     const vertex_type& vertex,
                     edge_type& edge) const {
    set_union_gather gather;
    graphlab::vertex_id_type otherid = edge.target().id() == vertex.id() ?
                                       edge.source().id() : edge.target().id();
    if (CUR_PHASE == 0) {
      // collect everything which matches the phase number.
      // An optimization if current has few neighbors, we collect it all
      bool cur_is_below_count =  vertex.num_in_edges() + vertex.num_out_edges() 
                                 <=  MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION;
      if (NUM_PHASES == 1 || 
          cur_is_below_count || otherid % NUM_PHASES == CUR_PHASE) {
        if (PER_VERTEX_COUNT || otherid > vertex.id()) {
          gather.vid_vec.push_back(otherid);
        } 
      }
    }
    else {
      if (NUM_PHASES == 1 || otherid % NUM_PHASES == CUR_PHASE) {
        if (PER_VERTEX_COUNT || otherid > vertex.id()) {
          gather.vid_vec.push_back(otherid);
        }
      }

    }
                         
    bool nbr_is_above_count = false;
    if (edge.target().id() != vertex.id()) {
      const vertex_type& othervtx = edge.target();
      nbr_is_above_count = othervtx.num_in_edges() + othervtx.num_out_edges() 
                              >  MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION;
    }
    // large neighbors optimization
    // we do not even issue the scatter if all out neighbors are small
    gather.large_neighbors = nbr_is_above_count;
    return gather;
  }

  /*
   * the gather result now contains the vertex IDs in the neighborhood.
   * store it on the vertex. 
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& neighborhood) {
    // we overwrite the local map only if a gather is performed
    // a gather is performed if this is phase 0 or this vertex
    // has alot of neighbors.
    bool gather_performed = CUR_PHASE == 0  ||
                  vertex.num_in_edges() + vertex.num_out_edges() 
                                    >  MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION;
    if (CUR_PHASE == 0) {
      vertex.data().has_large_neighbors = neighborhood.large_neighbors;
    }

    do_not_scatter = false;
    if (gather_performed) {
      vertex.data().vid_set.clear();
      vertex.data().vid_set.reserve(neighborhood.vid_vec.size() * 2);
      // we performed a gather
      foreach(graphlab::vertex_id_type vid, neighborhood.vid_vec) {
        vertex.data().vid_set.insert(vid);
      }
      do_not_scatter = vertex.data().vid_set.size() == 0;
    }
    else {
      // There is no need to scatter if I have no large neighbors
      // and I did not perform any gathers
      do_not_scatter = !vertex.data().has_large_neighbors;
    }
  } // end of apply

  /*
   * Scatter over all edges to compute the intersection.
   * I only need to touch each edge once, so if I scatter just on the
   * out edges, that is sufficient.
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    if (do_not_scatter) return graphlab::NO_EDGES;
    else return graphlab::OUT_EDGES;
  }


  /*
   * Computes the size of the intersection of two unordered sets
   */
  static uint32_t count_set_intersect(
               const hash_set& smaller_set,
               const hash_set& larger_set) {
    if (smaller_set.size() == 0) return 0;
    uint32_t count = 0;
    foreach(vertex_id_type vid, smaller_set) {
      count += larger_set.count(vid);
    }
    return count;
  }

  /*
   * For each edge, count the intersection of the neighborhood of the
   * adjacent vertices. This is the number of triangles this edge is involved
   * in.
   */
  void scatter(icontext_type& context,
              const vertex_type& vertex,
              edge_type& edge) const {

    vertex_type othervtx = edge.target();
    // 
    bool cur_is_above_count = vertex.num_in_edges() + vertex.num_out_edges() 
                              >  MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION;
    bool nbr_is_above_count = othervtx.num_in_edges() + othervtx.num_out_edges() 
                              >  MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION;

    if (CUR_PHASE == 0  || cur_is_above_count || nbr_is_above_count) {
      const vertex_data_type& srclist = edge.source().data();
      const vertex_data_type& targetlist = edge.target().data();
      if (srclist.vid_set.size() >= targetlist.vid_set.size()) {
        edge.data() += count_set_intersect(targetlist.vid_set, srclist.vid_set);
      }
      else {
        edge.data() += count_set_intersect(srclist.vid_set, targetlist.vid_set);
      }
    }
  }
};



/*
 * This class is used in a second engine call if per vertex counts are needed.
 * The number of triangles a vertex is involved in can be computed easily
 * by summing over the number of triangles each adjacent edge is involved in
 * and dividing by 2. 
 */
class get_per_vertex_count :
      public graphlab::ivertex_program<graph_type, size_t>,
      /* I have no data. Just force it to POD */
      public graphlab::IS_POD_TYPE  {
public:
  // Gather on all edges
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }
  // We gather the number of triangles each edge is involved in
  size_t gather(icontext_type& context,
                     const vertex_type& vertex,
                     edge_type& edge) const {
    return edge.data();
  }

  /* the gather result is the total sum of the number of triangles
   * each adjacent edge is involved in . Dividing by 2 gives the
   * desired result.
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& num_triangles) {
    vertex.data().num_triangles = num_triangles / 2;
  }

  // No scatter
  edge_dir_type scatter_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }


};

typedef graphlab::synchronous_engine<triangle_count> engine_type;

/* Used to sum over all the edges in the graph in a
 * map_reduce_edges call
 * to get the total number of triangles
 */
size_t get_edge_data(const graph_type::edge_type& e) {
  return e.data();
}

graphlab::empty signal_large_vertices(engine_type::icontext_type& context,
                                      const graph_type::vertex_type& vertex) {
  if (vertex.num_in_edges() + vertex.num_out_edges() >= 
                MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION || 
      vertex.data().has_large_neighbors ) {
    context.signal(vertex);
  }
  return graphlab::empty();
}


/*
 * A saver which saves a file where each line is a vid / # triangles pair
 */
struct save_triangle_count{
  std::string save_vertex(graph_type::vertex_type v) { 
    double nt = v.data().num_triangles();
    double n_followed = v.data().num_out_edges();
    double n_following = v.data().num_in_edges();

    return graphlab::tostr(v.id()) + "\t" +
           graphlab::tostr(v.data().num_triangles) + "\t" +
           graphlab::tostr(v.data().n_followed) + "\t" + 
           graphlab::tostr(v.data().n_following) + "\n";
  }
  std::string save_edge(graph_type::edge_type e) {
    return "";
  }
};


int main(int argc, char** argv) {
  std::cout << "This program counts the exact number of triangles in the "
            "provided graph.\n\n";

  graphlab::command_line_options clopts("Exact Triangle Counting. "
    "Given a graph, this program computes the total number of triangles "
    "in the graph. An option (per_vertex) is also provided which "
    "computes for each vertex, the number of triangles it is involved in."
    "The algorithm assumes that each undirected edge appears exactly once "
    "in the graph input. If edges may appear more than once, this procedure "
    "will over count.");
  std::string prefix, format;
  std::string per_vertex;
  clopts.attach_option("graph",
                       &prefix, prefix,
                       "Graph input. reads all graphs matching prefix*");
  clopts.attach_option("format",
                       &format, format,
                       "The graph format");
  clopts.attach_option("phases",
                       &NUM_PHASES, NUM_PHASES,
                       "cuts up the execution to run over multiple."
                       "phases. Useful if memory requirements are high");
  clopts.attach_option("phase_0_nbr",
                       &MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION , 
                       MINIMUM_NBR_SIZE_FOR_PHASE_COLLECTION,
                       "First phase will collect all vertices with at most"
                       "this number of neighbors. Only meaningful if"
                       "phases > 1");
  clopts.attach_option("per_vertex",
                       &per_vertex, per_vertex,
                       "If not empty, will count the number of "
                       "triangles each vertex belongs to and "
                       "save to file with prefix \"[per_vertex]\". "
                       "The algorithm used is slightly different "
                       "and thus will be a little slower");
  
  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (prefix == "") {
    std::cout << "--graph is not optional\n";
    return EXIT_FAILURE;
  }
  else if (format == "") {
    std::cout << "--format is not optional\n";
    return EXIT_FAILURE;
  }

  if (NUM_PHASES == 0) {
    std::cout << "phases cannot be 0\n";
    return EXIT_FAILURE;
  }

  if (per_vertex != "") PER_VERTEX_COUNT = true;
  // Initialize control plane using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  graphlab::launch_metric_server();
  // load graph
  graph_type graph(dc, clopts);
  graph.load_format(prefix, format);
  graph.finalize();
  dc.cout() << "Number of vertices: " << graph.num_vertices() << std::endl
            << "Number of edges:    " << graph.num_edges() << std::endl;

  graphlab::timer ti;
  
  // create engine to count the number of triangles
  dc.cout() << "Counting Triangles..." << std::endl;
  engine_type engine(dc, graph, clopts);
  for (CUR_PHASE = 0; CUR_PHASE < NUM_PHASES; ++CUR_PHASE) {
    if (CUR_PHASE == 0) engine.signal_all();
    else engine.map_reduce_vertices<graphlab::empty>(signal_large_vertices);
    engine.start();
  }

  dc.cout() << "Counted in " << ti.current_time() << " seconds" << std::endl;

  if (PER_VERTEX_COUNT == false) {
    size_t count = graph.map_reduce_edges<size_t>(get_edge_data);
    dc.cout() << count << " Triangles"  << std::endl;
  }
  else {
    graphlab::synchronous_engine<get_per_vertex_count> engine(dc, graph, clopts);
    engine.signal_all();
    engine.start();
    graph.save(per_vertex,
            save_triangle_count(),
            false, /* no compression */
            true, /* save vertex */
            false, /* do not save edge */
            1); /* one file per machine */

  }
  
  graphlab::stop_metric_server_on_eof();

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main

