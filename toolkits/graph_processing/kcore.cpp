#include <boost/unordered_set.hpp>
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>
/**
 *
 * In this program we implement the "k-core" decomposition algorithm.
 * We use a parallel variant of
 * 
 * V. Batagelj and M. Zaversnik, An O(m) algorithm for cores
 * decomposition of networks,
 *
 *  - Essentially, recursively remove everything with degree 1
 *  - Then recursively remove everything with degree 2
 *  - etc.
 */

/*
 * Each vertex maintains a "degree" count. If this value
 * is 0, the vertex is "deleted"
 */
typedef int vertex_data_type;

/*
 * Don't need any edges
 */
typedef graphlab::empty edge_data_type;

/*
 * Define the type of the graph
 */
typedef graphlab::distributed_graph<vertex_data_type,
                                    edge_data_type> graph_type;


size_t CURRENT_K;
                                    
class k_core :
      public graphlab::ivertex_program<graph_type,
                                       graphlab::empty, // gathers are integral
                                       int>,   // messages are integral
      public graphlab::IS_POD_TYPE  {
public:
  int msg;
  bool just_deleted;

  k_core():msg(0),just_deleted(false) {
  }
  
  void recv_message(icontext_type& context,
                    const vertex_type& vertex,
                    const message_type& message) {
    msg = message;
    just_deleted = false;
  }

  // gather is never invoked
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }

  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& unused) {
    if (vertex.data() > 0) {
      vertex.data() -= msg;
      if (vertex.data() < CURRENT_K) {
        just_deleted = true;
        vertex.data() = 0;
      }
    }
  } 

  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return just_deleted ?
                 graphlab::ALL_EDGES : graphlab::NO_EDGES;
  }

  void scatter(icontext_type& context,
              const vertex_type& vertex,
              edge_type& edge) const {
    vertex_type other = edge.source().id() == vertex.id() ?
                          edge.target() : edge.source();
    if (other.data() > 0) {
      context.signal(edge.source().id() == vertex.id() ?
                      edge.target() : edge.source(), 1);
    }
  }
  
};

typedef graphlab::synchronous_engine<k_core> engine_type;

void initialize_vertex_values(graph_type::vertex_type &v) {
  v.data() = v.num_in_edges() + v.num_out_edges();
}

graphlab::empty signal_vertices_at_k(engine_type::icontext_type& ctx,
                          graph_type::vertex_type &vertex) {
  if (vertex.data() > 0 && vertex.data() < CURRENT_K) {
    ctx.signal(vertex, 0);
  }
  return graphlab::empty();
}


size_t count_active_vertices(graph_type::vertex_type &vertex) {
  return vertex.data() > 0;
}

size_t double_count_active_edges(graph_type::vertex_type &vertex) {
  return (size_t) vertex.data();
}



struct save_core_at_k{
  std::string save_vertex(graph_type::vertex_type) { return ""; }
  std::string save_edge(graph_type::edge_type e) {
    if (e.source().data() > 0 && e.target().data() > 0) {
      return graphlab::tostr(e.source().id()) + "\t" +
             graphlab::tostr(e.target().id()) + "\n";
    }
    else return "";
  }
};
    
int main(int argc, char** argv) {
  std::cout << "Computes a k-core decomposition of a graph.\n";

  graphlab::command_line_options clopts("Exact Triangle Counting");
  std::string prefix, format;
  size_t kmin = 0;
  size_t kmax = 100;
  std::string savecores;
  clopts.attach_option("graph",
                       &prefix, prefix,
                       "Graph Prefix");
  clopts.attach_option("format",
                       &format, format,
                       "The graph format");
  clopts.attach_option("kmin",
                       &kmin, kmin,
                       "Compute the k-Core for k the range [kmin,kmax]");
  clopts.attach_option("kmax",
                       &kmax, kmax,
                       "Compute the k-Core for k the range [kmin,kmax]");
  clopts.attach_option("savecores",
                       &savecores, savecores,
                       "If non-empty, will save tsv of each core with prefix [savecores].K");

  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (prefix == "") {
    std::cout << "--graph is not optional\n";
    return EXIT_FAILURE;
  }
  else if (format == "") {
    std::cout << "--format is not optional\n";
    return EXIT_FAILURE;
  }
  else if (kmax < kmin) {
    std::cout << "kmax must be at least as large as kmin\n";
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

  graphlab::synchronous_engine<k_core> engine(dc, graph, clopts);
  graph.transform_vertices(initialize_vertex_values);
  // create engine to generate neighborhood
  for (CURRENT_K = kmin; CURRENT_K <= kmax; CURRENT_K++) {
    engine.map_reduce_vertices<graphlab::empty>(signal_vertices_at_k);
    engine.start();
    size_t numv = graph.map_reduce_vertices<size_t>(count_active_vertices);
    size_t nume = graph.map_reduce_vertices<size_t>(double_count_active_edges) / 2;
    dc.cout() << "K=" << CURRENT_K << ":  #V = "
              << numv << "   #E = " << nume << std::endl;
    if (savecores != "") {
      
      graph.save(savecores + "." + graphlab::tostr(CURRENT_K),
                 save_core_at_k(),
                 false, /* no compression */ 
                 false, /* do not save vertex */
                 true, /* save edge */ 
                 1); /* one file per machine */
    }
  }
  
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main

