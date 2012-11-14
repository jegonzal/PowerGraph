#include <graphlab.hpp>

namespace pagerank {
  typedef std::vector<int> gather_type;
}
namespace std {
// The gather type is vector of float = numdocs * <topic samples>
inline pagerank::gather_type& 
  operator+=(pagerank::gather_type& lvalue, const pagerank::gather_type& rvalue) {
  for (size_t i = 0; i < lvalue.size(); ++i)
    lvalue[i] += rvalue[i];
      return lvalue;
}
}

namespace pagerank {
// Global random reset probability
float RESET_PROB = 0.15;
float TOLERANCE = 1.0E-2;
int DEFAULT_NDOCS = 20;
int DEFAULT_TOPICVAL = 10;
int NPAGES = 0;
int NTOPICS = 50;

/**
 * \brief Create a pagerank changes event tracker which is reported in
 * the GraphLab metrics dashboard.
 */
DECLARE_EVENT(PAGERANK_CHANGES);

struct top_pages_type {
  graphlab::mutex lock;
  std::string json_string;
  top_pages_type() : 
    json_string("{\n" + json_header_string() + "\tvalues: [] \n }") { }
  inline std::string json_header_string() const {
    return
      "\t\"npages\": " + graphlab::tostr(NPAGES) + ",\n" + 
      "\t\"ntopics\": " + graphlab::tostr(NPAGES) + ",\n" + 
      "\t\"default_ndocs\": " + graphlab::tostr(DEFAULT_NDOCS) + ",\n" +
      "\t\"default_topicval\": " + graphlab::tostr(DEFAULT_TOPICVAL) + ",\n"; 
  } // end of json header string
} TOP_PAGES;

std::pair<std::string, std::string>
pagerank_callback(std::map<std::string, std::string>& varmap) {
  TOP_PAGES.lock.lock();
  const std::pair<std::string, std::string>
    pair("text/html",TOP_PAGES.json_string);
  TOP_PAGES.lock.unlock();
  return pair;
}


// Graph Types
// ============================================================================
struct user_feature {
   size_t join_key; // hash value of the user name, used as vertex join key
   float rank;
   int numdocs;
   std::vector<int> topics;
   user_feature() : join_key(0), rank(1.0),
                              numdocs(DEFAULT_NDOCS),
                              topics(std::vector<int>(NTOPICS, DEFAULT_TOPICVAL)) { }
   user_feature(size_t join_key) : join_key(join_key), rank(1.0),
                              numdocs(DEFAULT_NDOCS),
                              topics(std::vector<int>(NTOPICS, DEFAULT_TOPICVAL)) { }
  void save(graphlab::oarchive& oarc) const {
    oarc << join_key << rank << numdocs << topics;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> join_key >> rank >> numdocs >> topics;
  }
};

// The vertex data is its pagerank 
typedef user_feature vertex_data_type;
// The edge data is its the marginal transition probability 
typedef float edge_data_type;
// The gather type (for compute_transit vertex program) is the user topic vector. 
typedef std::vector<int> gather_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

void init_vertex(graph_type::vertex_type& vertex) { 
  vertex.data().numdocs = DEFAULT_NDOCS; 
  std::fill(vertex.data().topics.begin(), vertex.data().topics.end(), 1);
  vertex.data().rank = 1.0;
}

/*
 * Load a table of userid -> username, used for create vertex_join key
 */
bool vertex_line_parser(graph_type& graph,
                        const std::string& filename,
                        const std::string& textline) {
  ASSERT_FALSE(textline.empty());
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  graphlab::vertex_id_type vid(-1);
  std::string name;

  const bool success = qi::phrase_parse
    (textline.begin(), textline.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(vid) = qi::_1] >> -qi::char_(',') >>  name 
      )
     ,
     //  End grammar
     ascii::space); 
  if(!success) return false;  

  size_t join_key = boost::hash<std::string>()(name);

  std::cout << vid << ", " << name << "(" << join_key << ")" << std::endl;

  user_feature vdata(join_key); 
  ASSERT_FALSE(join_key == 0); // 0 is reserved for non_join_vertex
  graph.add_vertex(vid, vdata);
  return true;
}

inline bool is_join(const graph_type::vertex_type& vertex) {
  return vertex.data().join_key != 0;
}

bool load_and_initialize_graph(graphlab::distributed_control& dc,
                               graph_type& graph,
                               const std::string& edge_dir,
                               const std::string& vertex_dir) {  
  dc.cout() << "Loading pagerank graph." << std::endl;
  graphlab::timer timer; timer.start();

  graph.load_format(edge_dir, "snap"); 
  graph.load(vertex_dir, vertex_line_parser);

  dc.cout() << "Loading pagerank graph. Finished in "
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Finalizing pagerank graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing pagerank graph. Finished in "
            << timer.current_time() << " seconds." << std::endl;

  // must call finalize before querying the graph
  dc.cout() << " #vertices: " << graph.num_vertices() 
            << " #edges:" << graph.num_edges()
            << std::endl;

  size_t isjoin = graph.map_reduce_vertices<size_t>(is_join);
  ASSERT_EQ(isjoin, graph.num_vertices());
  TOP_PAGES.lock.lock();
  TOP_PAGES.json_string = "{\n" + TOP_PAGES.json_header_string() +
    "\t\"values\": [] \n }";
  TOP_PAGES.lock.unlock();
  return true;
}


////////// Vertex Program 1 //////////////
// One pass trasformation to initialize the topic-specific jump probability
class compute_transit_prob :
  public graphlab::ivertex_program<graph_type, gather_type>,
  public graphlab::IS_POD_TYPE {

    std::vector<int> normalizer;

public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::OUT_EDGES;
    }

    // gather_type is vector<float> of length = #topics.
    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      int numdocs = edge.target().data().numdocs; 
      const std::vector<int>& topics = edge.target().data().topics; 
      std::vector<int> ret(topics);
      for (size_t i = 0; i < topics.size(); ++i) ret[i] *= numdocs;
      return ret;
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
      normalizer = total;
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::OUT_EDGES;
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      float pij = 0.0; // the marginal (over topics) jump probability from this vertex to its target
      std::vector<int>& pk_ij = edge.target().data().topics;
      int nj = edge.target().data().numdocs;
      for (size_t k = 0; k < vertex.data().topics.size(); ++k) {
        pij += vertex.data().topics[k] * (nj * pk_ij[k] / (float)normalizer[k]);
      }
      pij /= std::accumulate(vertex.data().topics.begin(), vertex.data().topics.end(), 0);
      edge.data() = pij;
    }
  };

/////////////// Vertex Program 2 //////////////////
class compute_pagerank :
  public graphlab::ivertex_program<graph_type, float>,
  public graphlab::IS_POD_TYPE {
public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::IN_EDGES;
    }

  /* Gather the weighted rank of the adjacent page   */
  float gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    return ((1.0 - RESET_PROB) * edge.data() * edge.source().data().rank);
  }

  /* Use the total rank of adjacent pages to update this page */
  void apply(icontext_type& context, vertex_type& vertex,
             const float& total) {
    const double newval = total + RESET_PROB;
    vertex.data().rank = newval;
  }

  /* The scatter edges depend on whether the pagerank has converged */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::OUT_EDGES;
  }

  /* The scatter function just signal adjacent pages */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    context.signal(edge.target());
  }
}; // end of factorized_pagerank update functor
}
