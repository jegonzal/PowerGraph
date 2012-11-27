#include <graphlab.hpp>
#include <vector>
#include <algorithm>
#include <graphlab/macros_def.hpp>

namespace pagerank {
  typedef std::vector<float> gather_type;
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
int DEFAULT_NDOCS = 3;
int DEFAULT_TOPICVAL = 1;
int NPAGES = 0;
int NTOPICS = 50;
size_t TOPK = 50;
bool HAS_DOC_COUNT = false;
bool JOIN_ON_ID = true;
std::vector<float> w_personal; // personalization weights for topics

struct top_pages_type {
  graphlab::mutex lock;
  std::string json_string;
  top_pages_type() : 
    json_string("{\n" + json_header_string() + "\tvalues: [] \n }") { }
  inline std::string json_header_string() const {
    return
      "\t\"npages\": " + graphlab::tostr(NPAGES) + ",\n" + 
      "\t\"ntopics\": " + graphlab::tostr(NTOPICS) + ",\n" + 
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


std::pair<std::string, std::string>
weight_update_callback(std::map<std::string, std::string>& varmap) {
  std::string c = varmap["topic"];
  size_t topick = atoi(c.c_str());
  std::pair<std::string, std::string> pair("text/html","");
  if (topick < w_personal.size()) {
    w_personal[topick] *= 1.2;
    double s = 0;
    for (size_t i = 0;i < w_personal.size(); ++i) s += w_personal[i];
    for (size_t i = 0;i < w_personal.size(); ++i) w_personal[i] /= s;
    pair.second = "1";
  }
  else {
    pair.second = "0";
  }
  // add a boost to the topic 
  return pair;
}



// Graph Types
// ============================================================================
struct user_feature {
   size_t join_key; // hash value of the user name, used as vertex join key
   float rank;
   int numdocs;
   std::vector<float> topics;
   std::vector<float> out_normalizer;
   user_feature() : join_key(-1), rank(1.0),
                              numdocs(DEFAULT_NDOCS),
                              topics(std::vector<float>(NTOPICS, 1.0/NTOPICS)),
                              out_normalizer(topics) { }
   user_feature(size_t join_key) : join_key(join_key), rank(1.0),
                              numdocs(DEFAULT_NDOCS),
                              topics(std::vector<float>(NTOPICS, 1.0/NTOPICS)),
                              out_normalizer(topics) { }
   std::string topics_tostr() const {
     std::stringstream strm;
     for (size_t i = 0; i < topics.size(); ++i) {
       strm << topics[i];
       if (i < topics.size()-1)
         strm << ", ";
     }
     return strm.str();
   } 
  void save(graphlab::oarchive& oarc) const {
    oarc << join_key << rank << numdocs << topics << out_normalizer;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> join_key >> rank >> numdocs >> topics >> out_normalizer;
  }
};

// The vertex data is its pagerank 
typedef user_feature vertex_data_type;
// The edge data is its the marginal transition probability 
typedef float edge_data_type;
// The gather type (for compute_transit vertex program) is the user topic vector. 
typedef std::vector<float> gather_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

// void init_vertex(graph_type::vertex_type& vertex) { 
//   vertex.data().numdocs = DEFAULT_NDOCS; 
//   std::fill(vertex.data().topics.begin(), vertex.data().topics.end(), 1);
//   vertex.data().rank = 1.0;
// }

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
  int numdocs = DEFAULT_NDOCS;

  if (HAS_DOC_COUNT) {
    const bool success = qi::phrase_parse
      (textline.begin(), textline.end(),       
       //  Begin grammar
       (
        qi::ulong_[phoenix::ref(vid) = qi::_1] >> -qi::char_(',') >>
        qi::int_[phoenix::ref(numdocs) = qi::_1] >> -qi::char_(',') >>
        +qi::char_[phoenix::ref(name) = qi::_1]
        )
       ,
       //  End grammar
       ascii::space); 
    if(!success) return false;  
  } else {
    const bool success = qi::phrase_parse
      (textline.begin(), textline.end(),       
       //  Begin grammar
       (
        qi::ulong_[phoenix::ref(vid) = qi::_1] >> -qi::char_(',') >>
        +qi::char_[phoenix::ref(name) = qi::_1]
        )
       ,
       //  End grammar
       ascii::space); 
    if(!success) return false;  
  }

  size_t join_key = boost::hash<std::string>()(name);

  // std::cout << vid << ", " << name << "(" << join_key << ")" << std::endl;

  user_feature vdata(join_key); 
  vdata.numdocs = numdocs;
  ASSERT_FALSE(join_key == size_t(-1)); // -1 is reserved for non_join_vertex
  graph.add_vertex(vid, vdata);
  return true;
}

/* Global thread func for join vertices and recompute transit probs */
size_t left_emit_key (const graph_type::vertex_type& vertex) {
  return JOIN_ON_ID ? vertex.id() : vertex.data().join_key;
}
inline bool is_join(const graph_type::vertex_type& vertex) {
  return left_emit_key(vertex) != size_t(-1);
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
  NPAGES = graph.num_vertices();
  TOP_PAGES.lock.lock();
  TOP_PAGES.json_string = "{\n" + TOP_PAGES.json_header_string() +
    "\t\"values\": [] \n }";
  TOP_PAGES.lock.unlock();
  return true;
}


////////// Vertex Program 1 //////////////
// One pass trasformation to initialize the topic-specific jump probability
class compute_transit_prob :
  public graphlab::ivertex_program<graph_type, gather_type>
  {

    std::vector<float> normalizer;

public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::ALL_EDGES;
    }

    void save(graphlab::oarchive& arc) const {
      arc << normalizer;
    }

    void load(graphlab::iarchive& arc) {
      arc >> normalizer; 
    }


    // gather_type is vector<float> of length = #topics.
    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      vertex_type other = edge.target().id() == vertex.id() ? edge.source() : edge.target();
      int numdocs = other.data().numdocs; 
      const std::vector<float>& topics = other.data().topics; 
      std::vector<float> ret(topics);
      return ret;
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
      vertex.data().out_normalizer = total;
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::NO_EDGES;
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
      vertex_type other = edge.target().id() == vertex.id() ? edge.source() : edge.target();
      double pij = 0.0; // the marginal (over topics) jump probability from this vertex to its target
      std::vector<float>& pk_ij = other.data().topics;
      int numdocs = other.data().numdocs; 
      ASSERT_EQ(pk_ij.size(), w_personal.size());
      for (size_t k = 0; k < vertex.data().topics.size(); ++k) {
        if (normalizer[k] > 0) 
        pij +=  (w_personal[k] * pk_ij[k] / (float)normalizer[k]);  // p(i -> j | k)
        // pij +=  (1.0/NTOPICS)*(pk_ij[k] / ((float)normalizer[k]) );  // p(i -> j | k)
      }
      edge.data() = (float)pij; 
    }
  }; // end of compute_transit_prob

/////////////// Vertex Program 2 //////////////////
class compute_pagerank :
  public graphlab::ivertex_program<graph_type, float>,
  public graphlab::IS_POD_TYPE {
public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
      return graphlab::ALL_EDGES;
    }

  /* Gather the weighted rank of the adjacent page   */
  float gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    vertex_type other = edge.target().id() == vertex.id() ? edge.source() : edge.target();

    float pij = 0;
    for (size_t k = 0; k < vertex.data().topics.size(); ++k) {
        if (other.data().out_normalizer[k] > 0) {
          pij +=  (w_personal[k] * vertex.data().topics[k] / (float)other.data().out_normalizer[k]);  
        }
    }

    return ((1.0 - RESET_PROB) * pij * other.data().rank);
    // return ((1.0 - RESET_PROB) * edge.source().data().rank / edge.source().num_out_edges());
  }

  /* Use the total rank of adjacent pages to update this page */
  void apply(icontext_type& context, vertex_type& vertex,
             const float& total) {
    const double newval = total + RESET_PROB;
    vertex.data().rank = std::isnan(newval) ? 1 : newval;
    context.signal(vertex);
  }

  /* The scatter edges depend on whether the pagerank has converged */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }

  /* The scatter function just signal adjacent pages */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    //context.signal(edge.target());
  }
}; // end of factorized_pagerank update functor

typedef compute_pagerank::icontext_type icontext_type;
// ========================================================
// Aggregators
struct vinfo_type {
  graphlab::vertex_id_type vid;
  user_feature features;
  void save(graphlab::oarchive& arc) const {
    arc << vid << features;
  }
  void load(graphlab::iarchive& arc) {
    arc >> vid >> features; 
  }
  bool operator<(const vinfo_type& other) const {
      if (features.rank < other.features.rank) 
        return true;
      else if (features.rank == other.features.rank)
        return (vid < other.vid);
      else
        return false;
  }
};

/** Periodically compute top k rank pages. */
class topk_aggregator {
private:
  std::set<vinfo_type> top_pages;
public:
  topk_aggregator() {}
  void save(graphlab::oarchive& arc) const { arc << top_pages; }
  void load(graphlab::iarchive& arc) { arc >> top_pages; }

  topk_aggregator& operator+=(const topk_aggregator& other) {
    if(other.top_pages.empty()) return *this;
      // Merge the topk
      top_pages.insert(other.top_pages.begin(),
                          other.top_pages.end());
      // Remove excess elements
      while(top_pages.size() > TOPK)
        top_pages.erase(top_pages.begin());
    return *this;
  } // end of operator +=

  static topk_aggregator map(icontext_type& context,
                             const graph_type::vertex_type& vertex) {
    topk_aggregator ret_value;
    vinfo_type vinfo;
    vinfo.vid = vertex.id();
    vinfo.features = vertex.data();
    ret_value.top_pages.insert(vinfo);
    return ret_value;
  } // end of map function


  static void finalize(icontext_type& context,
                       const topk_aggregator& total) {
    if(context.procid() != 0) return;
    std::string json = "{\n"+ TOP_PAGES.json_header_string() +
      "\t\"values\": [\n";
      size_t counter = 0;
      std::set<vinfo_type>::iterator iter = total.top_pages.begin();
      while(iter != total.top_pages.end())  {
        json += "\t[\"" + graphlab::tostr(iter->vid) + "\", " +
          graphlab::tostr(iter->features.rank) + ", " +
          graphlab::tostr(iter->features.numdocs) + ", " +
          "[" + iter->features.topics_tostr() + "]"+
          "]";
        if(++counter < total.top_pages.size()) json += ", ";
        json += '\n';
        std::cout << graphlab::tostr(iter->vid)
                  << ": " << iter->features.rank
                  << " (" << iter->features.topics_tostr() << ")" << "\n";
        ++iter;
      }
      json += '\n';
      std::cout << std::endl;
    json += "]}";
    // Post the change to the global variable
    TOP_PAGES.lock.lock();
    TOP_PAGES.json_string.swap(json);
    TOP_PAGES.lock.unlock();
  } // end of finalize
}; // end of topk_aggregator
} // end of namespace pagerank
#include <graphlab/macros_undef.hpp>
