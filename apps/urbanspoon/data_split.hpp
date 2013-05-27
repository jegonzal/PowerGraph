#ifndef DATA_SPLIT_HPP 
#define DATA_SPLIT_HPP 
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

using namespace graphlab;
//add 2 to negative node id to prevent -0 and -1 which arenot allowed
const int SAFE_NEG_OFFSET = 2;

/**
 * \ingroup toolkit_matrix_factorization
 *
 * \brief the vertex data type which contains the latent factor.
 *
 * Each row and each column in the matrix corresponds to a different
 * vertex in the ALS graph.  Associated with each vertex is a factor
 * (vector) of latent parameters that represent that vertex.  The goal
 * of the ALS algorithm is to find the values for these latent
 * parameters such that the non-zero entries in the matrix can be
 * predicted by taking the dot product of the row and column factors.
 */
struct vertex_data : public IS_POD_TYPE{ }; // end of vertex data

size_t hash_value(vertex_data const& v) {
  return 0;
}

/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data also stores the most recent error estimate.
 */
struct edge_data {
  /** \brief the observed value for the edge */
  std::vector<std::string> fields;

  boost::gregorian::date d;

  /** \brief basic initialization */
  edge_data() { }
  edge_data(std::vector<std::string>& fields) : fields(fields) { 
    d = boost::gregorian::from_simple_string(fields[1]);
  }

  void load(iarchive& iarc) {
    iarc >> fields;
    d = boost::gregorian::from_simple_string(fields[1]);
  }
  void save(oarchive& oarc) const {
    oarc << fields;
  }
}; // end of edge data

/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */
typedef distributed_graph<vertex_data, edge_data> graph_type;

/**
 * \brief The engine type used by the ALS matrix factorization
 * algorithm.
 *
 * The ALS matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef gl3engine<graph_type> engine_type;


////////////////////////////  INPUT UTILITIES    ////////////////////////
/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph,
                         const std::string& filename,
                         const std::string& line) {
  ASSERT_FALSE(line.empty());
  if (boost::starts_with(line, "#") || boost::starts_with(line, "%")) {
    logstream(LOG_INFO) << line << std::endl;
    return true;
  }

  // Define parsers
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  qi::rule<std::string::const_iterator, unsigned long(), ascii::space_type> quoted_long;
  qi::rule<std::string::const_iterator, float(), ascii::space_type> quoted_float;
  qi::rule<std::string::const_iterator, std::string(), ascii::space_type> quoted_string;
  quoted_long %= qi::lexeme['"' >> qi::ulong_ >> '"'];
  quoted_float %= qi::lexeme['"' >> qi::float_ >> '"'];
  quoted_string %= qi::lexeme['"' >> +(qi::char_ - '"') >> '"'];

  // Parse the line
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  std::vector<std::string> fields;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),
     //  Begin grammar
     (
      quoted_long[phoenix::ref(source_id) = qi::_1] >> -qi::char_(',') >>   // userid
      quoted_long[phoenix::ref(target_id) = qi::_1] >> -qi::char_(',') >>  // restid
      (quoted_string[phoenix::push_back(phoenix::ref(fields), qi::_1)] % ',')
      ),
     //  End grammar
     ascii::space);

  if(!success) return false;

  // map target id into a separate number space
  target_id = -(vertex_id_type(target_id + SAFE_NEG_OFFSET));
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(fields));
  return true; // successful load
} // end of graph_loader

struct graph_edge_writer {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    return "";
  }
  std::string save_edge(const edge_type& edge) const {
    std::stringstream ss;
    ss << "\"" << edge.source().id() << "\",";
    ss << "\"" << (-edge.target().id()  - SAFE_NEG_OFFSET) << "\",";
    const std::vector<std::string>& fields = edge.data().fields;
    for (size_t i = 0; i < fields.size(); ++i) {
      ss << "\"" << fields[i] << "\"";
      if (i < fields.size()-1) {
        ss << ",";
      }
    }
    ss << "\n";
    return ss.str();
  }
};
#endif
