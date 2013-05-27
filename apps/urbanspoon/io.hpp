#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

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

  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::find_first(filename,"test").size()) role = edge_data::VALIDATE;

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
  std::string date_str;
  float obs(0);
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),
     //  Begin grammar
     (
      quoted_long[phoenix::ref(source_id) = qi::_1] >> -qi::char_(',') >>   // userid
      quoted_long[phoenix::ref(target_id) = qi::_1] >> -qi::char_(',') >>  // restid
      quoted_float[phoenix::ref(obs) = qi::_1] >> -qi::char_(',') >>
      quoted_string[phoenix::ref(date_str) = qi::_1] >> -qi::char_(',') >>  // date
      quoted_long // num comments
      ),
     //  End grammar
     ascii::space);

  if(!success) return false;

  // map target id into a separate number space
  target_id = -(vertex_id_type(target_id + SAFE_NEG_OFFSET));

  edge_data edata(obs, role);

  if (USE_DATE_FEATURE) {
    boost::hash<std::string> hash;
    std::vector<std::string> dt_splits;
    boost::algorithm::split(dt_splits,
                            date_str,
                            boost::is_any_of("-:"),
                            boost::token_compress_on); 
    // year feature
    size_t year_key = hash("year_"+dt_splits[0]);
    if (!feature_weights.count(year_key)) {
      std::cout << "add feature: " << "year_"+dt_splits[0] << std::endl;;
    }
    edata.features[year_key] = 1;
    feature_weights[year_key] = 0;
    // month feature
    size_t month_key = hash("month_"+dt_splits[1]);
    if (!feature_weights.count(month_key)) {
      std::cout << "add feature: " << "month_"+dt_splits[1] << std::endl;
    }
    edata.features[month_key] = 1;
    feature_weights[month_key] = 0;
  }
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edata);

  return true; // successful load
} // end of graph_loader

/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool rest_feature_loader(graph_type& graph,
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
  graph_type::vertex_id_type vid(-1);
  vertex_data vdata;

  using phoenix::push_back;
  using phoenix::ref;
  std::vector<std::string> value_strs;
  bool success = qi::phrase_parse
    (line.begin(), line.end(),
     //  Begin grammar
     (
      quoted_long[phoenix::ref(vid) = qi::_1] >> qi::char_(',') // userid 
      >> (quoted_string[push_back(ref(value_strs), qi::_1)] 
          | qi::lit("\\N")[push_back(ref(value_strs), "NA")]) % ","
      ),
     //  End grammar
     ascii::space);

  size_t nfields = 15;
  size_t nfeatures = 9;
  int offset = nfields - nfeatures;
  if(!success || value_strs.size() != nfields) return false;
  const char *feature_names_arr[] = {"rprice", "rnupvotes", "rnvotes", "rncomments", "rnphotos", "rnreviews", "rnblogs", "rcusine", "rtags"};

  std::vector<std::string> feature_names(feature_names_arr, feature_names_arr+nfeatures);
  boost::hash<std::string> hash;

  if (USE_NUMERIC_FEATURE) {
    // numerical feature
    for (size_t i = 0; i <  (nfeatures-2); ++i) { 
      size_t key = hash(feature_names[i]);
      if (!feature_weights.count(key)) {
        std::cout << "add feature: " << feature_names[i] << std::endl;;
      }
      if (value_strs[i+offset] != "NA") {
        double val = boost::lexical_cast<double>(value_strs[i+offset]);
        vdata.features[key] = val;
        feature_latent[key] = vec_type();
        feature_weights[key] = 0;
        if (vid == 1) {
          std::cout << feature_names[i] << ": " << val << std::endl;
        }
      }
    }
  }

  if (USE_CATEGORY_FEATURE) {
    // cusine feature
    if (value_strs[nfields-2] != "NA") {
      std::vector<std::string> splits;
      boost::algorithm::split(splits,
                              value_strs[nfields-2],
                              boost::is_any_of(","),
                              boost::token_compress_on); 
      foreach (std::string& token, splits) {
        size_t key = hash("cusine_"+token);
        if (!feature_weights.count(key)) {
          std::cout << "add feature: " << "cusion_"+token << std::endl;;
        }
        vdata.features[key] = 1;
        feature_latent[key] = vec_type();
        feature_weights[key] = 0;

        if (vid == 1) {
          std::cout << token << ", ";
        }
      }
      if (vid == 1) {
        std::cout << std::endl;
      }
    }

    // tag feature
    if (value_strs[nfields-1] != "NA") {
      std::vector<std::string> splits;
      boost::algorithm::split(splits,
                              value_strs[nfields-1],
                              boost::is_any_of(","),
                              boost::token_compress_on); 
      foreach (std::string& token, splits) {
        size_t key = hash("tag_"+token);
        if (!feature_weights.count(key)) {
          std::cout << "add feature: " << "tag_"+token << std::endl;;
        }
        vdata.features[key] = 1;
        feature_latent[key] = vec_type();
        feature_weights[key] = 0;
        if (vid == 1) {
          std::cout << token << ", ";
        }
      }
      if (vid == 1) {
        std::cout << std::endl;
      }
    }
  } //  end of load category feature

  vid = -(vertex_id_type(vid + SAFE_NEG_OFFSET));
  // Create an edge and add it to the graph
  graph.add_vertex(vid, vdata);
  return true; // successful load
} // end of rest_feature_loader

/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool user_feature_loader(graph_type& graph,
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
  graph_type::vertex_id_type vid(-1);
  vertex_data vdata;

  using phoenix::push_back;
  using phoenix::ref;
  std::vector<std::string> value_strs;
  bool success = qi::phrase_parse
    (line.begin(), line.end(),
     //  Begin grammar
     (
      quoted_long[phoenix::ref(vid) = qi::_1] >> qi::char_(',') // userid 
      >> (quoted_string[push_back(ref(value_strs), qi::_1)] 
          | qi::lit("\\N")[push_back(ref(value_strs), "NA")]) % ","
      ),
     //  End grammar
     ascii::space);


  size_t total_fields = 6;
  size_t nfeatures = 4;
  int offset = total_fields - nfeatures;
  if(!success || value_strs.size() != total_fields) return false;
  const char *feature_names_arr[] = {"u_upvotes", "u_nvotes", "u_ncomments", "u_nphotos"};
  std::vector<std::string> feature_names(feature_names_arr, feature_names_arr+nfeatures);
 
  boost::hash<std::string> hash;
  if (USE_NUMERIC_FEATURE) {
    for (size_t i = 0; i <  nfeatures; ++i) { 
      size_t key = hash(feature_names[i]);
      if (!feature_weights.count(key)) {
        std::cout << "add feature: " << feature_names[i] << std::endl;;
      }
      if (value_strs[i+offset] != "NA") {
        double val = boost::lexical_cast<double>(value_strs[i+offset]);
        vdata.features[key] = val;
        feature_latent[key] = vec_type();
        feature_weights[key] = 0;
        if (vid == 1) {
          std::cout << feature_names[i] << ": " << val << std::endl;
        }
      }
    }
  }
  // Create an edge and add it to the graph
  graph.add_vertex(vid, vdata);
  return true; // successful load
}

// /**
//  * \brief Load metadata for visualization 
//  */
// inline void load_meta() {
// }

////////////////////////////  OUTPUT UTILITIES    ////////////////////////
struct linear_model_saver_U {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
  */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() > 0){
      std::string ret = boost::lexical_cast<std::string>(vertex.id()) + " ";
      for (uint i=0; i< vertex_data::NLATENT; i++)
        ret += boost::lexical_cast<std::string>(vertex.data().factor[i]) + " ";
        ret += "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
};

struct linear_model_saver_V {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
  */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() == 0){
      std::string ret = boost::lexical_cast<std::string>(-vertex.id()-SAFE_NEG_OFFSET) + " ";
      for (uint i=0; i< vertex_data::NLATENT; i++)
        ret += boost::lexical_cast<std::string>(vertex.data().factor[i]) + " ";
        ret += "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
};

// struct prediction_saver {
//   typedef graph_type::vertex_type vertex_type;
//   typedef graph_type::edge_type   edge_type;
//   std::string save_vertex(const vertex_type& vertex) const {
//     return ""; //nop
//   }
//   std::string save_edge(const edge_type& edge) const {
//     if(edge.data().role == edge_data::PREDICT) {
//       std::stringstream strm;
//       const double prediction =
//         edge.source().data().factor.dot(edge.target().data().factor);
//       strm << edge.source().id() << '\t';
//       strm << (-edge.target().id() - SAFE_NEG_OFFSET) << '\t';
//       strm << prediction << '\n';
//       return strm.str();
//     } else return "";
//   }
// }; // end of prediction_saver

//inline bool genre_feature_loader(graph_type& graph,
//                          const std::string& filename,
//                          const std::string& line,
//                          bool is_rvertex = false) {
//   ASSERT_FALSE(line.empty());
//   if (boost::starts_with(line, "#") || boost::starts_with(line, "%")) {
//     logstream(LOG_INFO) << line << std::endl;
//     return true;
//   }
//   namespace qi = boost::spirit::qi;
//   namespace ascii = boost::spirit::ascii;
//   namespace phoenix = boost::phoenix;
// 
//   graph_type::vertex_id_type id(-1);
//   vertex_data vdata;
//   std::vector<std::string> genres;
// 
//   using phoenix::ref;
//   using phoenix::push_back;
//   // Parse the line
//   qi::rule<std::string::const_iterator,
//       std::string(), ascii::space_type> quoted_string;
//   quoted_string %= qi::lexeme['"' >> +(qi::char_ - '"') >> '"'];
// 
//   bool success = qi::phrase_parse(line.begin(), line.end(), 
//                    (
//                        qi::ulong_[ref(id) = qi::_1]
//                        >> "{" 
//                        >> (qi::omit[quoted_string] >> ":"
//                            >> quoted_string[push_back(ref(genres), qi::_1)]) % ","
//                        >> "}"), ascii::space);
//   
//   if (!success) {
//     logstream(LOG_WARNING) << "Error parsing vertex data with genre features:\n" 
//                            << line << std::endl;
//     return false;
//   }
// 
//   // map target id into a separate number space
//   if (is_rvertex) {
//     id = -(vertex_id_type(id + SAFE_NEG_OFFSET));
//   }
// 
//   for (size_t i = 0; i < genres.size(); ++i) {
//     vdata.features[boost::hash<std::string>()(genres[i])] = 1;
//   }
// 
//   // Create an edge and add it to the graph
//   graph.add_vertex(id, vdata);
//   return true; // successful load
// }
// 
// inline bool topic_feature_loader(graph_type& graph,
//                          const std::string& filename,
//                          const std::string& line,
//                          bool is_rvertex = false) {
//   ASSERT_FALSE(line.empty());
//   if (boost::starts_with(line, "#") || boost::starts_with(line, "%")) {
//     logstream(LOG_INFO) << line << std::endl;
//     return true;
//   }
//   namespace qi = boost::spirit::qi;
//   namespace ascii = boost::spirit::ascii;
//   namespace phoenix = boost::phoenix;
// 
//   graph_type::vertex_id_type id(-1);
//   vertex_data vdata;
// 
//   std::vector<int> topics;
//   using phoenix::ref;
//   using phoenix::push_back;
//   bool success = qi::phrase_parse(line.begin(), line.end(), 
//                    (
//                     qi::ulong_[ref(id) = qi::_1]
//                     >> *(qi::int_[push_back(ref(topics), qi::_1)])
//                     )
//                    , ascii::space);
//   if (!success) {
//     logstream(LOG_WARNING) << "Error parsing vertex data with topic features:\n" 
//                            << line << std::endl;
//     return false;
//   }
// 
//   // map target id into a separate number space
//   if (is_rvertex) {
//     id = -(vertex_id_type(id + SAFE_NEG_OFFSET));
//   }
// 
//   for (size_t i = 0; i < topics.size(); ++i) {
//     vdata.features[boost::hash<std::string>()(std::string("topic")+boost::lexical_cast<std::string>(i))] += topics[i];
//   }
//   // Create an edge and add it to the graph
//   graph.add_vertex(id, vdata);
//   return true; // successful load
// }
