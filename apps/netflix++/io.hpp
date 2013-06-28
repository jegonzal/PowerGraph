#ifndef GRAPHLAB_NETFLIXPP_IO_HPP 
#define GRAPHLAB_NETFLIXPP_IO_HPP 

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <graphlab/util/hdfs.hpp>

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
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
  else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;
  // Parse the line
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  const bool success = qi::phrase_parse
      (line.begin(), line.end(),
       //  Begin grammar
       (
           qi::ulong_[phoenix::ref(source_id) = qi::_1] >> -qi::char_(',') >>
           qi::ulong_[phoenix::ref(target_id) = qi::_1] >>
           -(-qi::char_(',') >> qi::float_[phoenix::ref(obs) = qi::_1])
       )
       ,
       //  End grammar
       ascii::space);

  if(!success) return false;


  // map target id into a separate number space
  target_id = movieid2vid(target_id);
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(obs, role));
  return true; // successful load
} // end of graph_loader


inline bool genre_feature_loader(graph_type& graph,
                                 const std::string& filename,
                                 const std::string& line) {
  ASSERT_FALSE(line.empty());
  if (boost::starts_with(line, "#") || boost::starts_with(line, "%")) {
    logstream(LOG_INFO) << line << std::endl;
    return true;
  }
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  graph_type::vertex_id_type id(-1);
  vertex_data vdata;
  std::vector<std::string> genres;

  using phoenix::ref;
  using phoenix::push_back;
  // Parse the line
  qi::rule<std::string::const_iterator,
      std::string(), ascii::space_type> quoted_string;
  quoted_string %= qi::lexeme['"' >> +(qi::char_ - '"') >> '"'];

  bool success = qi::phrase_parse(line.begin(), line.end(), 
                                  (  qi::ulong_[ref(id) = qi::_1]
                                   >> "{" 
                                   >> (qi::omit[quoted_string] >> ":"
                                       >> quoted_string[push_back(ref(genres), qi::_1)]) % ","
                                   >> "}"), ascii::space);
  if (!success) {
    logstream(LOG_WARNING) << "Error parsing vertex data with genre features:\n" 
                           << line << std::endl;
    return false;
  }

  // map target id into a separate number space
  id = movieid2vid(id);

  for (size_t i = 0; i < genres.size(); ++i) {
    std::string name = genres[i];
    feature_table.add_feature(name);
    vdata.features[feature_table.get_key(name)] = 1;
    if (id == 1) {
      std::cout << name << "\t" << feature_table.get_key(name) << "\t" 
                << std::endl;
    }
  }

  // Create an edge and add it to the graph
  graph.add_vertex(id, vdata);
  return true; // successful load
}

inline bool topic_feature_loader(graph_type& graph,
                                 const std::string& filename,
                                 const std::string& line) {
  ASSERT_FALSE(line.empty());
  if (boost::starts_with(line, "#") || boost::starts_with(line, "%")) {
    logstream(LOG_INFO) << line << std::endl;
    return true;
  }
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  graph_type::vertex_id_type id(-1);
  vertex_data vdata;

  std::vector<int> topics;
  using phoenix::ref;
  using phoenix::push_back;
  bool success = qi::phrase_parse(line.begin(), line.end(), 
                                  (
                                      qi::ulong_[ref(id) = qi::_1]
                                      >> *(qi::int_[push_back(ref(topics), qi::_1)])
                                  )
                                  , ascii::space);
  if (!success) {
    logstream(LOG_WARNING) << "Error parsing vertex data with topic features:\n" 
                           << line << std::endl;
    return false;
  }

  // map target id into a separate number space
  id = movieid2vid(id);

  for (size_t i = 0; i < topics.size(); ++i) {
    std::string name = std::string("topic")+boost::lexical_cast<std::string>(i);
    feature_table.add_feature(name);
    vdata.features[feature_table.get_key(name)] += topics[i];
    if (id == 1) {
      std::cout << name << "\t" << feature_table.get_key(name) << "\t" 
                << topics[i] << std::endl;
    }
  }

  // Create an edge and add it to the graph
  graph.add_vertex(id, vdata);
  return true; // successful load
}


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
      std::string ret = boost::lexical_cast<std::string>(vid2movieid(vertex.id())) + " ";
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

struct prediction_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    return ""; //nop
  }
  std::string save_edge(const edge_type& edge) const {
    if(edge.data().role == edge_data::PREDICT) {
      std::stringstream strm;
      const double prediction =
        edge.source().data().factor.dot(edge.target().data().factor);
      strm << edge.source().id() << '\t';
      strm << (vid2movieid(edge.target().id())) << '\t';
      strm << prediction << '\n';
      return strm.str();
    } else return "";
  }
}; // end of prediction_saver

void save_model(distributed_control& dc, graph_type& graph) {
  std::string prefix = saveprefix 
      + "_D=" + boost::lexical_cast<std::string>(vertex_data::NLATENT)
      + "_iter=" + boost::lexical_cast<std::string>(MAX_ITER)
      + "_lambda=" + boost::lexical_cast<std::string>(LAMBDA);

  // Save model parameters 
  dc.cout() << "Save model to " << prefix << "..." << std::endl;
  graph.save(prefix + ".user", linear_model_saver_U(),
             false, true, false);
  graph.save(prefix + ".movie", linear_model_saver_V(),
             false, true, false);

  // Save hyper parameters, and global results
  if (dc.procid() == 0) {
    std::string outfile = prefix + ".info";
    std::stringstream ss;
    ss << "\"NLATENT\": " << vertex_data::NLATENT << "\n"
       << "\"iter\": " << MAX_ITER << "\n"
       << "\"lambda\": " << LAMBDA << "\n" 
       << "\"NTRAIN\": "<< NTRAIN << "\n"
       << "\"NTEST\": "<< NTEST << "\n"
       << "\"Training RMSE\": " << TRAIN_RMSE << "\n"
       << "\"Testing RMSE\": " << TEST_RMSE << "\n"; 

    if (boost::starts_with(prefix, "hdfs://")) {
      typedef graphlab::hdfs::fstream base_fstream_type;
      if(!hdfs::has_hadoop()) {
        logstream(LOG_FATAL)
          << "\n\tAttempting to save a graph to HDFS but GraphLab"
          << "\n\twas built without HDFS."
          << std::endl;
      }
      hdfs& hdfs = hdfs::get_hdfs();
      // open the stream
      base_fstream_type fout(hdfs, outfile, true);
      fout << ss.str();
      fout.close();
     } else {
      std::ofstream fout;
      fout.open(outfile.c_str());
      fout << ss.str();
      fout.close();
    }
  }
  dc.cout() << "done" << std::endl;
}
#endif
