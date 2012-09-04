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


/**
 *
 * \brief This application performs MAP inference on Markov Nets 
 * provided in standard UAI file format via Dual-Decomposition. 
 *
 *
 *  \author Dhruv Batra
 */


#include "dd_grlab.hpp"





/**
 * \brief The vertex load is used by the graph loading API to parse
 * the lines of prior data in the vertex data file.
 *
 * This parser uses the boost::spirit library to parse the vertex data
 * file. As a consequence it is fairly flexible allowing both comma
 * and tab delimited files as well as vertices with different numbers
 * of states.
 */
bool vertex_loader(graph_type& graph, const std::string& fname, 
                   const std::string& line) {
  // If the line is empty simply skip it
  if(line.empty()) return true;
  // We use the boost spirit parser which requires (too) many separate
  // namespaces so to make things clear we shorten them here.
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  graphlab::vertex_id_type vid(-1);
  std::vector<double> values;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(vid) = qi::_1] >> -qi::char_(",") >>
      (qi::double_[phoenix::push_back(phoenix::ref(values), qi::_1)] % -qi::char_(",") )
      )
     ,
     //  End grammar
     ascii::space); 
  // Test to see if the boost parser was able to parse the line
  if(!success) {
    logstream(LOG_ERROR) << "Parse error in vertex prior parser." << std::endl;
    return false;
  }

  // Ensure that a prior was provided.  Technically this should not be
  // reached since the parser requires at least one prior entry
  if(values.empty()) {
    logstream(LOG_ERROR) << "Vertex has no prior." << std::endl;
    return false;
  }

  // Renormalize the vertex data. We require positive probabilities.
  double sum = 0;
  for(size_t i = 0; i < values.size(); ++i) {
    if(values[i] < 0) { 
      logstream(LOG_ERROR) << "Encountered negative probability." << std::endl;
      return false;
    }
    if(values[i] == 0) { 
      logstream(LOG_ERROR) 
        << "Zero probability assignments are not currently supported." << std::endl;
      return false;
    }
    sum += values[i]; 
  }
  ASSERT_GT(sum, 0);
  for(size_t i = 0; i < values.size(); ++i) values[i] /= sum;

  vertex_data vdata;
  vdata.potential.resize(values.size());
  for(size_t i = 0; i < values.size(); ++i) {
    ASSERT_GT(values[i], 0);
    vdata.potential(i) = std::log(values[i]);
  }
  graph.add_vertex(vid, vdata);
  return true;
} // end of vertex_loader;


/**
 * \brief The edge data loader is used by the GraphLab graph loading
 * API to parse lines in the edge data file. 
 */
bool edge_loader(graph_type& graph, const std::string& fname, 
                 const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  graphlab::vertex_id_type source(-1), target(-1);
  double weight = 1;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(source) = qi::_1] >>  -qi::char_(',') 
      >> qi::ulong_[phoenix::ref(target) = qi::_1] >>  
      -(-qi::char_(',') >> qi::double_[phoenix::ref(weight) = qi::_1])
      )
     ,
     //  End grammar
     ascii::space); 
  if(!success) return false;  
  graph.add_edge(source, target, edge_data(weight));
  return true;
} // end of edge loader



/**
 * \brief The edge initializer is used to allocate the messages along
 * each edge based on the number of states of the source and target
 * vertex.
 */
void edge_initializer(graph_type::edge_type& edge) {
  edge_data& edata = edge.data();
  const graphlab::vertex_id_type source_id = edge.source().id();
  const size_t nsource = edge.source().data().potential.size(); 
  const graphlab::vertex_id_type target_id = edge.target().id();
  const size_t ntarget = edge.target().data().potential.size();
  edata.initialize(source_id, nsource, target_id, ntarget);
} // end of edge initializer




/**
 * \brief The belief prediction saver is used to save the belief
 * predictions for each vertex.
 */
struct belief_prediction_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    std::stringstream strm;
    strm << vertex.id() << '\t';
    factor_type pred = vertex.data().belief;
    double sum = 0;
    for(int i = 0; i < pred.size(); ++i) 
      sum += (pred(i) = std::exp(pred(i)));
    pred.array() /= sum;
    for(int i = 0; i < pred.size(); ++i) 
      strm << pred(i) << (i+1 < pred.size()? '\t' : '\n');
    return strm.str();
  }
  std::string save_edge(const edge_type& edge) const {
    return ""; // nop
  }
}; // end of belief_prediction_saver


/**
 * \brief The MAP prediction saver is used to save the map estimated
 * for each vertex.  The MAP estimate is the most likely assignment
 */
struct map_prediction_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    std::stringstream strm;
    size_t prediction = 0;
    vertex.data().belief.maxCoeff(&prediction);
    strm << vertex.id() << '\t' << prediction << '\n';
    return strm.str();
  }
  std::string save_edge(const edge_type& edge) const {
    return ""; // nop
  }
}; // end of map prediction_saver





int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  // Parse command line options -----------------------------------------------
  // \todo update description string
  const std::string description = "Structure prediction solver";
  graphlab::command_line_options clopts(description);
  std::string prior_dir; 
  std::string graph_dir;
  std::string output_dir = "pred";
  std::string exec_type = "async";
  bool map = false;
  clopts.attach_option("prior", prior_dir,
                       "The directory containing the prior");
  clopts.add_positional("prior");
  clopts.attach_option("graph", graph_dir,
                       "The directory containing the adjacency graph");
  clopts.add_positional("graph");
  clopts.attach_option("output", output_dir,
                       "The directory in which to save the predictions");
  clopts.add_positional("output");
  clopts.attach_option("smoothing", SMOOTHING,
                       "The amount of smoothing (larger = more)");
  clopts.attach_option("damping", DAMPING,
                       "The amount of damping (0 -> no damping and 1 -> no progress)");
  clopts.attach_option("tol", TOLERANCE,
                       "The tolerance level for convergence.");
  clopts.attach_option("map", map,
                       "Return maximizing assignment instead of the posterior distribution.");
  clopts.attach_option("engine", exec_type,
                       "The type of engine to use {async, sync}.");
  if(!clopts.parse(argc, argv)) {
    graphlab::mpi_tools::finalize();
    return clopts.is_set("help")? EXIT_SUCCESS : EXIT_FAILURE;
  }

  if(prior_dir.empty()) {
    logstream(LOG_ERROR) << "No prior was provided." << std::endl;
    return EXIT_FAILURE;
  }

  if(graph_dir.empty()) {
    logstream(LOG_ERROR) << "No graph was provided." << std::endl;
    return EXIT_FAILURE;
  }

  // Start the webserver
  graphlab::launch_metric_server();

  ///! load the graph
  graph_type graph(dc, clopts);  


  ///! load the graph
  graph.load(prior_dir, vertex_loader);
  graph.load(graph_dir, edge_loader);
  graph.finalize();
  graph.transform_edges(edge_initializer);

  typedef graphlab::omni_engine<bp_vertex_program> engine_type;
  engine_type engine(dc, graph, exec_type, clopts);
  engine.signal_all();
  graphlab::timer timer;
  engine.start();  
  const double runtime = timer.current_time();
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;
    
    
  std::cout << "Saving predictions" << std::endl;
  const bool gzip_output = false;
  const bool save_vertices = true;
  const bool save_edges = false;
  const size_t threads_per_machine = 2;
  if(map) {
    graph.save(output_dir, map_prediction_saver(),
               gzip_output, save_vertices, 
               save_edges, threads_per_machine);
  } else { 
    graph.save(output_dir, belief_prediction_saver(),
               gzip_output, save_vertices, 
               save_edges, threads_per_machine);
  }


    
  //  graphlab::stop_metric_server_on_eof();
  graphlab::stop_metric_server();
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main





















