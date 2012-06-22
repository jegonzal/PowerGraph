
/**  
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
 * This application is designed to support structured prediction on a
 * graph.
 * 
 * \todo Finish documenting
 *
 *  \author Joseph Gonzalez
 */


#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>


#include <Eigen/Dense>
#include "eigen_serialization.hpp"



#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>





typedef Eigen::VectorXd factor_type;


double SMOOTHING = 2;
double DAMPING = 0.1;
double TOLERANCE = 0.01;


struct vertex_data {
  factor_type belief;
  factor_type potential;
  void load(graphlab::iarchive& arc) { arc >> belief >> potential; }
  void save(graphlab::oarchive& arc) const { arc << belief << potential; }
};


class edge_data {
  factor_type messages_[4];
  double weight_;
  size_t message_idx(size_t source_id, size_t target_id, bool is_new) {
    return size_t(source_id < target_id)  + 2 * size_t(is_new);
  }
public:
  edge_data(const double w = 1) : weight_(w) { }
  const double& weight() const { return weight_; }
  factor_type& message(size_t source_id, size_t target_id) { 
    return messages_[message_idx(source_id, target_id, true)];
  }
  factor_type& old_message(size_t source_id, size_t target_id) { 
    return messages_[message_idx(source_id, target_id, false)];
  }
  void update_old(size_t source_id, size_t target_id) { 
    old_message(source_id, target_id) = message(source_id, target_id);
  }
  void save(graphlab::oarchive& arc) const {
    for(size_t i = 0; i < 4; ++i) arc << messages_[i];
    arc << weight_;
  }
  void load(graphlab::iarchive& arc) {
    for(size_t i = 0; i < 4; ++i) arc >> messages_[i];
    arc >> weight_;
  }
}; // End of edge data


typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



/** 
 * Belief Propagation Vertex Program
 *
 */
struct bp_vertex_program : 
  public graphlab::ivertex_program< graph_type, factor_type,
                                    graphlab::messages::sum_priority >,
  public graphlab::IS_POD_TYPE {

  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /**
   * Update the old message to be the new message and collect the
   * message value.
   */
  factor_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const {
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    edge_data& edata = edge.data();
    edata.update_old(other_vertex.id(), vertex.id());
    return edata.old_message(other_vertex.id(), vertex.id());
  }; // end of gather function

  /**
   * Multiply message product by node potential and update the belief.
   */
  void apply(icontext_type& context, vertex_type& vertex, 
             const gather_type& total) {
    vertex_data& vdata = vertex.data();
    vdata.belief = vdata.potential + total;
    // Rescale the belief
    vdata.belief.array() -= double(vdata.belief.maxCoeff());
  }; // end of apply

  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  /**
   * Compute new message value for each edge.
   */
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {  
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    edge_data& edata = edge.data();
    // construct the cavity
    factor_type cavity = vertex.data().belief - 
      edata.old_message(other_vertex.id(), vertex.id());
    // compute the new message
    factor_type& new_message = 
      edata.message(vertex.id(), other_vertex.id());
    const factor_type& old_message = 
      edata.old_message(vertex.id(), other_vertex.id());
    ASSERT_NE(&new_message, &old_message);
    convolve(cavity, new_message);
    // Renormalize
    new_message.array() -= new_message.maxCoeff();
    // Dampen message
    new_message = DAMPING * old_message + (1-DAMPING) * new_message;
    // Compute message residual
    const double residual = (new_message - old_message).cwiseAbs().sum();
    context.clear_gather_cache(other_vertex);
    // Schedule the adjacent vertex
    if(residual > TOLERANCE) context.signal(other_vertex, residual);
 }; // end of scatter

private:
  static void convolve(const factor_type& cavity, factor_type& message) {
    for(size_t i = 0; i < message.size(); ++i) {
      double value = 0;
      for(size_t j = 0; j < cavity.size(); ++j) {
        value += cavity(j) * ( i == j? 0 : -SMOOTHING ); 
      }
      message(i) = value;
    }
  } // end of convolve

  /**
   * Return the other vertex
   */
  vertex_type get_other_vertex(edge_type& edge, 
                               const vertex_type& vertex) const {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
  }; // end of other_vertex

}; // end of class bp_vertex_program








bool vertex_loader(graph_type& graph, const std::string& fname, 
                   const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  std::cout << "Reading line: " << line << std::endl;
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  graphlab::vertex_id_type vid(-1);
  std::vector<double> values;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(vid) = qi::_1] >> qi::char_("\t,") >>
      (qi::double_[phoenix::push_back(phoenix::ref(values), qi::_1)] % qi::char_("\t,") )
      )
     ,
     //  End grammar
     ascii::space); 
  if(!success) return false;

  std::cout << "value: " << vid << ": ";
  for(size_t i = 0; i < values.size(); ++i) std::cout << values[i] << ", ";
  std::cout << std::endl;

  vertex_data vdata;
  vdata.potential.resize(values.size());
  for(size_t i = 0; i < values.size(); ++i) {
    ASSERT_GT(values[i], 0);
    vdata.potential(i) = std::log(values[i]);
  }
  graph.add_vertex(vid, vdata);
  return true;
} // end of vertex_loader;



bool edge_loader(graph_type& graph, const std::string& fname, 
                 const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  std::cout << "Reading edgeline: " << line << std::endl;
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  graphlab::vertex_id_type source(-1), target(-1);
  double weight = 1;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(source) = qi::_1] >>  qi::char_(',') 
      >> qi::ulong_[phoenix::ref(target) = qi::_1] >>  
      -(qi::char_(',') >> qi::double_[phoenix::ref(weight) = qi::_1])
      )
     ,
     //  End grammar
     ascii::space); 
  if(!success) return false;

  std::cout << source << ", " << target << ", " << weight << std::endl;
  
  graph.add_edge(source, target, edge_data(weight));
  return true;
} // end of edge loader



void edge_initializer(graph_type::edge_type& edge) {
  edge_data& edata = edge.data();
  const graphlab::vertex_id_type source_id = edge.source().id();
  const size_t nsource = edge.source().data().potential.size(); 
  const graphlab::vertex_id_type target_id = edge.target().id();
  const size_t ntarget = edge.target().data().potential.size(); 
  edata.message(source_id, target_id).resize(ntarget);
  edata.old_message(source_id, target_id).resize(ntarget);
  edata.message(target_id, source_id).resize(nsource);
  edata.old_message(target_id, source_id).resize(nsource);
} // end of edge initializer



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
  clopts.attach_option("prior", &prior_dir, prior_dir,
                       "The directory containing the prior");
  clopts.add_positional("prior");
  clopts.attach_option("graph", &graph_dir, graph_dir,
                       "The directory containing the adjacency graph");
  clopts.add_positional("graph");
  clopts.attach_option("smoothing", &SMOOTHING, SMOOTHING,
                       "The amount of smoothing (larger = more)");
  clopts.attach_option("damping", &DAMPING, DAMPING,
                       "The amount of damping (0 -> no damping and 1 -> no progress)");
  clopts.attach_option("tol", &TOLERANCE, TOLERANCE,
                       "The tolerance level for convergence.");

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

  // typedef graphlab::omni_engine<bp_vertex_program> engine_type;
  // engine_type engine(dc, graph, clopts, "asynchronous");
  // engine.signal_all();
  // graphlab::timer timer;
  // engine.start();  
  // const double runtime = timer.current_time();
  //   dc.cout() 
  //   << "----------------------------------------------------------" << std::endl
  //   << "Final Runtime (seconds):   " << runtime 
  //   << std::endl
  //   << "Updates executed: " << engine.num_updates() << std::endl
  //   << "Update Rate (updates/second): " 
  //   << engine.num_updates() / runtime << std::endl;

    
  graphlab::stop_metric_server_on_eof();
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main





















