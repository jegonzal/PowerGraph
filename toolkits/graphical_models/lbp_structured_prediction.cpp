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
 * \brief This application used for structured prediction on a graph.
 * For example, suppose you want to model the interests of users in a
 * social network.  
 *
 * Overview and Usage
 * ======================
 *
 * For simplicity lets suppose you want to know the users interest in
 * the categories movies, sports, and music.  After analyzing each
 * users profile you might be able to estimate a crude distribution
 * over her interests.  However you would like to leverage similarity
 * among friends to improve your estimates.  This application is
 * designed to do exactly that.
 *
 * As an input you provide two folders (or files) the first contains
 * the prior probabilities for each vertex in the form:
 *
 *   <vertexId> \t <Pr Category1> \t <Pr Category2> ... \n
 *
 * For example:
 *   
 *   1    0.2   0.2   0.6
 *   2    0.3   0.6   0.1
 *   3    0.3   0.3   0.4
 *           ... 
 * 
 * The second folder contains the graph structure in the form:
 *
 *   <sourceId> \t <targetId> \t [Optional Weight]
 *
 * For example:
 *
 *   1   2
 *   1   3  1.7
 *   3   2  0.3
 *
 * The default weight value is 1 (times the smoothing parameter passed
 * in as a command line argument).  Larger weight values imply
 * stronger relationships.  A negative weight implies a "repulsive"
 * relationship in which neighboring vertices would like to have
 * different assignments.
 *
 * We have provided a synthetic data generator which creates a
 * synthetic dataset for an simulated image denoising task. See the
 * synthetic_image_data application for details.
 *
 * As output the application produces another set files with a format
 * identical to the vertex prior file with each weight (probability)
 * corresponding to the posterior predictions.
 *
 * Technical Explanation
 * ========================
 *
 * This application creates a pair-wise Markov Random Field with
 * Ising-Potts edge factors and then uses residual loopy belief
 * propagation to compute posterior belief estimates for each vertex.
 *
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




/**
 * \brief Eigen library vectors are used to store factor in _LOG
 * SPACE_.
 */
typedef Eigen::VectorXd factor_type;

/**
 * \brief The Ising smoothing parameter which controls the coupling
 * between adjacent predictions in the graph.  Larger values imply
 * greater smoothing (stronger coupling). 
 *
 * \code
 * edge_factor(xi, xj) = exp( (xi == xj)? 0 : -SMOOTHING * edge_weight ); 
 * \endcode
 *
 * Not that the default edge weight is 1 however the graph file can
 * contain an additional edge weight column which allows per edge
 * control of the smoothing parameter.
 *
 * This parameter is set as a command line argument.
 */
double SMOOTHING = 2;

/**
 * \brief The Damping parameter which helps ensure stable convergence.
 * Larger damping values lead to slower but more stable convergence.
 *
 * Currently damping is implemented in log-space in the following
 * equation:
 *
 * \code
 * log(new_message) = DAMPING * log(old_message) + 
 *                         (1-DAMPING) * log(new_message);
 * \endcode
 *
 * This parameter is set as a command line argument.
 */
double DAMPING = 0.1;

/**
 * \brief The convergence threshold for each message.  Smaller values
 * imply tighter convergence but slower execution.
 *
 *
 * The algorithm convergence when:
 *   
 * \code
 * sum(abs(log(old_message) - log(new_message))) < TOLERANCE
 * \endcode
 *
 * The parameter is set as a command line argument
 */
double TOLERANCE = 0.01;


/**
 * \brief The vertex data contains the vertex potential as well as the
 * current belief estimate and represents a random variable in the
 * Markov Random Field.
 *
 * The vertex potential represents the prior and is obtained from the
 * vertex prior file (stored in log form).
 *
 * The belief represents the current posterior estimate.
 */
struct vertex_data {
  factor_type belief;
  factor_type potential;
  void load(graphlab::iarchive& arc) { arc >> belief >> potential; }
  void save(graphlab::oarchive& arc) const { arc << belief << potential; }
}; // end of vertex_data


/**
 * \brief The edge data represents an edge in the Markov Random Field
 * and contains the loopy belief propagation message in both
 * directions along that edge as well as the old message in each
 * direction.  In addition each edge contains the weight parameter
 * used to set edge specific smoothing (default value is 1).
 */
class edge_data {
  /**
   * \brief We store old and new messages in both directions as an
   * array of messages.  The particular message index is then computed
   * using the \ref message_idx function.
   */
  factor_type messages_[4];
  /**
   * \brief The weight associated with the edge (used to scale the
   * smoothing parameter)
   */
  double weight_;
  /**
   * \brief The function used to compute the message index in the edge
   * message array.
   */
  size_t message_idx(size_t source_id, size_t target_id, bool is_new) {
    return size_t(source_id < target_id)  + 2 * size_t(is_new);
  }

public:

  edge_data(const double w = 1) : weight_(w) { }
  const double& weight() const { return weight_; }

  /**
   * \brief Get the new message value from source_id to target_id
   */
  factor_type& message(size_t source_id, size_t target_id) { 
    return messages_[message_idx(source_id, target_id, true)];
  }
  /**
   * \brief Get the old message value from source_id to target_id
   */
  factor_type& old_message(size_t source_id, size_t target_id) { 
     return messages_[message_idx(source_id, target_id, false)];
  }

  /**
   * \brief Set the old message value equal to the new message value
   */
  void update_old(size_t source_id, size_t target_id) { 
    old_message(source_id, target_id) = message(source_id, target_id);
  }
  
  /**
   * \brief Initialize the edge data with source and target having the
   * appropriate number of states.
   *
   * \param source_id the vertex id of the source
   * \param nsource the number of states the source vertex takes
   * \param target_id the vertex id of the target
   * \param ntarget the number of states the target vertex takes
   */
  void initialize(size_t source_id, size_t nsource, size_t target_id, size_t ntarget) {
    ASSERT_GT(nsource, 0); ASSERT_GT(ntarget, 0);
    message(source_id, target_id).setZero(ntarget);
    old_message(source_id, target_id).setZero(ntarget);
    message(target_id, source_id).setZero(nsource);
    old_message(target_id, source_id).setZero(nsource);
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



/**
 * \brief The graph type used to store the Markov Random Field with
 * vertex data containing node potentials and beliefs and edge data
 * containing messages and weights.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



/** 
 * \brief The Loopy Belief Propagation Vertex Program which computes
 * the product of the inbound messages during the gather phase,
 * updates the belief during the apply phase, and then computes the
 * new out-bound messages during the scatter phase.
 *
 * Since the gather phase is computing the product of the inbound
 * messages and the messages are stored in log form the resulting sum
 * operation is actually a vector sum and so the gather type is simply
 * the factor type and the operator+= operation for the factor type is
 * sufficient.
 *
 */
struct bp_vertex_program : 
  public graphlab::ivertex_program< graph_type, factor_type,
                                    graphlab::messages::sum_priority >,
  public graphlab::IS_POD_TYPE {

  /**
   * \brief Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /**
   * \brief Update the old message to be the new message and collect the
   * message value.
   */
  factor_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const {
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    edge_data& edata = edge.data();
    // Update the old message with the value of the new Message.  We
    // then receive the old message during gather and then compute the
    // "cavity" during scatter (again using the old message).
    edata.update_old(other_vertex.id(), vertex.id());
    const factor_type& recv_message = 
      edata.old_message(other_vertex.id(), vertex.id());
    // Ensure that the received message has the correct size
    ASSERT_EQ(recv_message.size(), vertex.data().potential.size());
    return recv_message;
  }; // end of gather function

  /**
   * \brief Multiply message product by node potential and update the
   * belief.
   */
  void apply(icontext_type& context, vertex_type& vertex, 
             const factor_type& total) {
    // If we have no neighbors than the belief is equal to the
    // potential so simply update the belief
    if(vertex.num_in_edges() + vertex.num_out_edges() == 0) {
      vertex.data().belief = vertex.data().potential;
    } else {
      vertex_data& vdata = vertex.data();
      ASSERT_EQ(vdata.potential.size(), total.size());
      // Multiply (add in log space) the potential to compute the belief
      vdata.belief = vdata.potential + total;
      ASSERT_GT(vdata.belief.size(), 0);
      // Rescale the belief to ensure numerical stability.  (This is
      // essentially normalization in log-space.)
      vdata.belief.array() -= vdata.belief.maxCoeff();
    }
  }; // end of apply

  /**
   * \brief Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  /**
   * \brief Compute new message value for each edge.
   */
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {  
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    edge_data& edata = edge.data();
    // Divide (subtract in log space) out of the belief the old in
    // message to construct the cavity
    const factor_type& old_in_message = 
      edata.old_message(other_vertex.id(), vertex.id());
    ASSERT_EQ(old_in_message.size(), vertex.data().belief.size());
    factor_type cavity = vertex.data().belief - old_in_message;
    // compute the new message by convolving with the Ising-Potts Edge
    // factor.
    factor_type& new_out_message = 
      edata.message(vertex.id(), other_vertex.id());
    const factor_type& old_out_message = 
      edata.old_message(vertex.id(), other_vertex.id());
    convolve(cavity, edata.weight(), new_out_message);
    // Renormalize (done in log space)
    new_out_message.array() -= new_out_message.maxCoeff();
    // Apply damping to the message to stabilize convergence.
    new_out_message = DAMPING * old_out_message + 
      (1-DAMPING) * new_out_message;
    // Compute message residual
    const double residual = 
      (new_out_message - old_out_message).cwiseAbs().sum();
    context.clear_gather_cache(other_vertex);
    // Schedule the adjacent vertex
    if(residual > TOLERANCE) context.signal(other_vertex, residual);
 }; // end of scatter

private:

  /**
   * \brief Compute the convolution of the cavity with the Ising-Potts
   * edge potential and store the result in the message
   *
   * \param cavity the belief minus the in-bound message
   * \param weight the edge weight used to scale the smoothing parameter
   * \param [out] message The message in which to store the result of
   * the convolution.
   */
  inline void convolve(const factor_type& cavity, const double& weight, 
                       factor_type& message) const {
    for(int i = 0; i < message.size(); ++i) {
      double sum = 0;
      for(int j = 0; j < cavity.size(); ++j) {
        sum += std::exp( cavity(j)  + ( i == j? 0 : -(SMOOTHING*weight) ) ); 
      }
      // To try and ensure numerical stability we do not allow
      // messages to underflow in log-space
      message(i) = (sum > 0)? std::log(sum) : 
        std::numeric_limits<double>::min();
    }
  } // end of convolve
  
  /**
   * \brief Given an edge and a vertex return the other vertex along
   * that edge. 
   */
  inline vertex_type get_other_vertex(edge_type& edge, 
                                      const vertex_type& vertex) const {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
  }; // end of other_vertex

}; // end of class bp_vertex_program







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
    clopts.print_description();
    return EXIT_FAILURE;
  }

  if(graph_dir.empty()) {
    logstream(LOG_ERROR) << "No graph was provided." << std::endl;
    clopts.print_description();
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





















