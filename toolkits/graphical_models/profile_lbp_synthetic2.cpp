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
 * \brief This application is almost identical to LBP structured
 * prediction except that it generates an artificial field for every
 * vertex enabling the studying of distributed LBP on various graph
 * structures without the need for actual data.  
 *
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


double FIELD = 2;
size_t NSTATES = 5;

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



bool USE_CACHE = false;


/**
 * Make a synthetic node potential
 */
factor_type make_node_potential(size_t vid) {
  // const size_t obs = vid % NSTATES;
  const size_t obs = 0;
  factor_type factor;
  factor.setZero(NSTATES);
  if(vid % 101 < 1) {
    for(size_t i = 0; i <  NSTATES; ++i) {
      factor(i) = (i == obs)? 0 : -FIELD;
    }
  } 
  return factor;
}

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
  void load(graphlab::iarchive& arc) { arc >> belief; }
  void save(graphlab::oarchive& arc) const { arc << belief; }
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
  factor_type messages_[2];
  /**
   * \brief The weight associated with the edge (used to scale the
   * smoothing parameter)
   */
  double weight_;
  /**
   * \brief The function used to compute the message index in the edge
   * message array.
   */
  size_t message_idx(size_t source_id, size_t target_id) {
    return size_t(source_id < target_id);
  }

public:

  edge_data(const double w = 1) : weight_(w) { }
  const double& weight() const { return weight_; }

  /**
   * \brief Get the new message value from source_id to target_id
   */
  factor_type& message(size_t source_id, size_t target_id) { 
    return messages_[message_idx(source_id, target_id)];
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
    message(target_id, source_id).setZero(nsource);
  }
  void save(graphlab::oarchive& arc) const {
    for(size_t i = 0; i < 2; ++i) arc << messages_[i];
    arc << weight_;
  }
  void load(graphlab::iarchive& arc) {
    for(size_t i = 0; i < 2; ++i) arc >> messages_[i];
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
  float residual;
  bp_vertex_program() : residual(0) { }

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
    const factor_type old_out_message = edata.message(vertex.id(), other_vertex.id());
    const factor_type old_in_message = edata.message(other_vertex.id(), vertex.id());
    factor_type cavity;
    if(other_vertex.data().belief.size()  == old_out_message.size()) {
      cavity = other_vertex.data().belief - old_out_message;
    } else { cavity = make_node_potential(other_vertex.id()); }
    factor_type new_in_message(old_in_message.size()); 
    convolve(cavity, edata.weight(), new_in_message);
    new_in_message.array() -= new_in_message.maxCoeff();
    new_in_message = DAMPING * old_in_message + (1-DAMPING) * new_in_message;
    new_in_message.array() -= new_in_message.maxCoeff();
    edata.message(other_vertex.id(), vertex.id()) = new_in_message;
    return new_in_message;
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
      vertex.data().belief = make_node_potential(vertex.id());
    } else {
      vertex_data& vdata = vertex.data();
      // Multiply (add in log space) the potential to compute the belief
      factor_type new_belief = make_node_potential(vertex.id()) + total;
      ASSERT_GT(new_belief.size(), 0);
      // Rescale the belief to ensure numerical stability.  (This is
      // essentially normalization in log-space.)
      new_belief.array() -= new_belief.maxCoeff();
      if(vdata.belief.size() != new_belief.size()) { residual = 1; }
      else { residual = (new_belief - vdata.belief).cwiseAbs().sum();}
      vdata.belief = new_belief;
    }
  }; // end of apply

  /**
   * \brief Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    if(USE_CACHE || residual > TOLERANCE)
      return graphlab::ALL_EDGES; 
    else return graphlab::NO_EDGES;
  }; // end of scatter edges

  /**
   * \brief Compute new message value for each edge.
   */
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {  
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    if(USE_CACHE) {
      context.clear_gather_cache(other_vertex);
    }
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
  if(source == target) return true;
  else {
    graph.add_edge(source, target, edge_data(weight));
    return true;
  }
} // end of edge loader



/**
 * \brief The edge initializer is used to allocate the messages along
 * each edge based on the number of states of the source and target
 * vertex.
 */
void edge_initializer(graph_type::edge_type& edge) {
  edge_data& edata = edge.data();
  const graphlab::vertex_id_type source_id = edge.source().id();
  const size_t nsource = NSTATES;
  const graphlab::vertex_id_type target_id = edge.target().id();
  const size_t ntarget = NSTATES;
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
 
  std::string graph_dir;
  std::string output_dir = "pred";
  std::string exec_type = "async";
  std::string format = "tsv";
  bool map = false;
  clopts.attach_option("graph", graph_dir,
                       "The directory containing the adjacency graph");
  clopts.add_positional("graph"); 
  clopts.attach_option("field", FIELD, 
                       "The background field used to construct the node potentials");
  clopts.attach_option("nstates", NSTATES, 
                       "The number of states for each variable");
  clopts.attach_option("cache", USE_CACHE, "use gather caching");
  clopts.attach_option("output", output_dir,
                       "The directory in which to save the predictions");
  clopts.attach_option("format", format, "The graph file format.");
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

  clopts.get_engine_args().set_option("use_cache", USE_CACHE);

  if(graph_dir.empty()) {
    logstream(LOG_ERROR) << "No graph was provided." << std::endl;
    return EXIT_FAILURE;
  }

  // Start the webserver
  graphlab::launch_metric_server();

  ///! load the graph
  graph_type graph(dc, clopts);  


  ///! load the graph
  graph.load_format(graph_dir, format);
  graph.finalize();
  dc.cout() << "Initializing edge data" << std::endl;
  graph.transform_edges(edge_initializer);

  typedef graphlab::omni_engine<bp_vertex_program> engine_type;
  engine_type engine(dc, graph, exec_type, clopts);
  engine.signal_all();
  graphlab::timer timer;
  dc.cout() << "Running engine" << std::endl;
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





















