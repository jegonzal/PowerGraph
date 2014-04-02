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
 * \file
 * This code iplements the NMF algorithm described in the paper:
 * Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix
 * Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>

#include <Eigen/Dense>
#include "eigen_serialization.hpp"
#include <graphlab/macros_def.hpp>
#include <graphlab/util/timer.hpp>
#include "stats.hpp"

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat_type;

//when using negative node id range, we are not allowed to use
//0 and 1 so we add 2.
const static int SAFE_NEG_OFFSET=2;
const double epsilon = 1e-16;
static bool debug;
int iter = 0;

bool isuser(uint node){
  return ((int)node) >= 0;
}
/** 
 * \ingroup toolkit_matrix_pvecization
 *
 * \brief the vertex data type which contains the latent pvec.
 *
 * Each row and each column in the matrix corresponds to a different
 * vertex in the SGD graph.  Associated with each vertex is a pvec
 * (vector) of latent parameters that represent that vertex.  The goal
 * of the SGD algorithm is to find the values for these latent
 * parameters such that the non-zero entries in the matrix can be
 * predicted by taking the dot product of the row and column pvecs.
 */
struct vertex_data {
  /**
   * \brief A shared "constant" that specifies the number of latent
   * values to use.
   */
  static size_t NLATENT;
  /** \brief The latent pvec for this vertex */
  vec pvec;

  double train_rmse;
  double validation_rmse;
  /** 
   * \brief Simple default constructor which randomizes the vertex
   *  data 
   */
  vertex_data() { if (debug) pvec = vec::Ones(NLATENT); else randomize(); train_rmse = validation_rmse = 0; } 
  /** \brief Randomizes the latent pvec */
  void randomize() { pvec.resize(NLATENT); pvec.setRandom(); }
  /** \brief Save the vertex data to a binary archive */
  void save(graphlab::oarchive& arc) const { 
    arc << pvec << train_rmse << validation_rmse;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) { 
    arc >> pvec >>train_rmse >> validation_rmse;
  }
}; // end of vertex data


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data nmfo stores the most recent error estimate.
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  /**
   * \brief The type of data on the edge;
   *
   * \li *Train:* the observed value is correct and used in training
   * \li *Validate:* the observed value is correct but not used in training
   * \li *Predict:* The observed value is not correct and should not be
   *        used in training.
   */
  enum data_role_type { TRAIN, VALIDATE, PREDICT  };

  /** \brief the observed value for the edge */
  float weight;

  /** \brief The train/validation/test designation of the edge */
  data_role_type role;

  /** \brief basic initialization */
  edge_data(float weight = 0, data_role_type role = PREDICT) :
    weight(weight), role(role) { }

}; // end of edge data



/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



double extract_l2_error(const graph_type::edge_type & edge);


/**
 * \brief Given a vertex and an edge return the other vertex in the
 * edge.
 */
inline graph_type::vertex_type
get_other_vertex(graph_type::edge_type& edge, 
    const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex


class gather_type {
  public:
    vec pvec;
    double training_rmse;
    double validation_rmse;
    gather_type() { training_rmse = validation_rmse = 0; }
    gather_type(const vec & _pvec, double _train_rmse, double _validation_rmse){ pvec = _pvec; training_rmse = _train_rmse; validation_rmse = _validation_rmse; }
    void reset(){ pvec = vec::Zero(vertex_data::NLATENT); training_rmse = 0; validation_rmse = 0; }
    void save(graphlab::oarchive& arc) const { arc << pvec << training_rmse << validation_rmse; }
    void load(graphlab::iarchive& arc) { arc >> pvec >> training_rmse >> validation_rmse; }  
    gather_type& operator+=(const gather_type& other) {
      pvec += other.pvec;
      training_rmse += other.training_rmse;
      validation_rmse += other.validation_rmse;
      return *this;
    } 

};

gather_type x1;
gather_type x2;
gather_type * px;


bool isuser_node(const graph_type::vertex_type& vertex){
  return isuser(vertex.id());
}
/**
 * SGD vertex program type
 */ 
class nmf_vertex_program :
  public graphlab::ivertex_program<graph_type, gather_type, gather_type>,
  public graphlab::IS_POD_TYPE{
    public:
      /** The convergence tolerance */
      static double TOLERANCE;
      static double MAXVAL;
      static double MINVAL;
      static bool debug;
      static size_t MAX_UPDATES;

      /** compute a missing value based on NMF algorithm */
      static float nmf_predict(const vertex_data& user, 
          const vertex_data& movie, 
          const float rating, 
          double & prediction){

        prediction = user.pvec.dot(movie.pvec);
        //truncate prediction to allowed values
        prediction = std::min((double)prediction, nmf_vertex_program::MAXVAL);
        prediction = std::max((double)prediction, nmf_vertex_program::MINVAL);
        //return the squared error
        float err = rating - prediction;
        assert(!std::isnan(err));
        return err*err; 

      }


      /** The set of edges to gather along */
      edge_dir_type gather_edges(icontext_type& context, 
          const vertex_type& vertex) const { 
        //UNUSED 
        return graphlab::ALL_EDGES; 
      }; // end of gather_edges 

      /** The gather function computes XtX and Xy */
      gather_type gather(icontext_type& context, const vertex_type& vertex, 
          edge_type& edge) const {

        if (edge.data().role == edge_data::TRAIN || edge.data().role == edge_data::VALIDATE){
          const vertex_type other_vertex = get_other_vertex(edge, vertex);
          double prediction = 0;
          double rmse = nmf_predict(vertex.data(), other_vertex.data(), edge.data().weight, prediction);
          if (prediction == 0)
            logstream(LOG_FATAL)<<"Got into numerical error!" << std::endl;
          if (edge.data().role == edge_data::TRAIN)
            return gather_type(other_vertex.data().pvec * (edge.data().weight / prediction), rmse, 0);
          else //validation
            return gather_type(vec::Zero(vertex_data::NLATENT), 0, rmse);

        }
        return gather_type(vec::Zero(vertex_data::NLATENT), 0, 0);
      } // end of gather function

      void apply(icontext_type& context, vertex_type& vertex,
          const gather_type& sum) {
        vertex_data& vdata = vertex.data();  
        if (vdata.pvec.sum() != 0){
          for (uint i=0; i< vertex_data::NLATENT; i++){
            vdata.pvec[i] *= sum.pvec[i] / px->pvec[i];
            ASSERT_NE(px->pvec[i] , 0);
            if (vdata.pvec[i] < epsilon)
              vdata.pvec[i] = epsilon;
          }
        }
        vdata.train_rmse = sum.training_rmse;
        vdata.validation_rmse = sum.validation_rmse;
      }

      edge_dir_type scatter_edges(icontext_type& context,
          const vertex_type& vertex) const { 
        //UNUSED 
        return graphlab::ALL_EDGES; 
      }; // end of scatter edges

      void scatter(icontext_type& context, const vertex_type& vertex, 
          edge_type& edge) const {
        //we do not schedule any more neighbors to run
      } 
      static void verify_rows(graph_type::vertex_type& vertex){
        if (isuser(vertex.id()) && vertex.num_out_edges() == 0)
          logstream(LOG_FATAL)<<"NMF algorithm can not work when the row " << vertex.id() << " of the matrix contains all zeros" << std::endl;
      }

      static gather_type pre_iter(const graph_type::vertex_type & vertex){
        gather_type ret;
        ret.pvec = vertex.data().pvec;
        ret.training_rmse = vertex.data().train_rmse;
        ret.validation_rmse = vertex.data().validation_rmse;
        return ret;
      }


      static graphlab::empty signal_left(icontext_type& context,
          const vertex_type& vertex) {
        if(vertex.num_out_edges() > 0) context.signal(vertex);
        return graphlab::empty();
      } // end of signal_left 

      static graphlab::empty signal_right(icontext_type& context,
          const vertex_type& vertex) {
        if(vertex.num_in_edges() > 0) context.signal(vertex);
        return graphlab::empty();
      } // end of signal_left 

  }; // end of nmf vertex program

gather_type count_edges(nmf_vertex_program::icontext_type & context, const graph_type::edge_type& edge) {
  gather_type ret;
  if (edge.data().role == edge_data::TRAIN){
    ret.training_rmse = 1;
  }
  else if (edge.data().role == edge_data::VALIDATE){
    ret.validation_rmse = 1;
  }
  if (edge.data().weight < 0)
    logstream(LOG_FATAL)<<"Found a negative entry in matirx row " << edge.source().id() << " with value: " << edge.data().weight << std::endl;
  return ret;
}


struct prediction_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
     */
  std::string save_vertex(const vertex_type& vertex) const {
    return "";
  }
  std::string save_edge(const edge_type& edge) const {
    if (edge.data().role != edge_data::PREDICT)
      return "";

    std::stringstream strm;
    const double prediction = 
      edge.source().data().pvec.dot(edge.target().data().pvec);
    strm << edge.source().id() << '\t' 
      << -edge.target().id()-SAFE_NEG_OFFSET << '\t'
      << prediction << '\n';
    return strm.str();
  }
}; // end of prediction_saver

struct linear_model_saver_U {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
     */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() > 0){
      std::string ret = boost::lexical_cast<std::string>(vertex.id()) + ") ";
      for (uint i=0; i< vertex_data::NLATENT; i++)
        ret += boost::lexical_cast<std::string>(vertex.data().pvec[i]) + " ";
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
      std::string ret = boost::lexical_cast<std::string>(-vertex.id()-SAFE_NEG_OFFSET) + ") ";
      for (uint i=0; i< vertex_data::NLATENT; i++)
        ret += boost::lexical_cast<std::string>(vertex.data().pvec[i]) + " ";
      ret += "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
}; 



/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph, 
    const std::string& filename,
    const std::string& line) {
  ASSERT_FALSE(line.empty()); 

 // Parse the line
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float weight(0);
  strm >> source_id >> target_id;

  if (source_id == graph_type::vertex_id_type(-1) || target_id == graph_type::vertex_id_type(-1)){
    logstream(LOG_WARNING)<<"Failed to read input line: "<< line << " in file: "  << filename << " (or node id is -1). " << std::endl;
    return true;
  }

  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
  else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;
 
  // for test files (.predict) no need to read the actual rating value.
  if(role == edge_data::TRAIN || role == edge_data::VALIDATE){
    strm >> weight;
    if (weight < nmf_vertex_program::MINVAL || weight > nmf_vertex_program::MAXVAL)
      logstream(LOG_FATAL)<<"Rating values should be between " << nmf_vertex_program::MINVAL << " and " << nmf_vertex_program::MAXVAL << ". Got value: " << weight << " [ user: " << source_id << " to item: " <<target_id << " ] " << std::endl; 
  }
  target_id = -(graphlab::vertex_id_type(target_id + SAFE_NEG_OFFSET));

  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(weight, role)); 
  return true; // successful load
} // end of graph_loader






size_t vertex_data::NLATENT = 20;
double nmf_vertex_program::TOLERANCE = 1e-3;
size_t nmf_vertex_program::MAX_UPDATES = -1;
double nmf_vertex_program::MAXVAL = 1e+100;
double nmf_vertex_program::MINVAL = -1e+100;
bool nmf_vertex_program::debug = false;


/**
 * \brief The engine type used by the ALS matrix factorization
 * algorithm.
 *
 * The ALS matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef graphlab::omni_engine<nmf_vertex_program> engine_type;

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir;
  std::string predictions;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
      "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D", vertex_data::NLATENT,
      "Number of latent parameters to use.");
  clopts.attach_option("engine", exec_type, 
      "The engine type synchronous or asynchronous");
  clopts.attach_option("max_iter", nmf_vertex_program::MAX_UPDATES,
      "The maxumum number of udpates allowed for a vertex");
  clopts.attach_option("debug", nmf_vertex_program::debug, 
      "debug - additional verbose info"); 
  clopts.attach_option("maxval", nmf_vertex_program::MAXVAL, "max allowed value");
  clopts.attach_option("minval", nmf_vertex_program::MINVAL, "min allowed value");
  clopts.attach_option("predictions", predictions,
      "The prefix (folder and filename) to save predictions.");

  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  debug = nmf_vertex_program::debug;
  
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; 
  graph_type graph(dc, clopts);  
  graph.load(input_dir, graph_loader); 
  dc.cout() << "Loading graph. Finished in " 
    << timer.current_time() << std::endl;
  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " 
    << timer.current_time() << std::endl;

  if (!graph.num_edges() || !graph.num_vertices())
     logstream(LOG_FATAL)<< "Failed to load graph. Check your input path: " << input_dir << std::endl;     



  dc.cout() 
    << "========== Graph statistics on proc " << dc.procid() 
    << " ==============="
    << "\n Num vertices: " << graph.num_vertices()
    << "\n Num edges: " << graph.num_edges()
    << "\n Num replica: " << graph.num_replicas()
    << "\n Replica to vertex ratio: " 
    << float(graph.num_replicas())/graph.num_vertices()
    << "\n --------------------------------------------" 
    << "\n Num local own vertices: " << graph.num_local_own_vertices()
    << "\n Num local vertices: " << graph.num_local_vertices()
    << "\n Replica to own ratio: " 
    << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
    << "\n Num local edges: " << graph.num_local_edges()
    //<< "\n Begin edge id: " << graph.global_eid(0)
    << "\n Edge balance ratio: " 
    << float(graph.num_local_edges())/graph.num_edges()
    << std::endl;

  dc.cout() << "Creating engine" << std::endl;
  engine_type engine(dc, graph, exec_type, clopts);


  // Run the NMF ---------------------------------------------------------
  dc.cout() << "Running NMF" << std::endl;
  dc.cout() << "(C) Code by Danny Bickson, CMU " << std::endl;
  dc.cout() << "Please send bug reports to danny.bickson@gmail.com" << std::endl;
  dc.cout() << "Time   Training    Validation" <<std::endl;
  dc.cout() << "       RMSE        RMSE " <<std::endl;
  timer.start();

  gather_type edge_count = engine.map_reduce_edges<gather_type>(count_edges);
  dc.cout()<<"Training edges: " << edge_count.training_rmse << " validation edges: " << edge_count.validation_rmse << std::endl;

  graphlab::vertex_set left = graph.select(isuser_node);
  graphlab::vertex_set right = ~left;
  graph.transform_vertices(nmf_vertex_program::verify_rows, left);

  graphlab::timer mytimer; mytimer.start();

  for (uint j=0; j< nmf_vertex_program::MAX_UPDATES; j++){
    x1 = graph.map_reduce_vertices<gather_type>(nmf_vertex_program::pre_iter,right);
    px = &x1;
    for (int i=0; i< (int)vertex_data::NLATENT; i++)
      ASSERT_NE(px->pvec[i], 0);
    
    dc.cout()<< std::setw(8) << mytimer.current_time() << " " << sqrt(x1.training_rmse/edge_count.training_rmse);
    if (edge_count.validation_rmse > 0)
      dc.cout() << " " << std::setw(8) << sqrt(x1.validation_rmse/edge_count.validation_rmse) << std::endl;
    else dc.cout() << std::endl;
    engine.map_reduce_vertices<graphlab::empty>(nmf_vertex_program::signal_left);
    engine.start();
    x1.reset();

    x2 = graph.map_reduce_vertices<gather_type>(nmf_vertex_program::pre_iter,left);
    px = &x2;

    engine.map_reduce_vertices<graphlab::empty>(nmf_vertex_program::signal_right);
    engine.start();
    x2.reset();
  }

  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
    << std::endl
    << "Final Runtime (seconds):   " << runtime 
                                        << std::endl
                                        << "Updates executed: " << engine.num_updates() << std::endl
                                        << "Update Rate (updates/second): " 
                                          << engine.num_updates() / runtime << std::endl;


  // Make predictions ---------------------------------------------------------
  if(!predictions.empty()) {
    std::cout << "Saving predictions" << std::endl;
    const bool gzip_output = false;
    const bool save_vertices = false;
    const bool save_edges = true;
    const size_t threads_per_machine = 1;
    //save the predictions
    graph.save(predictions, prediction_saver(),
        gzip_output, save_vertices, 
        save_edges, threads_per_machine);
    //save the linear model
    graph.save(predictions + ".U", linear_model_saver_U(),
        gzip_output, save_edges, save_vertices, threads_per_machine);
    graph.save(predictions + ".V", linear_model_saver_V(),
        gzip_output, save_edges, save_vertices, threads_per_machine);

  }

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



