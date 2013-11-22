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
 *      http://graphlab.org
 *
 */


/**
 * \file
 * 
 * \brief The main file for the BIAS-SGD matrix factorization algorithm.
 *
 * This file contains the main body of the BIAS-SGD matrix factorization
 * algorithm. 
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>
#include "eigen_serialization.hpp"
#include <Eigen/Dense>
#include <graphlab/macros_def.hpp>



typedef Eigen::VectorXd vec_type;
typedef Eigen::MatrixXd mat_type;

//when using negative node id range, we are not allowed to use
//0 and 1 so we add 2.
const static int SAFE_NEG_OFFSET=2;
static bool debug;
int iter = 0;
float itmBiasStep = 1e-4;
float itmBiasReg = 1e-4;
float usrBiasStep = 1e-4;
float usrBiasReg = 1e-4;
float usrFctrStep = 1e-4;
float usrFctrReg = 1e-4;
float itmFctrStep = 1e-4;
float itmFctrReg = 1e-4; //gamma7
float itmFctr2Step = 1e-4;
float itmFctr2Reg = 1e-4;
/** 
 * \ingroup toolkit_matrix_pvecization
 *
 * \brief the vertex data type which contains the latent pvec.
 *
 * Each row and each column in the matrix corresponds to a different
 * vertex in the BIASSGD graph.  Associated with each vertex is a pvec
 * (vector) of latent parameters that represent that vertex.  The goal
 * of the BIASSGD algorithm is to find the values for these latent
 * parameters such that the non-zero entries in the matrix can be
 * predicted by taking the dot product of the row and column pvecs.
 */
struct vertex_data {
  /**
   * \brief A shared "constant" that specifies the number of latent
   * values to use.
   */
  static size_t NLATENT;
  /** \brief The number of times this vertex has been updated. */
  uint32_t nupdates;
  /** \brief The latent pvec for this vertex */
  vec_type pvec;
  vec_type weight;
  double bias;
  /** 
   * \brief Simple default constructor which randomizes the vertex
   *  data 
   */
  vertex_data() : nupdates(0) { randomize(); } 
  /** \brief Randomizes the latent pvec */
  void randomize() { pvec.resize(NLATENT); pvec.setRandom(); weight.resize(NLATENT); weight.setRandom(); }
  /** \brief Save the vertex data to a binary archive */
  void save(graphlab::oarchive& arc) const { 
    arc << nupdates << pvec << weight << bias;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) { 
    arc >> nupdates >> pvec >> weight >> bias;
  }
}; // end of vertex data


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data svdppo stores the most recent error estimate.
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
  float obs;

  /** \brief The train/validation/test designation of the edge */
  data_role_type role;

  /** \brief basic initialization */
  edge_data(float obs = 0, data_role_type role = PREDICT) :
    obs(obs), role(role) { }

}; // end of edge data


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

#include "implicit.hpp"

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





/**
 * \brief The gather type used to construct XtX and Xty needed for the BIASSGD
 * update
 *
 * To compute the ALS update we need to compute the sum of 
 * \code
 *  sum: XtX = nbr.pvec.transpose() * nbr.pvec 
 *  sum: Xy  = nbr.pvec * edge.obs
 * \endcode
 * For each of the neighbors of a vertex. 
 *
 * To do this in the Gather-Apply-Scatter model the gather function
 * computes and returns a pair consisting of XtX and Xy which are then
 * added. The gather type represents that tuple and provides the
 * necessary gather_type::operator+= operation.
 *
 */
class gather_type {
  public:
    /**
     * \brief Stores the current sum of nbr.pvec.transpose() *
     * nbr.pvec
     */

    /**
     * \brief Stores the current sum of nbr.pvec * edge.obs
     */
    vec_type pvec;
    vec_type weight;
    double bias;

    /** \brief basic default constructor */
    gather_type() { }

    /**
     * \brief This constructor computes XtX and Xy and stores the result
     * in XtX and Xy
     */
    gather_type(const vec_type& X, const vec_type & _weight, double _bias) {
      pvec = X;
      bias = _bias;
      weight = _weight;
    } // end of constructor for gather type

    /** \brief Save the values to a binary archive */
    void save(graphlab::oarchive& arc) const { arc << pvec << bias << weight; }

    /** \brief Read the values from a binary archive */
    void load(graphlab::iarchive& arc) { arc >> pvec >> bias >> weight; }  

    /** 
     * \brief Computes XtX += other.XtX and Xy += other.Xy updating this
     * tuples value
     */
    gather_type& operator+=(const gather_type& other) {
      if (pvec.size() == 0){
        pvec = other.pvec;
        bias = other.bias;
        weight = other.weight;
        return *this;
      }
      else if (other.pvec.size() == 0)
        return *this;
      pvec += other.pvec;
      bias += other.bias;
      weight += other.weight;
      return *this;
    } // end of operator+=

}; // end of gather type

//typedef gather_type message_type;


enum{
  PHASE1 = 0, PHASE2 = 1
};

/**
 * BIASSGD vertex program type
 */ 
class svdpp_vertex_program : 
  public graphlab::ivertex_program<graph_type, gather_type,
  gather_type> {
    public:
      /** The convergence tolerance */
      static double TOLERANCE;
      static double LAMBDA;
      static double GAMMA;
      static double MAXVAL;
      static double MINVAL;
      static double STEP_DEC;
      static bool debug;
      static size_t MAX_UPDATES;
      static double GLOBAL_MEAN;
      static size_t NUM_TRAINING_EDGES;
      static uint   USERS;

      gather_type pmsg;
      void save(graphlab::oarchive& arc) const { 
        arc << pmsg;
      }
      /** \brief Load the vertex data from a binary archive */
      void load(graphlab::iarchive& arc) { 
        arc >> pmsg;
      }

      /** The set of edges to gather along */
      edge_dir_type gather_edges(icontext_type& context, 
          const vertex_type& vertex) const { 
        return graphlab::ALL_EDGES; 
      }; // end of gather_edges 

      gather_type gather(icontext_type& context, const vertex_type& vertex, 
          edge_type& edge) const {
        vec_type step = vec_type::Zero(vertex_data::NLATENT);
        double bias =0, other_bias = 0;
        vec_type delta, other_delta;

        //user node
        if (vertex.num_in_edges() == 0){
          vertex_type other_vertex(get_other_vertex(edge, vertex));
          vertex_type my_vertex(vertex);

          int phase = my_vertex.data().nupdates % 2; 
          if (phase == PHASE1){
            //my_vertex.data().weight += movie.weight;
            context.signal(other_vertex, gather_type(vec_type::Zero(vertex_data::NLATENT), vec_type::Zero(vertex_data::NLATENT), 0));
            return gather_type(vec_type::Zero(vertex_data::NLATENT), other_vertex.data().weight, 0);
          }
          else if (phase == PHASE2){
            //vertex_data & my_data = my_vertex.data();
            double pred = svdpp_vertex_program::GLOBAL_MEAN + 
              my_vertex.data().bias + other_vertex.data().bias + my_vertex.data().pvec.dot(other_vertex.data().pvec+other_vertex.data().weight);
            pred = std::min(pred, svdpp_vertex_program::MAXVAL);
            pred = std::max(pred, svdpp_vertex_program::MINVAL); 
            const float err = edge.data().obs - pred;
            if (debug)
              std::cout<<"entering edge " << (int)edge.source().id() << ":" << (int)edge.target().id() << " err: " << err << " rmse: " << err*err <<std::endl;
            if (std::isnan(err))
              logstream(LOG_FATAL)<<"Got into numeric errors.. try to tune step size and regularization using command line flags" << std::endl;
            if (edge.data().role == edge_data::TRAIN){
              vec_type itmFctr = other_vertex.data().pvec;
              vec_type usrFctr = my_vertex.data().pvec;

              bias = usrBiasStep*(err - usrBiasReg*bias);
              other_bias = itmBiasStep*(err - itmBiasReg*other_bias);

              delta = usrFctrStep*(err*(itmFctr - usrFctrReg *usrFctr));
              other_delta = itmFctrStep*(err*(usrFctr+my_vertex.data().weight) - itmFctrReg*other_vertex.data().pvec);

              step = err*itmFctr;
              float usrNorm = double(1.0/sqrt(my_vertex.num_out_edges()));
              step *= itmFctr2Step*usrNorm;

              double mult = itmFctr2Step*itmFctr2Reg;
              step -= mult*other_vertex.data().weight; 
              //A HACK: update memory cached values to reflect new vals 
              /*my_vertex.data().bias += bias;
                other_vertex.data().bias += other_bias;
                my_vertex.data().pvec += delta;
                other_vertex.data().pvec += other_delta;*/

              if (debug)
                std::cout<<"new val:" << (int)edge.source().id() << ":" << (int)edge.target().id() << " U " << my_vertex.data().pvec.transpose() << " V " << other_vertex.data().pvec.transpose() << std::endl;

              if(other_vertex.data().nupdates < MAX_UPDATES) 
                context.signal(other_vertex, gather_type(other_delta, step, other_bias));
            }
            return gather_type(delta, step, bias);

          } //end of PHASE2
        }
        return gather_type(delta, step, bias);
      }  

      //typedef vec_type message_type;
      void init(icontext_type& context,
          const vertex_type& vertex,
          const message_type& msg) {

        int phase = vertex.data().nupdates % 2;
        //movie node receives updates here
        if (vertex.num_in_edges() > 0){
          if (phase == PHASE1){
            pmsg = msg;
          }
          else if (phase == PHASE2){
            pmsg = msg;
          }
        }
      }

      /** apply collects the sum of XtX and Xy */
      void apply(icontext_type& context, vertex_type& vertex,
          const gather_type& sum) {
        vertex_data& vdata = vertex.data(); 
        int phase = vdata.nupdates %2;

        if (phase == PHASE1){
          //user node receives the sum of movie weights
          if (vertex.num_out_edges() > 0){
            vertex.data().weight = sum.pvec;
            float usrNorm = double(1.0/sqrt(vertex.num_out_edges()));
            vertex.data().weight *= usrNorm;
          }
          //movie node doe nothing
          else {}
        }
        else if (phase == PHASE2){
          //user node update gradients and bias
          if (vertex.num_in_edges() == 0){
            vdata.pvec += sum.pvec;
            vdata.bias += sum.bias;
            //does not update weight here (since was done in phase1)
          }
          //movie node
          else {
            vdata.weight += pmsg.weight; //step
            vdata.pvec += pmsg.pvec;
            vdata.bias += pmsg.bias;
          }
        }
        ++vdata.nupdates;
      } // end of apply

      /** The edges to scatter along */
      edge_dir_type scatter_edges(icontext_type& context,
          const vertex_type& vertex) const { 
        return graphlab::ALL_EDGES; 
      }; // end of scatter edges

      /** Scatter reschedules neighbors */  
      void scatter(icontext_type& context, const vertex_type& vertex, 
          edge_type& edge) const {
        edge_data& edata = edge.data();
        if(edata.role == edge_data::TRAIN) {
          const vertex_type other_vertex = get_other_vertex(edge, vertex);
          // Reschedule neighbors ------------------------------------------------
          if(other_vertex.data().nupdates < MAX_UPDATES) 
            context.signal(other_vertex, gather_type(vec_type::Zero(vertex_data::NLATENT),vec_type::Zero(vertex_data::NLATENT),0));
        }
      } // end of scatter function


      /**
       * \brief Signal all vertices on one side of the bipartite graph
       */
      static graphlab::empty signal_left(icontext_type& context,
          vertex_type& vertex) {
        if(vertex.num_out_edges() > 0) context.signal(vertex, gather_type(vec_type::Zero(vertex_data::NLATENT),vec_type::Zero(vertex_data::NLATENT),0));
        return graphlab::empty();
      } // end of signal_left 

  }; // end of svdpp vertex program


struct error_aggregator : public graphlab::IS_POD_TYPE {
  typedef svdpp_vertex_program::icontext_type icontext_type;
  typedef graph_type::edge_type edge_type;
  double train_error, validation_error;
  size_t ntrain, nvalidation;
  error_aggregator() : 
    train_error(0), validation_error(0), ntrain(0), nvalidation(0) { }
  error_aggregator& operator+=(const error_aggregator& other) {
    train_error += other.train_error;
    assert(!std::isnan(train_error));
    validation_error += other.validation_error;
    ntrain += other.ntrain;
    nvalidation += other.nvalidation;
    return *this;
  }
  static error_aggregator map(icontext_type& context, const graph_type::edge_type& edge) {
    error_aggregator agg;
    if (edge.data().role == edge_data::TRAIN){
      agg.train_error = extract_l2_error(edge); agg.ntrain = 1;
      assert(!std::isnan(agg.train_error));
    }
    else if (edge.data().role == edge_data::VALIDATE){
      agg.validation_error = extract_l2_error(edge); agg.nvalidation = 1;
    }
    return agg;
  }


  static void finalize(icontext_type& context, const error_aggregator& agg) {
    iter++;
    if (iter%2 == 0)
      return; 
    ASSERT_GT(agg.ntrain, 0);
    const double train_error = std::sqrt(agg.train_error / agg.ntrain);
    assert(!std::isnan(train_error));
    context.cout() << std::setw(8) << context.elapsed_seconds() << std::setw(8) << train_error;
    if(agg.nvalidation > 0) {
      const double validation_error = 
        std::sqrt(agg.validation_error / agg.nvalidation);
      context.cout() << std::setw(8) << validation_error; 
    }
    context.cout() << std::endl;
    usrBiasStep *= svdpp_vertex_program::STEP_DEC;
    itmBiasStep *= svdpp_vertex_program::STEP_DEC;
    usrFctrStep  *= svdpp_vertex_program::STEP_DEC;
    itmFctrStep  *= svdpp_vertex_program::STEP_DEC;
    itmFctr2Step *= svdpp_vertex_program::STEP_DEC;

  }
}; // end of error aggregator

/**
 * \brief Given an edge compute the error associated with that edge
 */
double extract_l2_error(const graph_type::edge_type & edge) {
  double pred = svdpp_vertex_program::GLOBAL_MEAN + 
    edge.source().data().bias +
    edge.target().data().bias + 
    edge.source().data().pvec.dot(edge.target().data().pvec);
  pred = std::min(svdpp_vertex_program::MAXVAL, pred);
  pred = std::max(svdpp_vertex_program::MINVAL, pred);
  double rmse = (edge.data().obs - pred) * (edge.data().obs - pred);
  assert(rmse <= pow(svdpp_vertex_program::MAXVAL-svdpp_vertex_program::MINVAL,2));
  return rmse;
} // end of extract_l2_error


struct prediction_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    return ""; //nop
  }
  std::string save_edge(const edge_type& edge) const {
    if (edge.data().role != edge_data::PREDICT)
      return "";

    std::stringstream strm;
    double pred = svdpp_vertex_program::GLOBAL_MEAN +
      edge.target().data().bias + edge.source().data().bias + edge.source().data().pvec.dot(edge.target().data().pvec+edge.target().data().weight);
      pred = std::min(pred, svdpp_vertex_program::MAXVAL);
      pred = std::max(pred, svdpp_vertex_program::MINVAL);
    strm << edge.source().id() << '\t' 
      << -edge.target().id()-SAFE_NEG_OFFSET << '\t'
      << pred << '\n';
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
      std::string ret = boost::lexical_cast<std::string>(vertex.id()) + " ";
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
      std::string ret = boost::lexical_cast<std::string>(-vertex.id()-SAFE_NEG_OFFSET) + " ";
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

struct linear_model_saver_bias_U {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
     */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() > 0){
      std::string ret = boost::lexical_cast<std::string>(vertex.id()) + " ";
      ret += boost::lexical_cast<std::string>(vertex.data().bias) + "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
}; 
struct linear_model_saver_bias_V {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
     */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() == 0){
      std::string ret = boost::lexical_cast<std::string>(-vertex.id()-SAFE_NEG_OFFSET) + " ";
      ret += boost::lexical_cast<std::string>(vertex.data().bias) + "\n";
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

 // Parse the line
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  strm >> source_id >> target_id;

  if (source_id == graph_type::vertex_id_type(-1) || target_id == graph_type::vertex_id_type(-1)){
    logstream(LOG_WARNING)<<"Failed to read input line: "<< line << " in file: "  << filename << " (or node id is -1). " << std::endl;
    return true;
  }

  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
  else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;
 
  if(role == edge_data::TRAIN || role == edge_data::VALIDATE){
    strm >> obs;
    if (obs < svdpp_vertex_program::MINVAL || obs > svdpp_vertex_program::MAXVAL){
      logstream(LOG_WARNING)<<"Rating values should be between " << svdpp_vertex_program::MINVAL << " and " << svdpp_vertex_program::MAXVAL << ". Got value: " << obs << " [ user: " << source_id << " to item: " <<target_id << " ] " << std::endl; 
      assert(false); 
    }
  }
  target_id = -(graphlab::vertex_id_type(target_id + SAFE_NEG_OFFSET));
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(obs, role)); 
  return true; // successful load
} // end of graph_loader






size_t vertex_data::NLATENT = 20;
double svdpp_vertex_program::TOLERANCE = 1e-3;
double svdpp_vertex_program::LAMBDA = 0.001;
double svdpp_vertex_program::GAMMA = 0.001;
size_t svdpp_vertex_program::MAX_UPDATES = -1;
double svdpp_vertex_program::MAXVAL = 1e+100;
double svdpp_vertex_program::MINVAL = -1e+100;
double svdpp_vertex_program::STEP_DEC = 0.9;
bool svdpp_vertex_program::debug = false;
double svdpp_vertex_program::GLOBAL_MEAN = 0;
size_t svdpp_vertex_program::NUM_TRAINING_EDGES = 0;

/**
 * \brief The engine type used by the ALS matrix factorization
 * algorithm.
 *
 * The ALS matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef graphlab::omni_engine<svdpp_vertex_program> engine_type;

  double calc_global_mean(const graph_type::edge_type & edge){
    if (edge.data().role == edge_data::TRAIN)
      return edge.data().obs;
    else return 0;
  }

  size_t count_edges(const graph_type::edge_type & edge){
    if (edge.data().role == edge_data::TRAIN)
      return 1;
    else return 0;
  }


int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir, output_dir;
  std::string predictions;
  size_t interval = 0;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
      "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D", vertex_data::NLATENT,
      "Number of latent parameters to use.");
  clopts.attach_option("engine", exec_type, 
      "The engine type synchronous or asynchronous");
  clopts.attach_option("max_iter", svdpp_vertex_program::MAX_UPDATES,
      "The maxumum number of udpates allowed for a vertex");
  clopts.attach_option("lambda", svdpp_vertex_program::LAMBDA, 
      "SGD regularization weight"); 
  clopts.attach_option("gamma", svdpp_vertex_program::GAMMA, 
      "SGD step size"); 
  clopts.attach_option("debug", svdpp_vertex_program::debug, 
      "debug - additional verbose info"); 
  clopts.attach_option("tol", svdpp_vertex_program::TOLERANCE,
      "residual termination threshold");
  clopts.attach_option("maxval", svdpp_vertex_program::MAXVAL, "max allowed value");
  clopts.attach_option("minval", svdpp_vertex_program::MINVAL, "min allowed value");
  clopts.attach_option("step_dec", svdpp_vertex_program::STEP_DEC, "multiplicative step decrement");
  clopts.attach_option("user_bias_step", usrBiasStep, "user_bias_step");
  clopts.attach_option("user_bias_reg", usrBiasReg, "user_bias_reg");
  clopts.attach_option("item_bias_step",itmBiasStep, "item_bias_step");
  clopts.attach_option("item_bias_reg", itmBiasReg, "item_bias_reg");
  clopts.attach_option("user_factor_step", usrFctrStep, "user_factor_step");
  clopts.attach_option("user_factor_reg", usrFctrReg, "user_factor_reg");
  clopts.attach_option("item_factor_step", itmFctrStep, "item_factor_step");
  clopts.attach_option("item_factor_reg", itmFctrReg, "item_factor_reg");
  clopts.attach_option("item_factor2_step", itmFctr2Step, "item_factor2_step");
  clopts.attach_option("item_factor2_reg", itmFctr2Reg, "item_factor2_reg");
  clopts.attach_option("interval", interval, "The time in seconds between error reports");
  clopts.attach_option("predictions", predictions,
      "The prefix (folder and filename) to save predictions.");
  clopts.attach_option("output", output_dir, "Output results");

  parse_implicit_command_line(clopts);

  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  debug = svdpp_vertex_program::debug;
  //  omp_set_num_threads(clopts.get_ncpus());
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; 
  graph_type graph(dc, clopts);  
  graph.load(input_dir, graph_loader); 
  dc.cout() << "Loading graph. Finished in " 
    << timer.current_time() << std::endl;
 
  if (dc.procid() == 0) 
    add_implicit_edges<edge_data>(implicitratingtype, graph, dc);

  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " 
    << timer.current_time() << std::endl;


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

  // Add error reporting to the engine
  const bool success = engine.add_edge_aggregator<error_aggregator>
    ("error", error_aggregator::map, error_aggregator::finalize) &&
    engine.aggregate_periodic("error", interval);
  ASSERT_TRUE(success);


  svdpp_vertex_program::GLOBAL_MEAN = graph.map_reduce_edges<double>(calc_global_mean);
  svdpp_vertex_program::NUM_TRAINING_EDGES = graph.map_reduce_edges<size_t>(count_edges);
  svdpp_vertex_program::GLOBAL_MEAN /= svdpp_vertex_program::NUM_TRAINING_EDGES;
  dc.cout() << "Global mean is: " <<svdpp_vertex_program::GLOBAL_MEAN << std::endl;

  // Signal all vertices on the vertices on the left (libersgd) 
  engine.map_reduce_vertices<graphlab::empty>(svdpp_vertex_program::signal_left);


  // Run the PageRank ---------------------------------------------------------
  dc.cout() << "Running Bias-SGD" << std::endl;
  dc.cout() << "(C) Code by Danny Bickson, CMU " << std::endl;
  dc.cout() << "Please send bug reports to danny.bickson@gmail.com" << std::endl;
  dc.cout() << "Time   Training    Validation" <<std::endl;
  dc.cout() << "       RMSE        RMSE " <<std::endl;
  timer.start();
  engine.start();  

  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
    << std::endl
    << "Final Runtime (seconds):   " << runtime 
                                        << std::endl
                                        << "Updates executed: " << engine.num_updates() << std::endl
                                        << "Update Rate (updates/second): " 
                                          << engine.num_updates() / runtime << std::endl;

  // Compute the final training error -----------------------------------------
  dc.cout() << "Final error: " << std::endl;
  engine.aggregate_now("error");

  // Make predictions ---------------------------------------------------------
  if(!predictions.empty()) {
    std::cout << "Saving predictions" << std::endl;
    const bool gzip_output = false;
    const bool save_vertices = false;
    const bool save_edges = true;
    const size_t threads_per_machine = 1;
    graph.save(predictions, prediction_saver(),
        gzip_output, save_vertices, 
        save_edges, threads_per_machine);
    //save the linear model
    graph.save(predictions + ".U", linear_model_saver_U(),
        gzip_output, save_edges, save_vertices, threads_per_machine);
    graph.save(predictions + ".V", linear_model_saver_V(),
        gzip_output, save_edges, save_vertices, threads_per_machine);
    graph.save(predictions + ".bias.U", linear_model_saver_bias_U(),
        gzip_output, save_edges, save_vertices, threads_per_machine);
    graph.save(predictions + ".bias.V", linear_model_saver_bias_V(),
        gzip_output, save_edges, save_vertices, threads_per_machine);

  }



  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



