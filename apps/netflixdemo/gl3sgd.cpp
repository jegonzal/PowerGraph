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
 *
 * \brief The main file for the ALS matrix factorization algorithm.
 *
 * This file contains the main body of the ALS matrix factorization
 * algorithm.
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>
#include <graphlab/engine/gl3engine.hpp>
#include <Eigen/Dense>
#include "eigen_serialization.hpp"
#include <graphlab/macros_def.hpp>

#define VERTEX_DELTA_TASK_ID 0

typedef Eigen::VectorXd vec_type;
typedef Eigen::MatrixXd mat_type;

//when using negative node id range, we are not allowed to use
//0 and 1 so we add 2.
const static int SAFE_NEG_OFFSET=2;
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
  vec_type pvec;

  int nupdates;
  /**
   * \brief Simple default constructor which randomizes the vertex
   *  data
   */
  vertex_data() : nupdates(0) { if (debug) pvec = vec_type::Ones(NLATENT); else randomize(); }
  /** \brief Randomizes the latent pvec */
  void randomize() { pvec.resize(NLATENT); pvec.setRandom(); }
  /** \brief Save the vertex data to a binary archive */
  void save(graphlab::oarchive& arc) const {
    arc << nupdates << pvec;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) {
    arc >> nupdates >> pvec;
  }
}; // end of vertex data

std::size_t hash_value(vertex_data const& b) {
  return b.nupdates;
}


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data sgdo stores the most recent error estimate.
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

std::size_t hash_value(edge_data const& b) {
  return boost::hash_value(b.obs);
}


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::gl3engine<graph_type> engine_type;


/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph,
                         const std::string& filename,
                         const std::string& line) {
  if (boost::starts_with(line, "#") || boost::starts_with(line, "%")) {
    logstream(LOG_INFO) << line << std::endl;
    return true;
  }
  ASSERT_FALSE(line.empty());
  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
  else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;
  // Parse the line
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  strm >> source_id >> target_id;

  // for test files (.predict) no need to read the actual rating value.
  if(role == edge_data::TRAIN || role == edge_data::VALIDATE) {
    strm >> obs;
  }
  target_id = -(graphlab::vertex_id_type(target_id + SAFE_NEG_OFFSET));

  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(obs, role));
  return true; // successful load
} // end of graph_loader

double LAMBDA = 0.001;
double GAMMA = 0.001;


void vertex_delta(graph_type::vertex_type& vtx, const vec_type& delta) {
  vtx.data().pvec += delta;
  vtx.data().nupdates++;
}

void sgd_function(engine_type::context_type& context,
                  graph_type::edge_type& edge) {
  double pred = edge.source().data().pvec.dot(edge.target().data().pvec);

  const float err = edge.data().obs - pred;
  vec_type delta;
  delta = GAMMA*(err*edge.source().data().pvec - LAMBDA*edge.target().data().pvec);
  context.send_delta(VERTEX_DELTA_TASK_ID, edge.target(), delta);


  delta = GAMMA*(err*edge.target().data().pvec - LAMBDA*edge.source().data().pvec);
  context.send_delta(VERTEX_DELTA_TASK_ID, edge.source(), delta);
}


void sync_function(engine_type::context_type& context,
                   graph_type::vertex_type& vertex) {
  context.synchronize(vertex);
}

/**
 * \brief Given an edge compute the error associated with that edge
 */
double extract_l2_error(const graph_type::edge_type & edge) {
  double pred =
      edge.source().data().pvec.dot(edge.target().data().pvec);
  double rmse = (edge.data().obs - pred) * (edge.data().obs - pred);
  return rmse;
} // end of extract_l2_error



size_t vertex_data::NLATENT = 20;

struct map_join {
  boost::unordered_map<graphlab::vertex_id_type, double> data; // rating than movie
  void save(graphlab::oarchive& oarc) const {
    oarc << data;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> data;
  }

  map_join& operator+=(const map_join& other) {
    std::copy(other.data.begin(), other.data.end(),
              std::inserter(data, data.end()));
    return *this;
  }

  std::vector<std::pair<double, graphlab::vertex_id_type> > get_top_k(size_t n) {
    std::vector<std::pair<double, graphlab::vertex_id_type> > ret;
    boost::unordered_map<graphlab::vertex_id_type, double>::const_iterator iter = data.begin();
    std::vector<std::pair<double, graphlab::vertex_id_type> > all_copy;
    while (iter != data.end()) {
      all_copy.push_back(std::make_pair(iter->second, iter->first));
      ++iter;
    }
    std::sort(all_copy.rbegin(), all_copy.rend());
    size_t limit = all_copy.size() < 10 ? all_copy.size() : 10;
    std::copy(all_copy.begin(), all_copy.begin() + limit,
              std::inserter(ret, ret.end()));
    return ret;
  }
};


int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description =
      "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir;
  std::string movielist_dir;
  size_t interval = 0;
  size_t ITERATIONS = 10;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D", vertex_data::NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("interval", interval,
                       "The time in seconds between error reports");
  clopts.attach_option("lambda", LAMBDA,
                       "SGD regularization weight");
  clopts.attach_option("gamma", GAMMA,
                       "SGD step size");
  clopts.attach_option("iterations", ITERATIONS,
                       "number of SGD iterations");
  clopts.attach_option("movielist", movielist_dir,
            "Movie List");
  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  boost::unordered_map<graph_type::vertex_id_type, std::string> mlist;
  if (!movielist_dir.empty()) {
    std::ifstream fin(movielist_dir.c_str());
    size_t id = 1;
    while(fin.good()) {
      std::string name;
      std::getline(fin, name);
      graphlab::vertex_id_type gid = -(graphlab::vertex_id_type(id + SAFE_NEG_OFFSET));
      mlist[gid] = name;
      id++;
    }
    fin.close();
  }
  ///! Initialize control plain using mpi
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


  engine_type engine(dc, graph, clopts);

  engine.register_vertex_delta<vec_type>(VERTEX_DELTA_TASK_ID,
                                         vertex_delta);

  dc.cout() << "Running SGD" << std::endl;

  timer.start();
  for (size_t i = 0;i < ITERATIONS; ++i) {
    engine.parfor_all_local_edges(sgd_function);
    engine.parfor_all_local_vertices(sync_function);
    engine.wait();
    double rmse = graph.map_reduce_edges<double>(extract_l2_error);
    dc.cout() << "RMSE = " << sqrt(rmse / graph.num_edges()) << std::endl;
    GAMMA = GAMMA * 0.9;
  }




  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
            << std::endl
            << "Final Runtime (seconds):   " << runtime;

  // Compute the final training error -----------------------------------------
  dc.cout() << "Final error: " << std::endl;

  while(1) {
    int uid;
    if (dc.procid() == 0) {
      std::cout << "Enter User ID (-1) to quit: " ;
      std::cin >> uid;
    }
    dc.broadcast(uid, dc.procid() == 0);
    if (uid == -1) break;
    // does someone own uid?
    int tuid = graph.contains_vertex(uid);
    dc.all_reduce(tuid);
    if (tuid == 0) {
      if (dc.procid() == 0) {
        std::cout << "User " << uid << " does not exist\n";
      }
      continue;
    }
    // every search for user ID uid.
    map_join all_training;
    if (graph.contains_vertex(uid)) {
      graph_type::vertex_type vtx(graph.vertex(uid));
      graph_type::local_vertex_type lvtx(vtx);
      foreach(graph_type::local_edge_type edge, lvtx.out_edges()) {
        graph_type::local_vertex_type target = edge.target();
        graph_type::vertex_id_type gid = target.global_id();
        all_training.data[gid] = edge.data().obs;
      }
    }
    dc.all_reduce(all_training);
    // print the training data
    if (dc.procid() == 0) {
      std::cout << "Top 10 rated movies:\n";
      std::vector<std::pair<double, graphlab::vertex_id_type> > top10 = all_training.get_top_k(10) ;
      for(size_t i = 0;i < top10.size(); ++i) {
        graphlab::vertex_id_type gid = top10[i].second;
        int printingid = - gid - SAFE_NEG_OFFSET;
        std::cout << "\t" << printingid;
        if (mlist.find(gid) != mlist.end()) {
          std::cout << ": " << mlist[gid];
        }
        std::cout << " = " << top10[i].first << "\n";
      }
    }


    map_join all_predict;
    // now for the recommendations.
    bool is_master = graph.contains_vertex(uid) &&
        graph_type::local_vertex_type(graph.vertex(uid)).owned();
    // broadcast the user vector
    vec_type factor;
    if (is_master) {
      factor = graph.vertex(uid).data().pvec;
    }
    dc.broadcast(factor, is_master);
    // now loop through all the vertices.
    for (size_t i = 0;i < graph.num_local_vertices(); ++i) {
      graph_type::local_vertex_type lvtx(graph.l_vertex(i));
      if (lvtx.owned() && (int)(lvtx.global_id()) < 0) {
        double pred = lvtx.data().pvec.dot(factor);
        all_predict.data[lvtx.global_id()] = pred;
      }
    }
    dc.all_reduce(all_predict);
    if (dc.procid() == 0) {
      std::cout << "Top 10 predicted movies:\n";
      std::vector<std::pair<double, graphlab::vertex_id_type> > top10 = all_predict.get_top_k(10) ;
      for(size_t i = 0;i < top10.size(); ++i) {
        graphlab::vertex_id_type gid = top10[i].second;
        int printingid = - gid - SAFE_NEG_OFFSET;
        std::cout << "\t" << printingid;
        if (mlist.find(gid) != mlist.end()) {
          std::cout << ": " << mlist[gid];
        }
        std::cout << " = " << top10[i].first << "\n";
      }
    }
  }
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



