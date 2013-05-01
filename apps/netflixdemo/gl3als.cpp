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
#include <graphlab.hpp>
#include <graphlab/engine/gl3engine.hpp>
#include <graphlab/util/random.hpp>
using namespace graphlab;

/**
 * Includes vertex_data, edge_data and helper functions for serialization
 */
#include "gl3als.hpp"
#include <graphlab/macros_def.hpp>

/**
 * \brief The engine type used by the ALS matrix factorization
 * algorithm.
 *
 * The ALS matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef gl3engine<graph_type> engine_type;



////////////////// Initialization Task ////////////////////////////////

// // Initialize the graph, sample 20% movies each user rated as test data
#define INIT_TASK 0
int init_map(const graph_type::vertex_type& center,
               graph_type::edge_type& edge,
               const graph_type::vertex_type& other) {
  // sample 20% of movies as validate set
  bool flip = random::fast_bernoulli(TEST_PERCENT);
  if (flip) {
    edge.data().role = edge_data::VALIDATE;
  }
  return 0;
}
int init_sum(int& v1, const int& v2) { return 0; }
void init_function(engine_type::context_type& context,
                   graph_type::vertex_type& vertex,
                   const engine_type::message_type& unused) {
  context.map_reduce<int>(INIT_TASK, OUT_EDGES);
}


////////////////// ALS Compute Task ////////////////////////////////
/**
 * \brief The gather type used to construct XtX and Xty needed for the ALS
 * update
 *
 * To compute the ALS update we need to compute the sum of
 * \code
 *  sum: XtX = nbr.factor.transpose() * nbr.factor
 *  sum: Xy  = nbr.factor * edge.obs
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
   * \brief Stores the current sum of nbr.factor.transpose() *
   * nbr.factor
   */
  mat_type XtX;

  /**
   * \brief Stores the current sum of nbr.factor * edge.obs
   */
  vec_type Xy;

  /** \brief basic default constructor */
  gather_type() { }

  /**
   * \brief This constructor computes XtX and Xy and stores the result
   * in XtX and Xy
   */
  gather_type(const vec_type& X, const double y) :
    XtX(X.size(), X.size()), Xy(X.size()) {
    XtX.triangularView<Eigen::Upper>() = X * X.transpose();
    Xy = X * y;
  } // end of constructor for gather type

  /** \brief Save the values to a binary archive */
  void save(graphlab::oarchive& arc) const { arc << XtX << Xy; }

  /** \brief Read the values from a binary archive */
  void load(graphlab::iarchive& arc) { arc >> XtX >> Xy; }

  /**
   * \brief Computes XtX += other.XtX and Xy += other.Xy updating this
   * tuples value
   */
  gather_type& operator+=(const gather_type& other) {
    if(other.Xy.size() == 0) {
      ASSERT_EQ(other.XtX.rows(), 0);
      ASSERT_EQ(other.XtX.cols(), 0);
    } else {
      if(Xy.size() == 0) {
        ASSERT_EQ(XtX.rows(), 0);
        ASSERT_EQ(XtX.cols(), 0);
        XtX = other.XtX; Xy = other.Xy;
      } else {
        XtX.triangularView<Eigen::Upper>() += other.XtX;
        Xy += other.Xy;
      }
    }
    return *this;
  } // end of operator+=
}; // end of gather type


#define ALS_MAP_REDUCE 1
/** The gather function computes XtX and Xy */
gather_type als_map(const graph_type::vertex_type& center,
                    graph_type::edge_type& edge,
                    const graph_type::vertex_type& other) {
  if(edge.data().role == edge_data::TRAIN) {
    return gather_type(other.data().factor, edge.data().obs);
  } else return gather_type();
} // end of gather function

gather_type als_sum(gather_type& v1, const gather_type& v2) {
  v1 += v2;
  return v1;
}

void als_update_function(engine_type::context_type& context,
                         graph_type::vertex_type& vertex,
                         const engine_type::message_type& unused) {
  // Get and reset the vertex data
  vertex_data& vdata = vertex.data();

  gather_type sum = context.map_reduce<gather_type> (ALS_MAP_REDUCE, ALL_EDGES);

  // if (vertex.id() == 1) {
  //   std::cout<< "vertex " << vertex.id() << ": " << std::endl
  //            << vdata.factor << std::endl;
  //   std::cout<< "sum: " << std::endl << sum.Xy << std::endl;
  // }

  // Determine the number of neighbors.  Each vertex has only in or
  // out edges depending on which side of the graph it is located
  if(sum.Xy.size() == 0) { vdata.residual = 0; ++vdata.nupdates; return; }
  mat_type XtX = sum.XtX;
  vec_type Xy = sum.Xy;

  // Add regularization
  double regularization = LAMBDA;
  if (REGNORMAL)
    regularization = LAMBDA*vertex.num_out_edges();

  for(int i = 0; i < XtX.rows(); ++i)
    XtX(i,i) += regularization;

  // Solve the least squares problem using eigen ----------------------------
  const vec_type old_factor = vdata.factor;
  vdata.factor = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xy);
  vdata.Wu = XtX;

  vdata.residual = (vdata.factor - old_factor).cwiseAbs().sum() / XtX.rows();
  ++vdata.nupdates;

  // no scatter yet
}

/////////////////////////// Collect Top Prediction //////////////////////
#define COLLECT_TASK 2
/// Collect, Prediction, Recommendation
typedef map_join<graphlab::vertex_id_type, double> map_id_rating_t;
typedef std::pair<map_id_rating_t, map_id_rating_t> map_join_pair;
map_join_pair collect_map (const graph_type::vertex_type& center,
                           graph_type::edge_type& edge,
                           const graph_type::vertex_type& other) {
  map_join_pair ret;
  if (edge.data().role == edge_data::TRAIN) {
    ret.first.data[other.id()] = edge.data().obs; // save the old rating
  } else {
    // use prediction
    double pred = center.data().factor.dot(other.data().factor);
    ret.second.data[other.id()] = pred; // save the prediction
  }
  return ret;
}

map_join_pair collect_sum (map_join_pair& v1, const map_join_pair& v2) {
  v1.first += v2.first;
  v1.second += v2.second;
  return v1;
}



std::string rank_list_to_string(const std::vector<std::pair<double, graphlab::vertex_id_type> >& ls) {
  std::stringstream sstream;
    for(size_t i = 0;i < ls.size(); ++i) {
      graphlab::vertex_id_type gid = ls[i].second;
      sstream << "\t" << id2movieid(gid)
              << ": " << mlist[gid]
              << " = " << ls[i].first << "\n";
    }
  return sstream.str();
}

void collect_function (engine_type::context_type& context,
                       graph_type::vertex_type& vertex) {
  if (is_user(vertex)) {
    map_join_pair sum = context.map_reduce<map_join_pair>(COLLECT_TASK, ALL_EDGES);
    vertex.data().top_rated = sum.first.get_top_k(10);
    vertex.data().top_pred = sum.second.get_top_k(10);
  }
}

class recommendation_writer {
 public:
   std::string save_vertex(const graph_type::vertex_type& v) {
    std::stringstream sstream;
    if (is_user(v)) {
      const std::vector<std::pair<double, graphlab::vertex_id_type> >& top_rated 
          = v.data().top_rated;
      const std::vector<std::pair<double, graphlab::vertex_id_type> >& top_pred
          = v.data().top_pred;
      if (top_rated.size() < 10 || top_pred.size() == 0) {
        return "";
      }
      // save top rated
      sstream << v.id() << " "; 
      sstream << pair2str(top_rated[0]); 
      for (size_t i = 1; i < top_rated.size(); ++i) {
        sstream << "," <<  (pair2str(top_rated[i]));
      }
      // save top pred
      sstream << " ";
      sstream << pair2str(top_pred[0]);
      for (size_t i = 1; i < top_pred.size(); ++i) {
        sstream << "," << (pair2str(top_pred[i]));
      }
      sstream << "\n";
      return sstream.str();
    } else {
      return "";
    }
   }
   std::string save_edge(const graph_type::edge_type& e) {
     return "";
   }

 private:
   std::string pair2str(const std::pair<double, graphlab::vertex_id_type>& pair) {
     int printingid = -pair.second - SAFE_NEG_OFFSET;
     return boost::lexical_cast<std::string>(printingid) + ":" + boost::lexical_cast<std::string>(pair.first);
   }
};

// the explanation module
#include "gl3als_explain.hpp"

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description =
      "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir;
  std::string predictions;
  size_t interval = 10;
  std::string movielist_dir;
  std::string saveprefix="result";
  clopts.attach_option("matrix", input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D",  vertex_data::NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("max_iter", MAX_ITER,
                       "The maxumum number of udpates allowed for a vertex");
  clopts.attach_option("lambda", LAMBDA,
                       "ALS regularization weight");
  clopts.attach_option("tol", TOLERANCE,
                       "residual termination threshold");
  clopts.attach_option("maxval", MAXVAL, "max allowed value");
  clopts.attach_option("minval", MINVAL, "min allowed value");
  clopts.attach_option("testpercent", TEST_PERCENT, "percentage of movies used for test");
  clopts.attach_option("saveprefix", saveprefix, "prefix for result files");
  clopts.attach_option("interval", interval,
                       "The time in seconds between error reports");
  clopts.attach_option("predictions", predictions,
                       "The prefix (folder and filename) to save predictions.");
  clopts.attach_option("regnormal", REGNORMAL,
                       "regularization type. 1 = weighted according to neighbors num. 0 = no weighting - just lambda");
  clopts.attach_option("movielist", movielist_dir,
                       "Movie List");
  parse_implicit_command_line(clopts);

  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }

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
  dc.barrier();
  graphlab::launch_metric_server();
  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer;
  graph_type graph(dc, clopts);
  graph.load(input_dir, graph_loader);
  dc.cout() << "Loading graph. Finished in "
            << timer.current_time() << std::endl;

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
  engine_type engine(dc, graph, clopts);

  engine.register_map_reduce(ALS_MAP_REDUCE,
                             als_map,
                             als_sum);
  engine.register_map_reduce(INIT_TASK,
                             init_map,
                             init_sum);
  engine.register_map_reduce(COLLECT_TASK,
                             collect_map,
                             collect_sum);
  engine.register_map_reduce(FACTOR_GATHER_TASK,
                             factor_gather,
                             factor_combine);

  vertex_set user_set = graph.select(is_user);
  vertex_set movie_set = graph.select(is_movie);

  // initialize graph. for each user random choose 20% movies as test 
  engine.set_vertex_program(init_function);
  engine.signal_vset(user_set); engine.wait();
  // start als
  engine.set_vertex_program(als_update_function);
  for (size_t i = 0; i < MAX_ITER; ++i) {
    timer.start();
    engine.signal_vset(user_set); engine.wait();
    engine.signal_vset(movie_set); engine.wait();
    const double runtime = timer.current_time();
    dc.cout() << "Complete iteration: " << i << std::endl;
    dc.cout() << "----------------------------------------------------------" 
              << std::endl
              << "Final Runtime (seconds):   " << runtime
              << std::endl
              << "Updates executed: " << engine.num_updates() << std::endl
              << "Update Rate (updates/second): "
              << engine.num_updates() / runtime << std::endl;
    // Compute the final training error -----------------------------------------
    error_aggregator errors = graph.map_reduce_edges<error_aggregator>(extract_l2_error);
    double l2norm = graph.map_reduce_vertices<double>(extract_l2_norm);
    dc.cout() << "Factor squared norm: " << l2norm/graph.num_vertices() << std::endl;
    dc.cout() << "Training RMSE: " << sqrt(errors.train/errors.ntrain) << std::endl;
    dc.cout() << "Test RMSE: " << sqrt(errors.test/errors.ntest) << std::endl;

    while(1) {
      int uid;
      if (dc.procid() == 0) {
        std::cout << "Enter User ID (-1:continue, -2:collect, -3:save):  ";
        std::cin >> uid;
      }
      dc.broadcast(uid, dc.procid() == 0);
      if (uid == -1) {
        break;
      } else if (uid == -2) {
        timer.start();
        dc.cout() << "Collecting prediction ..." << std::endl;
        engine.parfor_all_local_vertices(collect_function);
        engine.wait();
        dc.cout() << "Finish collecting prediction in " << timer.current_time() << " secs" << std::endl;

        timer.start();
        dc.cout() << "Collecting explanation..." << std::endl;
        engine.parfor_all_local_vertices(exp_collect_function);
        engine.wait();
        dc.cout() << "Finish collecting explanation in " << timer.current_time() << " secs" << std::endl;
        continue;
      } else if (uid == -3) {
        dc.cout() << "Save results to " << saveprefix << "..." << std::endl;
        graph.save(saveprefix, recommendation_writer(),
                   false, // no gzip
                   true, // save vertices
                   false // save edges
                   );
        dc.cout() << "Save explains to " << saveprefix+".explain" << "..." << std::endl;
        graph.save(saveprefix + ".explain", explain_writer(),
                   false,
                   true,
                   false);
        dc.cout() << "done" << std::endl;
        continue;
      }
      // does someone own uid?
      int tuid = graph.contains_vertex(uid);
      dc.all_reduce(tuid);
      if (tuid == 0) {
        if (dc.procid() == 0) {
          std::cout << "User " << uid << " does not exist\n";
        }
        continue;
      }
      bool is_master = graph.contains_vertex(uid) &&
          graph_type::local_vertex_type(graph.vertex(uid)).owned();
      if (is_master) {
        graph_type::vertex_type vtx(graph.vertex(uid));
        std::cout << " Top rated:\n" 
                  << rank_list_to_string(vtx.data().top_rated) << std::endl;
        std::cout << " Top pred:\n"
                  << rank_list_to_string(vtx.data().top_pred) << std::endl;
        std::cout << " Top explanation:\n"
                  << explain_list_to_string(vtx.data().top_pred, vtx.data().top_explain) << std::endl;
      }
    }
  }
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
  } // end of main
#include <graphlab/macros_undef.hpp>
  //   // Add error reporting to the engine
  //   const bool success = engine.add_edge_aggregator<error_aggregator>
  //     ("error", error_aggregator::map, error_aggregator::finalize) &&
  //     engine.aggregate_periodic("error", interval);
  //   ASSERT_TRUE(success);
  // 	
  //   size_t initial_max_updates = als_vertex_program::MAX_UPDATES;
  //   while(1) {
  //     // Signal all vertices on the vertices on the left (liberals)
  //     engine.map_reduce_vertices<graphlab::empty>(als_vertex_program::signal_left);
  //     info = graph.map_reduce_edges<stats_info>(count_edges);
  //     dc.cout()<<"Training edges: " << info.training_edges << " validation edges: " << info.validation_edges << std::endl;
  // 
  //     // Run ALS ---------------------------------------------------------
  //     dc.cout() << "Running ALS" << std::endl;
  //     timer.start();
  //     engine.start();
  // 
  //     const double runtime = timer.current_time();
  //     dc.cout() << "----------------------------------------------------------"
                   //               << std::endl
                                    //               << "Final Runtime (seconds):   " << runtime
                                                     //               << std::endl
                                                                      //               << "Updates executed: " << engine.num_updates() << std::endl
                                                                                       //               << "Update Rate (updates/second): "
                                                                                                        //               << engine.num_updates() / runtime << std::endl;
                                                                                                        // 
                                                                                                        //     // Compute the final training error -----------------------------------------
                                                                                                        //     dc.cout() << "Final error: " << std::endl;
                                                                                                        //     engine.aggregate_now("error");
                                                                                                        // 
                                                                                                        //     // Make predictions ---------------------------------------------------------
                                                                                                        //     if(!predictions.empty()) {
                                                                                                        //       std::cout << "Saving predictions" << std::endl;
                                                                                                        //       const bool gzip_output = false;
                                                                                                        //       const bool save_vertices = false;
                                                                                                        //       const bool save_edges = true;
                                                                                                        //       const size_t threads_per_machine = 2;
                                                                                                        // 
                                                                                                        //       //save the predictions
                                                                                                        //       graph.save(predictions, prediction_saver(),
                                                                                                        //                  gzip_output, save_vertices,
                                                                                                        //                  save_edges, threads_per_machine);
                                                                                                        //       //save the linear model
                                                                                                        //       graph.save(predictions + ".U", linear_model_saver_U(),
                                                                                                        //                  gzip_output, save_edges, save_vertices, threads_per_machine);
                                                                                                        //       graph.save(predictions + ".V", linear_model_saver_V(),
                                                                                                        //                  gzip_output, save_edges, save_vertices, threads_per_machine);
                                                                                                        // 
                                                                                                        //     }
                                                                                                        // 
                                                                                                        //   }
