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

/**
 * Includes vertex_data, edge_data and helper functions for serialization
 */

#include <graphlab/macros_def.hpp>
#include "main.hpp"
#include "io.hpp"
#include "compute_tasks.hpp"
#include "aggregate_tasks.hpp"

using namespace graphlab;

distributed_control dc;

void update_predictions(engine_type& engine, graph_type& graph, bool truncate) {
    //  Update cached prediction on edges
    engine.parfor_all_local_edges(boost::bind(compute_prediction, _1, _2, truncate)); engine.wait();
    {
      error_aggregator errors = graph.map_reduce_edges<error_aggregator>(extract_l2_error);
      dc.cout() << "Training RMSE: " << sqrt(errors.train/errors.ntrain) << std::endl;
      dc.cout() << "Test RMSE: " << sqrt(errors.test/errors.ntest) << std::endl;
    }
}


int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description =
      "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir;
  std::string saveprefix="result";
  std::string user_feature_dir;
  std::string rest_feature_dir;

  clopts.attach_option("matrix", input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("user_feature", user_feature_dir, "The directory containing the user feature file");
  clopts.attach_option("rest_feature", rest_feature_dir, "The directory containing the restaurant feature file");

  // clopts.attach_option("movielist", movielist_dir, "Movie List");
  clopts.attach_option("D",  vertex_data::NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("max_iter", MAX_ITER,
                       "The maxumum number of udpates allowed for a vertex");
  clopts.attach_option("interactive", INTERACTIVE,
                       "Use interactive session");
  clopts.attach_option("lambda", LAMBDA, "ALS regularization weight");
  clopts.attach_option("lambda2", LAMBDA2, "ALS regularization weight");
  clopts.attach_option("truncate", TRUNCATE, "Wether to truncate the prediction");
  clopts.attach_option("use_feature_latent", USE_FEATURE_LATENT, "Add feature latent vector to the model");
  clopts.attach_option("use_als", USE_BIAS_LATENT, "Add feature latent vector to the model");
  clopts.attach_option("use_bias", USE_BIAS, "Add global bias and user/rest bias to the model");
  clopts.attach_option("use_feature_weights", USE_FEATURE_WEIGHTS, "Add numerical feature weights to the model");
  clopts.attach_option("use_category_feature", USE_CATEGORY_FEATURE, "Add categorical feature weights to the model");
  clopts.attach_option("use_date_feature", USE_DATE_FEATURE, "Add date feature weights to the model");
  clopts.attach_option("use_numeric_feature", USE_NUMERIC_FEATURE, "Add numeric feature weights to the model");
  clopts.attach_option("minval", MIN_VAL, "MIN value for the prediction");
  clopts.attach_option("maxval", MAX_VAL, "MAX value for the prediction");
  clopts.attach_option("tol", TOLERANCE, "residual termination threshold");
  clopts.attach_option("testpercent", TEST_PERCENT, "percentage of movies used for test");
  clopts.attach_option("saveprefix", saveprefix, "prefix for result files");

  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  dc.barrier();
  graphlab::launch_metric_server();

  if (!USE_BIAS_LATENT) 
    vertex_data::NLATENT=0;

  dc.cout() << "Loading graph from " << input_dir << std::endl;
  graphlab::timer timer;
  graph_type graph(dc, clopts);
  if (!user_feature_dir.empty()) {
      dc.cout() << "Loading user feature from " << user_feature_dir << std::endl;
      graph.load(user_feature_dir, boost::bind(user_feature_loader, _1, _2, _3));
  } 
  if (!rest_feature_dir.empty()) {
      dc.cout() << "Loading restaurant feature from " << rest_feature_dir << std::endl;
      graph.load(rest_feature_dir, boost::bind(rest_feature_loader,_1,_2,_3));
  }

  graph.load(input_dir, graph_loader);

  dc.cout() << "Loading graph. Finished in " << timer.current_time() << std::endl;
  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " << timer.current_time() << std::endl;
  NEDGES = graph.num_edges();

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


  logstream(LOG_EMPH) << "Number of features =  " << feature_weights.size() << std::endl;

  dc.cout() << "Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts);

  // #define INIT_MAP_REDUCE 0
  engine.register_map_reduce(ALS_MAP_REDUCE,
                             als_map,
                             als_sum);
  engine.register_map_reduce(BIAS_MAP_REDUCE,
                             bias_map,
                             pair_sum<double, size_t>);
  vertex_set user_set = graph.select(is_user);
  vertex_set rest_set = graph.complete_set()-user_set; 

  // initialize graph. for each user random choose 20% movies as test 
  init_global_vars(dc);
  engine.parfor_all_local_vertices(init_vertex); 
  engine.wait();

  NTRAIN = graph.map_reduce_edges<int>(init_edge_hold_out); 
  NTEST = NEDGES - NTRAIN;

  int iter = 0;
  bool quit = false;
  while (!quit) {
    timer.start();
    double residual = 0;
    if (USE_BIAS_LATENT) {
      /**************************************************************************/
      /*                                                                        */
      /*                               start als                                */
      /*                                                                        */
      /**************************************************************************/
      logstream(LOG_EMPH) << "ALS on user/rest factors" << std::endl;
      engine.set_vertex_program(boost::bind(als_update_function, _1,_2,_3));
      engine.signal_vset(user_set); engine.wait();
      engine.signal_vset(rest_set); engine.wait();
      //  Update cached prediction on edges
      update_predictions(engine, graph, TRUNCATE);
    }

    if (USE_BIAS) {
      /**************************************************************************/
      /*                                                                        */
      /*                          compute global bias                           */
      /*                                                                        */
      /**************************************************************************/
      logstream(LOG_EMPH) << "Compute w0" << std::endl;
      double bias = graph.map_reduce_edges<double>(compute_bias);
      w0 = bias / NTRAIN;
      logstream(LOG_EMPH) << "w0 = " << w0 << std::endl;
      update_predictions(engine, graph, TRUNCATE);

      /**************************************************************************/
      /*                                                                        */
      /*                         compute local bias                             */
      /*                                                                        */
      /**************************************************************************/
      logstream(LOG_EMPH) << "Compute wu, wv" << std::endl;
      engine.parfor_all_local_vertices(compute_vertex_bias); engine.wait();
      update_predictions(engine, graph, TRUNCATE);
    } // end of use bias


    if (USE_FEATURE_WEIGHTS) {
      /**************************************************************************/
      /*                                                                        */
      /*                        compute feature weights                         */
      /*                                                                        */
      /**************************************************************************/
      logstream(LOG_EMPH) << "Compute w1,w2..." << std::endl;
      als_gather_type sum = graph.map_reduce_edges<als_gather_type>(regression_edge_map);
      logstream(LOG_EMPH) << "Finish gather" << std::endl;
      mat_type XtX = sum.XtX;
      vec_type Xy = sum.Xy;
      if (sum.Xy.size() > 0) { 
        for(int i = 0; i < XtX.rows(); ++i)
          XtX(i,i) += LAMBDA2;
        // Update weights 
        vec_type w = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xy);
        for (size_t i= 0; i < feature_keys.size(); ++i) {
          feature_weights[feature_keys[i]] = w[i];
        }
        logstream(LOG_EMPH) << "Finish compute" << std::endl;
        engine.parfor_all_local_vertices(init_vertex); 
        engine.parfor_all_local_edges(init_edge); 
        engine.wait();
        logstream(LOG_EMPH) << "Finish scatter" << std::endl;
        update_predictions(engine, graph, TRUNCATE);
      }
    } // end of update feature weights

    if (USE_FEATURE_LATENT) {
      /**************************************************************************/
      /*                                                                        */
      /*                     compute feature latent vectors                     */
      /*                                                                        */
      /**************************************************************************/
      typedef boost::unordered_map<size_t, vec_type>::value_type kv_type;
      // coordinate descent for each feature
      foreach(kv_type& kv, feature_latent) {
        logstream(LOG_EMPH) << "Update global feature " << kv.first << std::endl;
        ASSERT_TRUE(kv.first != 0);
        size_t fid = kv.first;

        // Gather sufficient statistics on vertices
        als_gather_type sum = graph.map_reduce_edges<als_gather_type> (boost::bind(als_map_edge, _1, fid));
        mat_type XtX = sum.XtX;
        vec_type Xy = sum.Xy;
        if (sum.Xy.size() == 0) { 
          for(int i = 0; i < XtX.rows(); ++i)
            XtX(i,i) += LAMBDA2;
          vec_type old = kv.second;

          // Update coordinate
          kv.second = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xy);
          residual += (kv.second - old).cwiseAbs().sum() / XtX.rows();
          logstream(LOG_EMPH) << "feature " <<  kv.first << "\n" 
                              << "norm = " << kv.second.norm() << std::endl;

          //  Update cached prediction on vertices
          engine.parfor_all_local_vertices(init_vertex); 
          engine.parfor_all_local_edges(init_edge); 
          engine.wait();

          //  Update cached prediction on edges
          update_predictions(engine, graph, TRUNCATE);
        } // end of update
      } // end of foreach feature
    } // end of update feature feature latent 

    residual += graph.map_reduce_vertices<double>(extract_residual); 
    dc.cout() << "ALS residual: " << residual << std::endl;

    const double runtime = timer.current_time();
    dc.cout() << "Complete iteration: " << iter << std::endl;
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

    if (iter < MAX_ITER) {
      ++iter;
      continue;
    } else {
      // ++iter;
      break;
    }
  }
  logstream(LOG_EMPH) << "Finish." << std::endl;
    graphlab::mpi_tools::finalize();
    return EXIT_SUCCESS;
  } // end of main
#include <graphlab/macros_undef.hpp>

  // if (!movielist_dir.empty()) {
  //   dc.cout() << "Loading movie names from " << movielist_dir << std::endl; 
  //   std::ifstream fin(movielist_dir.c_str());
  //   size_t id = 1;
  //   while(fin.good()) {
  //     std::string name;
  //     std::getline(fin, name);
  //     graphlab::vertex_id_type gid = -(graphlab::vertex_id_type(id + SAFE_NEG_OFFSET));
  //     mv_names[gid] = name;
  //     id++;
  //   }
  //   fin.close();
  // }


    // // begin iteractive session
    // if (INTERACTIVE) {
    //   while(1) {
    //     int uid;
    //     if (dc.procid() == 0) {
    //       std::cout << "Enter User ID (-1:continue, -2:collect, -3:save):  ";
    //       std::cin >> uid;
    //     }
    //     dc.broadcast(uid, dc.procid() == 0);
    //     if (uid == -1) {
    //       break;
    //     } else if (uid == -2) {
    //       // timer.start();
    //       // dc.cout() << "Collecting prediction ..." << std::endl;
    //       // engine.parfor_all_local_vertices(collect_function);
    //       // engine.wait();
    //       // dc.cout() << "Finish collecting prediction in " << timer.current_time() << " secs" << std::endl;

    //       // timer.start();
    //       // dc.cout() << "Collecting explanation..." << std::endl;
    //       // engine.parfor_all_local_vertices(exp_collect_function);
    //       // engine.wait();
    //       // dc.cout() << "Finish collecting explanation in " << timer.current_time() << " secs" << std::endl;
    //       // continue;
    //     } else if (uid == -3) {
    //       dc.cout() << "Save results to " << saveprefix << "..." << std::endl;
    //       // graph.save(saveprefix+".recommend", recommendation_writer(),
    //       //            false, // no gzip
    //       //            true, // save vertices
    //       //            false // save edges
    //       //            );
    //       dc.cout() << "done" << std::endl;
    //       continue;
    //     } else if (uid == -4) {
    //       dc.cout() << "Save model to " << saveprefix << "..." << std::endl;
    //       graph.save(saveprefix + ".user", linear_model_saver_U(),
    //                  false, true, false);
    //       graph.save(saveprefix + ".movie", linear_model_saver_V(),
    //                  false, true, false);
    //       continue;
    //     } else if (uid == -5) {
    //       quit = true;
    //       break;
    //     }
    //     // does someone own uid?
    //     int tuid = graph.contains_vertex(uid);
    //     dc.all_reduce(tuid);
    //     if (tuid == 0) {
    //       if (dc.procid() == 0) {
    //         std::cout << "User " << uid << " does not exist\n";
    //       }
    //       continue;
    //     }
    //     bool is_master = graph.contains_vertex(uid) &&
    //         graph_type::local_vertex_type(graph.vertex(uid)).owned();
    //     // if (is_master) {
    //       // graph_type::vertex_type vtx(graph.vertex(uid));
    //       // std::cout << " Top rated:\n" 
    //       //           << rank_list_to_string(vtx.data().top_rated) << std::endl;
    //       // std::cout << " Top pred:\n"
    //       //           << rank_list_to_string(vtx.data().top_pred) << std::endl;
    //     // }
    //   }
    // } else {
    //   timer.start();
    //   // dc.cout() << "Collecting prediction ..." << std::endl;
    //   // engine.parfor_all_local_vertices(collect_function);
    //   // engine.wait();
    //   // dc.cout() << "Finish collecting prediction in " << timer.current_time() << " secs" << std::endl;
    //   // dc.cout() << "done" << std::endl;
    //   // dc.cout() << "Save results to " << saveprefix << "..." << std::endl;
    //   // graph.save(saveprefix+".recommend", recommendation_writer(),
    //   //            false, // no gzip
    //   //            true, // save vertices
    //   //            false // save edges);
    //   dc.cout() << "done" << std::endl;
    //   dc.cout() << "Save model to " << saveprefix << "..." << std::endl;
    //   graph.save(saveprefix + ".user", linear_model_saver_U(),
    //              false, true, false);
    //   graph.save(saveprefix + ".movie", linear_model_saver_V(),
    //              false, true, false);
    //   dc.cout() << "done" << std::endl;
    //   quit = true;
    // }

