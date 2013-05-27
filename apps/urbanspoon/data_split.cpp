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
#include "data_split.hpp"
#include <graphlab/macros_def.hpp>

using namespace graphlab;

distributed_control dc;
double subsample = 1;
double split_ratio = 0.2;
bool timeaware = true;

bool compare_edge_by_date_fun(const std::pair<graph_type::vertex_id_type,edge_data>& v1,
                              const std::pair<graph_type::vertex_id_type,edge_data>& v2) {
  return v1.second.d < v2.second.d;
}

#define VSPLIT_MAP_REDUCE 0
typedef std::vector< std::pair<graph_type::vertex_id_type, edge_data> > edge_vec_type;
edge_vec_type vsplit_map_fun(const graph_type::vertex_type& center,
                   graph_type::edge_type& edge,
                   const graph_type::vertex_type& other) {
  edge_vec_type ret;
  ret.push_back(std::pair<graph_type::vertex_id_type, edge_data> (other.id(), edge.data()));
  return ret;
}

edge_vec_type vsplit_sum_fun(edge_vec_type& v1, const edge_vec_type& v2) {
  if (timeaware) {
    edge_vec_type ret(v1.size() + v2.size());
    std::merge(v1.begin(), v1.end(),
               v2.begin(), v2.end(),
               ret.begin(), compare_edge_by_date_fun);
    v1 = ret;
  } else {
    v1.insert(v1.end(), v2.begin(), v2.end());
  }
  return v1;
}

void vsplit_update_function(engine_type::context_type& context,
                                graph_type::vertex_type& vertex,
                                graph_type* training_graph,
                                graph_type* testing_graph) {
  edge_vec_type sum = context.map_reduce<edge_vec_type>(VSPLIT_MAP_REDUCE, OUT_EDGES);
 
  // split the edgeset into training/testing. 
  int split = sum.size() * (1-split_ratio);
  for (size_t i = 0; i < split; ++i) {
      // use subsample
      bool flip = random::fast_bernoulli(subsample);
      if (flip) {
        std::pair<graph_type::vertex_id_type, edge_data> val = sum[i];
        training_graph->add_edge(vertex.id(), val.first, val.second);
      }
  }
  for (size_t i = split; i < sum.size(); ++i) {
    bool flip = random::fast_bernoulli(subsample);
    if (flip) {
      std::pair<graph_type::vertex_id_type, edge_data> val = sum[i];
      testing_graph->add_edge(vertex.id(), val.first, val.second);
    }
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
  clopts.attach_option("saveprefix", saveprefix, "prefix for result files");
  clopts.attach_option("subsample", subsample, "subsampling the raw data");
  clopts.attach_option("timeaware", timeaware, "subsampling the raw data");
  clopts.attach_option("split_ratio", split_ratio, "ratio for training/testing splits");
  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  dc.barrier();
  graphlab::launch_metric_server();
  dc.cout() << "Loading graph from " << input_dir << std::endl;
  graphlab::timer timer;
  graph_type graph(dc, clopts);

  // Load all edges
  graph.load(input_dir, graph_loader);
  dc.cout() << "Loading graph. Finished in " << timer.current_time() << std::endl;
  dc.cout() << "Finalizing graph." << std::endl;

  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " << timer.current_time() << std::endl;
  // Split into training and test files
  graph_type training_graph (dc, clopts);
  graph_type testing_graph (dc, clopts);


  // Create engine
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

  engine.register_map_reduce(VSPLIT_MAP_REDUCE,
                             vsplit_map_fun,
                             vsplit_sum_fun);

  logstream(LOG_EMPH) << "Splitting graph into training and testing set\n"
                      << "Split ratio = " << split_ratio << "\n"
                      << "Subsample ratio = " << subsample << "\n"
                      << "Time aware = " << timeaware << std::endl;

  engine.parfor_all_local_vertices(
      boost::bind(vsplit_update_function, _1, _2,
                                         &training_graph,
                                         &testing_graph));

  engine.wait();
  graph.clear();
  // Finalize and save the graph
  training_graph.finalize();
  testing_graph.finalize();

  size_t NTRAINING = training_graph.num_edges();
  size_t NTESTING = testing_graph.num_edges();
  logstream(LOG_EMPH) << "training_graph.num_edges() = " << NTRAINING << std::endl;
  logstream(LOG_EMPH) << "testing_graph.num_edges() = " << NTESTING << std::endl;

  training_graph.save(saveprefix + ".training",
                      graph_edge_writer(),
                      false, // no gzip
                      false, // save vertices
                      true); // save edges
  testing_graph.save(saveprefix + ".testing",
                     graph_edge_writer(),
                     false,
                     false,
                     true);
  logstream(LOG_EMPH) << "Finish." << std::endl;
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main
#include <graphlab/macros_undef.hpp>
