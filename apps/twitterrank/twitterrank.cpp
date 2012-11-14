#include <graphlab.hpp>
#include <graphlab/graph/graph_vertex_join.hpp>
#include <boost/functional/hash.hpp>
#include "pagerank.hpp"
#include "cgs_lda.hpp"

size_t NTOPICS = 20;
graphlab::distributed_control* dc_ptr;
graphlab::thread_group thgroup;
graphlab::omni_engine<pagerank::compute_transit_prob>* transit_prob_engine; 
graphlab::omni_engine<pagerank::compute_pagerank>* pagerank_engine; 
graphlab::omni_engine<lda::cgs_lda_vertex_program>* lda_engine; 
pagerank::graph_type* pagerank_graph;
lda::graph_type* lda_graph;
graphlab::graph_vertex_join<pagerank::graph_type, lda::graph_type>* vjoinptr;

/*  Global thread func for running pagerank */
void fn_run_pagerank () {
  static graphlab::mutex m;
  if (m.try_lock()) {
    pagerank_engine->signal_all();
    pagerank_engine->start();
    m.unlock();
  }
  // const float runtime = engine.elapsed_seconds();
  // dc.cout() << "Finished Running engine in " << runtime << " seconds." << std::endl;
}

/* Global thread func for running lda */
void fn_run_lda() {
  static graphlab::mutex m;
  if (m.try_lock()) {
    lda_engine->signal_all();
    lda_engine->start();
    m.unlock();
  }
}

/* Global thread func for join vertices and recompute transit probs */
size_t left_emit_key (const pagerank::graph_type::vertex_type& vertex) {
  return vertex.data().join_key;
}
size_t right_emit_key (const lda::graph_type::vertex_type& vertex) {
  return vertex.data().join_key;
}
void vertex_join_op (pagerank::graph_type::vertex_type& lvertex,
                     const lda::graph_type::vertex_data_type& rvertex_data) {
  lvertex.data().topics.clear();
  for (size_t i = 0; i < rvertex_data.factor.size(); i++) {
    lvertex.data().topics[i] = rvertex_data.factor[i];
  }
}
void fn_compute_transit_prob () {
  static graphlab::mutex m;
  if (m.try_lock()) {
    transit_prob_engine->signal_all();
    transit_prob_engine->start();
    m.unlock();
  }
}
void fn_join_vertex() {
  vjoinptr->left_injective_join(vertex_join_op);
  fn_compute_transit_prob();
}

void launch_server() {
  graphlab::launch_metric_server();

  // Start the lda webserver 
  graphlab::add_metric_server_callback("wordclouds", lda::word_cloud_callback);

  // Start the pagerank webserver 
  graphlab::add_metric_server_callback("pagerank", pagerank::pagerank_callback);
}

void load_pagerankgraph(const std::string& edge_dir, const std::string& vertex_dir) {
  pagerank::load_and_initialize_graph(*dc_ptr, *pagerank_graph, edge_dir, vertex_dir);
}

void load_ldagraph(const std::string& edge_dir,
                   const std::string& vertex_dir,
                   const std::string& dictionary_dir) {
  lda::load_and_initialize_graph(*dc_ptr, *lda_graph, edge_dir, vertex_dir);
  lda::load_dictionary(dictionary_dir);
}

int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  dc_ptr = &dc;
  global_logger().set_log_level(LOG_INFO);
  
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Twitter Rank.");
  std::string pagerank_edges;
  std::string pagerank_vertices;
  std::string lda_edges;
  std::string lda_vertices;
  std::string lda_dictionary;
  std::string execution_type = "synchronous";

  // The dir of the link graph 
  clopts.attach_option("pagerank_edges", pagerank_edges, "The pagerank graph (edges). Required ");
  clopts.attach_option("pagerank_vertices", pagerank_vertices, "The pagerank graph (vertices). Required ");
  clopts.attach_option("lda_edges", lda_edges, "The lda graph (edges). Required ");
  clopts.attach_option("lda_vertices", lda_vertices, "The lda graph (vertices). Required ");
  clopts.attach_option("lda_dictionary", lda_dictionary, "The lda word dictionary. Required ");
  clopts.attach_option("execution", execution_type, "Execution type (synchronous or asynchronous)");
  clopts.attach_option("ntopics", NTOPICS, "Number of topics to use");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (pagerank_edges == "" || pagerank_vertices == "") {
    dc.cout() << "Pagerank Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }
  if (lda_edges == "" || lda_vertices == "") {
    dc.cout() << "LDA Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }
  if (lda_dictionary == "") {
    dc.cout() << "LDA Dictionary not specified. Cannot continue";
    return EXIT_FAILURE;
  }
  pagerank::NTOPICS = NTOPICS;
  lda::NTOPICS = NTOPICS;

  // Build the pagerank (left) graph ----------------------------------------------------------
  // The global pagerank_graph points to lgraph 
  pagerank::graph_type lgraph(dc);
  pagerank_graph = &lgraph;
  load_pagerankgraph(pagerank_edges, pagerank_vertices);

  // Set up engine1 to compute transition probability
  graphlab::omni_engine<pagerank::compute_transit_prob> engine0(dc, lgraph, execution_type);
  transit_prob_engine = &engine0;
  fn_compute_transit_prob();

  // Build the lda (right) graph
  // The global lda_graph points to rgraph
  lda::graph_type rgraph(dc);
  lda_graph = &rgraph;
  load_ldagraph(lda_edges, lda_vertices, lda_dictionary);

  // create the graph join object
  // prepare the join
  // assign the global pointer to the join object
  graphlab::graph_vertex_join<pagerank::graph_type, lda::graph_type> vjoin(dc, lgraph, rgraph);
  vjoin.prepare_injective_join(left_emit_key, right_emit_key);
  vjoinptr = &vjoin;

 
  // Run pagerank ------------------------------------------------------------
  //
  graphlab::omni_engine<pagerank::compute_pagerank> engine1(dc, lgraph, execution_type); 
  pagerank_engine = &engine1;
  thgroup.launch(fn_run_pagerank);
  
  // Run lda -----------------------------------------------------------------
  graphlab::omni_engine<lda::cgs_lda_vertex_program> engine2(dc, rgraph, execution_type);
  lda_engine = &engine2;
  thgroup.launch(fn_run_lda);
  
  // Run background vertex join
  thgroup.launch(fn_join_vertex);


  // --------------------------------------
  // --------------------------------------
  // --------------------------------------


  graphlab::stop_metric_server_on_eof();
  graphlab::mpi_tools::finalize();


  // Run pagerank ------------------------------------------------------------
  return EXIT_SUCCESS;
}
