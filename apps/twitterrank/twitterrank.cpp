// #include <graphlab.hpp>
#include <graphlab/graph/graph_vertex_join.hpp>
#include <boost/functional/hash.hpp>
#include "pagerank.hpp"
#include "cgs_lda.hpp"

size_t NTOPICS = 20;
int JOIN_INTERVAL = 5;
bool JOIN_ON_ID = true;
int RUN_PAGERANK = 1;
int RUN_LDA = 2;



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
    dc_ptr->cout() << "Running pagerank" << std::endl;
    pagerank_engine->signal_all();
    pagerank_engine->start();
    dc_ptr->cout() << "Pagerank finished" << std::endl;
    m.unlock();
  }
  // const float runtime = engine.elapsed_seconds();
  // dc.cout() << "Finished Running engine in " << runtime << " seconds." << std::endl;
}

/* Global thread func for running lda */
void fn_run_lda() {
  static graphlab::mutex m;
  if (m.try_lock()) {
    dc_ptr->cout() << "Running The Collapsed Gibbs Sampler" << std::endl;
    lda_engine->map_reduce_vertices<graphlab::empty>(lda::signal_only::docs);
    // Enable sampling
    lda::cgs_lda_vertex_program::DISABLE_SAMPLING = false;
    // Run the engine
    lda_engine->start();
    // Finalize the counts
    // lda::cgs_lda_vertex_program::DISABLE_SAMPLING = true;
    // lda_engine->signal_all();
    // lda_engine->start();
    m.unlock();
  }
}


void vertex_join_op (pagerank::graph_type::vertex_type& lvertex,
                     const lda::graph_type::vertex_data_type& rvertex_data) {
  // lvertex.data().topics.clear();
  for (size_t i = 0; i < rvertex_data.factor.size(); i++) {
    // lvertex.data().topics.push_back(rvertex_data.factor[i]);
    lvertex.data().topics[i] = rvertex_data.factor[i];
  }
}
void fn_compute_transit_prob () {
  static graphlab::mutex m;
  if (m.try_lock()) {
    dc_ptr->cout() << "Compute transition prob" << std::endl;
    transit_prob_engine->signal_all();
    transit_prob_engine->start();
    m.unlock();
  }
}
void fn_join_vertex() {
  while(true) {
    dc_ptr->cout() << "Join vertex data" << std::endl;
    vjoinptr->left_injective_join(vertex_join_op);
    fn_compute_transit_prob();
    graphlab::timer::sleep(JOIN_INTERVAL);
  }
}

void launch_metric_server(int algorithm) {
  graphlab::launch_metric_server();

  // Start the pagerank webserver 
  if (algorithm & RUN_PAGERANK) 
    graphlab::add_metric_server_callback("pagerank", pagerank::pagerank_callback);

  // Start the lda webserver 
  if (algorithm & RUN_LDA) 
  graphlab::add_metric_server_callback("wordclouds", lda::word_cloud_callback);
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
  clopts.attach_option("wordid_offset", lda::WORDID_OFFSET, "The starting id of the words in the lda input.");
  clopts.attach_option("execution", execution_type, "Execution type (synchronous or asynchronous)");
  clopts.attach_option("ntopics", NTOPICS, "Number of topics to use");
  clopts.attach_option("topkpr", pagerank::TOPK, "Top k pages(users) for display");
  clopts.attach_option("topklda", lda::TOPK, "Top k words in each topic for display");
  clopts.attach_option("default_ndocs", pagerank::DEFAULT_NDOCS, "Top k pages(users) to display");
  clopts.attach_option("default_topicval", pagerank::DEFAULT_TOPICVAL, "Top k pages(users) to display");
  clopts.attach_option("join_on_id", JOIN_ON_ID, "If true, use the vertex id as join key. Otherwise, use vertex data as join key, so the vertex input for both graph must be provided.");
  clopts.attach_option("has_doc_count", pagerank::HAS_DOC_COUNT, "Whether or not the pagerank vertex data has a field for document count. This could be true for author-topic graph, depending on the input data.");

  bool pagerank_only = false;
  bool lda_only = false;
  clopts.attach_option("pagerank_only", pagerank_only, "Run pagerank only.");
  clopts.attach_option("lda_only", lda_only, "Run lda only.");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  if (pagerank_only && lda_only) { 
    dc.cout () << "pagerank_only and lda_only are mutually exclusive";
    return EXIT_FAILURE;
  }

  int algorithm = RUN_PAGERANK | RUN_LDA;
  if (pagerank_only) {
    dc.cout () << " Running pagerank only. " << std::endl;
    algorithm &= RUN_PAGERANK; 
  } else if (lda_only) {
    dc.cout () << " Running LDA only. " << std::endl;
    algorithm &= RUN_LDA; 
  } else {
    dc.cout () << " Running Joint Pagerank and LDA. " << std::endl;
  }

  if (!lda_only && pagerank_edges == "") {
    dc.cout() << "Pagerank Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }
  if (!pagerank_only && lda_edges == "") {
    dc.cout() << "LDA Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }

  if ((algorithm & RUN_PAGERANK) && (algorithm & RUN_LDA) && !JOIN_ON_ID && (pagerank_vertices == "" || lda_vertices == "")) {
    dc.cout() << "JOIN_ON_ID is false, vertex input for both graph must be provided. Please provide pagerank_vertices and lda_vertices.";
    return EXIT_FAILURE;
      }
  if ((algorithm & RUN_LDA) && lda_dictionary == "") {
    dc.cout() << "LDA Dictionary not specified. Cannot continue";
    return EXIT_FAILURE;
  }
  pagerank::NTOPICS = NTOPICS;
  lda::NTOPICS = NTOPICS;
  pagerank::JOIN_ON_ID = JOIN_ON_ID;
  lda::JOIN_ON_ID = JOIN_ON_ID;

  // Build the pagerank (left) graph ----------------------------------------------------------
  // The global pagerank_graph points to lgraph 
  pagerank::graph_type lgraph(dc);
  pagerank_graph = &lgraph;
  if (algorithm & RUN_PAGERANK) {
    load_pagerankgraph(pagerank_edges, pagerank_vertices);
  }

  // Build the lda (right) graph
  // The global lda_graph points to rgraph
  lda::graph_type rgraph(dc);
  lda_graph = &rgraph;
  if (algorithm & RUN_LDA) {
    load_ldagraph(lda_edges, lda_vertices, lda_dictionary);
  }

  // create the graph join object
  // prepare the join
  // assign the global pointer to the join object
  graphlab::graph_vertex_join<pagerank::graph_type, lda::graph_type> vjoin(dc, lgraph, rgraph);
  if ((algorithm & RUN_PAGERANK) && (algorithm & RUN_LDA)) {
    vjoin.prepare_injective_join(pagerank::left_emit_key, lda::right_emit_key);
    vjoinptr = &vjoin;
  }

 
  // Run pagerank ------------------------------------------------------------
  // Set up engine0 to compute transition probability
  graphlab::omni_engine<pagerank::compute_transit_prob> engine0(dc, lgraph, execution_type);
  if (algorithm  & RUN_PAGERANK)  {
    transit_prob_engine = &engine0;
    fn_compute_transit_prob();
  }

  // Set up engine1 to compute pagerank 
  graphlab::omni_engine<pagerank::compute_pagerank> engine1(dc, lgraph, execution_type); 
  if (algorithm & RUN_PAGERANK) {
    bool success = engine1.add_vertex_aggregator<pagerank::topk_aggregator>(
        "toppr", pagerank::topk_aggregator::map, pagerank::topk_aggregator::finalize) && engine1.aggregate_periodic("toppr", 5); 
    ASSERT_TRUE(success);
    pagerank_engine = &engine1;
    thgroup.launch(fn_run_pagerank);
  }
  
  // Run lda -----------------------------------------------------------------
  graphlab::omni_engine<lda::cgs_lda_vertex_program> engine2(dc, rgraph, execution_type);
  if (algorithm & RUN_LDA) {
    lda_engine = &engine2;
    lda::initialize_global();
    { 
      const bool success =
        engine2.add_vertex_aggregator<lda::topk_aggregator>
        ("topk", lda::topk_aggregator::map, lda::topk_aggregator::finalize) &&
        engine2.aggregate_periodic("topk", 5);
      ASSERT_TRUE(success);
    }

    { // Add the Global counts aggregator
      const bool success =
        engine2.add_vertex_aggregator<lda::factor_type>
        ("global_counts", 
         lda::global_counts_aggregator::map, 
         lda::global_counts_aggregator::finalize) &&
        engine2.aggregate_periodic("global_counts", 5);
      ASSERT_TRUE(success);
    }

    { // Add the likelihood aggregator
      const bool success =
        engine2.add_vertex_aggregator<lda::likelihood_aggregator>
        ("likelihood", 
         lda::likelihood_aggregator::map, 
         lda::likelihood_aggregator::finalize) &&
        engine2.aggregate_periodic("likelihood", 10);
      ASSERT_TRUE(success);
    }

    thgroup.launch(fn_run_lda);
  }

  // Run background vertex join
  if ((algorithm & RUN_PAGERANK) && (algorithm & RUN_LDA)) {
    thgroup.launch(fn_join_vertex);
  }

  thgroup.launch(boost::bind(launch_metric_server, algorithm));
  // --------------------------------------
  // ------------FINALIZE------------------
  // --------------------------------------
  graphlab::stop_metric_server_on_eof();
  thgroup.join();
  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
}
