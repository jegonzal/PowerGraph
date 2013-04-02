#include "cgs_lda.hpp"

size_t NTOPICS = 20;
graphlab::distributed_control* dc_ptr;


/* Global thread func for running lda */
void fn_run_lda(lda::engine_type& lda_engine) {
  static graphlab::mutex m;
  if (m.try_lock()) {
    dc_ptr->cout() << "Running The Collapsed Gibbs Sampler" << std::endl;
    lda_engine.map_reduce_vertices<graphlab::empty>(lda::signal_only::docs);
    // Enable sampling
    lda::cgs_lda_vertex_program::DISABLE_SAMPLING = false;
    // Run the engine
    lda_engine.start();
    // Finalize the counts
    // lda::cgs_lda_vertex_program::DISABLE_SAMPLING = true;
    // lda_engine->signal_all();
    // lda_engine->start();
    m.unlock();
  }
}

void launch_metric_server() {
  graphlab::launch_metric_server();
  // Start the lda webserver 
  graphlab::add_metric_server_callback("ldaparam", lda::set_param_callback);
  graphlab::add_metric_server_callback("lockword", lda::lock_word_callback);
  graphlab::add_metric_server_callback("wordclouds", lda::word_cloud_callback);
  graphlab::add_metric_server_callback("addtopic", lda::add_topic_callback);
}


void load_ldagraph(lda::graph_type& graph,
                   const std::string& edge_dir,
                   const std::string& dictionary_dir) {
  lda::load_and_initialize_graph(*dc_ptr, graph, edge_dir);
  lda::load_dictionary(dictionary_dir);
}


int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  dc_ptr = &dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("./ldademo");
  std::string lda_edges;
  std::string lda_dictionary;

  bool USE_SYNC=false;

  // The dir of the link graph 
  clopts.attach_option("lda_edges", lda_edges, "The lda graph (edges). Required ");
  clopts.attach_option("lda_dictionary", lda_dictionary, "The lda word dictionary. Required ");
  clopts.attach_option("wordid_offset", lda::WORDID_OFFSET, "The starting id of the words in the lda input.");
  clopts.attach_option("ntopics", lda::NTOPICS, "Number of topics to use");
  clopts.attach_option("topklda", lda::TOPK, "Top k words in each topic for display");
  clopts.attach_option("force_lock", lda::FORCE_LOCK, "force locked words");
  clopts.attach_option("use_sync", USE_SYNC, "Use Synchronous LDA");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  if(!clopts.is_set("lda_edges")) {
    std::cout << "LDA edge file not provided" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  if(!clopts.is_set("lda_dictionary")) {
    std::cout << "LDA dictionary file not provided" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }

  lda::INITIAL_NTOPICS = NTOPICS;

  // Build the lda (right) graph
  // The global lda_graph points to rgraph
  lda::graph_type ldagraph(dc);
  load_ldagraph(ldagraph, lda_edges, lda_dictionary);

  // Run lda -----------------------------------------------------------------
  launch_metric_server();

  graphlab::graphlab_options opts;
  std::string engine_type = "sync";
  if (!USE_SYNC) {
    opts.get_engine_args().set_option("factorized",true);
    engine_type = "async";
  }
  //opts2.get_engine_args().set_option("handler_intercept",true);
  graphlab::omni_engine<lda::cgs_lda_vertex_program> engine(dc, ldagraph, engine_type, opts);
  lda::initialize_global();
  { 
    const bool success =
        engine.add_vertex_aggregator<lda::topk_aggregator>
        ("topk", lda::topk_aggregator::map, lda::topk_aggregator::finalize) &&
        engine.aggregate_periodic("topk", 1);
    ASSERT_TRUE(success);
  }

  { // Add the Global counts aggregator
    const bool success =
        engine.add_vertex_aggregator<lda::factor_type>
        ("global_counts", 
         lda::global_counts_aggregator::map, 
         lda::global_counts_aggregator::finalize) &&
        engine.aggregate_periodic("global_counts", 1);
    ASSERT_TRUE(success);
  }
  fn_run_lda(engine);

// ------------FINALIZE------------------
// --------------------------------------
graphlab::stop_metric_server_on_eof();
graphlab::mpi_tools::finalize();

return EXIT_SUCCESS;
}
