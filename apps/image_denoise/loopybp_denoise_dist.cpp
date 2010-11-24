/**
 * This file contains an example of graphlab used for discrete loopy
 * belief propagation in a pairwise markov random field to denoise a
 * synthetic noisy image.
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>

#include <boost/program_options.hpp>
#include <graphlab/distributed/graph/cloned_graph.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
#include <graphlab.hpp>

#include "image.hpp"


// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>


enum constants {EDGE_FACTOR_ID, BOUND_ID, DAMPING_ID};



// STRUCTS (Edge and Vertex data) =============================================>

/**
 * The data associated with each directed edge in the pairwise markov
 * random field
 */
struct edge_data {
  graphlab::unary_factor message;
  graphlab::unary_factor old_message;
  void save(graphlab::oarchive &oarc) const{
    oarc << message;
    oarc << old_message;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> message;
    iarc >> old_message;
  }
}; // End of edge data


/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data {
  graphlab::unary_factor potential;
  graphlab::unary_factor belief;
  void save(graphlab::oarchive &oarc) const{
    oarc << potential;
    oarc << belief;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> potential;
    iarc >> belief;
  }
}; // End of vertex data


typedef graphlab::cloned_graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;


// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl_types::graph& graph);

/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data);
               

// Command Line Parsing =======================================================>

struct options {
  size_t ncpus;
  double bound;
  double damping;
  size_t num_rings;
  size_t rows;
  size_t cols;
  double sigma;
  double lambda;
  size_t splash_size;
  std::string smoothing;
  std::string engine;
  std::string scope;
  std::string scheduler;
  std::string orig_fn;
  std::string noisy_fn;
  std::string pred_fn;
  std::string pred_type;
  std::string visualizer;
  std::string partmethod;
  size_t clustersize;
};


/**
 * Parse the command line arguments.  Returns false if there was a
 * problem in parsing command line arguments
 */    
bool parse_command_line(int argc, char** argv, options& opts);

/**
 * Display the program options
 */
void display_options(options& opts);

// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using loopy belief propagation inside " << std::endl
            << "the graphlab framework." << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::distributed_control dc(&argc, &argv);
  dc.init_message_processing(4);
  dc.barrier();
  // Create the graph --------------------------------------------------------->
  gl_types::graph graph;
  options opts; 
  bool success = parse_command_line(argc, argv, opts);

  if (dc.procid() == 0) {
    // Parse command line arguments --------------------------------------------->
    if(!success) {
      return EXIT_FAILURE;
    }
    display_options(opts);
  

    // Create synthetic images -------------------------------------------------->
    // Creating image for denoising
    std::cout << "Creating a synethic image. " << std::endl;
    image img(opts.rows, opts.cols);
    img.paint_sunset(opts.num_rings);
    std::cout << "Saving image. " << std::endl;
    img.save(opts.orig_fn.c_str());
    std::cout << "Corrupting Image. " << std::endl;
    img.corrupt(opts.sigma);
    std::cout << "Saving corrupted image. " << std::endl;
    img.save(opts.noisy_fn.c_str());

    std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
    construct_graph(img, opts.num_rings, opts.sigma, graph);
  }
  dc.barrier();
  graph.distributed_partition(dc, graphlab::partition_method::PARTITION_METIS,2);
  graph.distribute(dc);
  dc.barrier();
  
  // Setup global shared variables -------------------------------------------->
  graphlab::distributed_shared_data<gl_types::graph> sdm(dc);
  dc.barrier();
  if (dc.procid() == 0) {
    // Initialize the edge agreement factor 
    std::cout << "Initializing shared edge agreement factor. " << std::endl;

    // dummy variables 0 and 1 and num_rings by num_rings
    graphlab::binary_factor edge_potential(0, opts.num_rings,
                                           0, opts.num_rings);
    // Set the smoothing type
    if(opts.smoothing == "square") {
      edge_potential.set_as_agreement(opts.lambda);
    } else if (opts.smoothing == "laplace") {
      edge_potential.set_as_laplace(opts.lambda);
    } else {
      std::cout << "Invalid smoothing stype!" << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << edge_potential << std::endl;
    
    sdm.set_constant(EDGE_FACTOR_ID, graphlab::any(edge_potential));
    sdm.set_constant(BOUND_ID, graphlab::any(opts.bound));
    sdm.set_constant(DAMPING_ID, graphlab::any(opts.damping));
  }
  dc.barrier();
  
  // Create the engine -------------------------------------------------------->
  gl_types::iengine* engine =
    new graphlab::distributed_engine<gl_types::graph,
    graphlab::distributed_scheduler_wrapper<gl_types::graph, 
    graphlab::fifo_scheduler<gl_types::graph> > >(dc, graph, 2);

  dc.barrier();

  if(engine == NULL) {
    std::cout << "Unable to construct engine!" << std::endl;
    return EXIT_FAILURE;   
  }
  
  engine->set_shared_data_manager(&sdm);
  engine->set_default_scope(graphlab::scope_range::EDGE_CONSISTENCY);
  graphlab::timer timer; 
  if (dc.procid() == 0) {
    std::cout << "Running the engine. " << std::endl;
  
    // Add the bp update to all vertices
    std::vector<graphlab::vertex_id_t> vec;
    for (graphlab::vertex_id_t i = 0; i < graph.num_vertices(); ++i) {
      vec.push_back(i);
    }
    std::random_shuffle(vec.begin(), vec.end());
    engine->add_tasks(vec, bp_update, 100.0);
    //engine->get_scheduler().add_task_to_all(bp_update, 100.0);
  }
  dc.barrier();
  timer.start();
  // Starte the engine
  engine->start();
  double runtime = timer.current_time();
  dc.barrier();
  for (size_t i = 0;i < graph.my_vertices().size(); ++i) {
    graph.update_vertex(graph.my_vertices()[i]);
  }
  dc.barrier();

  if (dc.procid() == 0) {

    size_t update_count = engine->last_update_count();
    std::cout << "Finished Running engine in " << runtime 
              << " seconds." << std::endl
              << "Total updates: " << update_count << std::endl
              << "Efficiency: " << (double(update_count) / runtime)
              << " updates per second "
              << std::endl;  
  
    image img(opts.rows, opts.cols);
    // Saving the output -------------------------------------------------------->
    std::cout << "Rendering the cleaned image. " << std::endl;
    if(opts.pred_type == "map") {
      for(size_t v = 0; v < graph.num_vertices(); ++v) {
        const vertex_data& vdata = graph.vertex_data(v);
        img.pixel(v) = vdata.belief.max_asg();    
      }
    } else if(opts.pred_type == "exp") {
      for(size_t v = 0; v < graph.num_vertices(); ++v) {
        const vertex_data& vdata = graph.vertex_data(v);
        img.pixel(v) = vdata.belief.expectation();
      }
    } else {
      std::cout << "Invalid prediction type! : " << opts.pred_type
                << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Saving cleaned image. " << std::endl;

    std::stringstream ss; ss << dc.procid(); 
    opts.pred_fn = opts.pred_fn + ss.str() + ".pgm";
    img.save(opts.pred_fn.c_str());
  }
  dc.barrier();
  std::cout << "Done!" << std::endl;
  
  if(engine != NULL) delete engine;

  return EXIT_SUCCESS;
} // End of main

bool parse_command_line(int argc, char** argv, options& opts) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Denoise a randomly generated image using Gibbs Sampling.");
  // Set the program options
  desc.add_options()
    ("help",   "produce this help message")
    ("ncpus",  boost_po::value<size_t>(&(opts.ncpus))->default_value(2),
     "Number of cpus to use.")
    ("bound",  boost_po::value<double>(&(opts.bound))->default_value(1E-15),
     "Residual termination bound")
    ("damping",  boost_po::value<double>(&(opts.damping))->default_value(0),
     "Residual termination bound")
    ("rings",  boost_po::value<size_t>(&(opts.num_rings))->default_value(10),
     "Number of rings in the noisy image")
    ("rows",  boost_po::value<size_t>(&(opts.rows))->default_value(300),
     "Number of rows in the noisy image")
    ("cols",  boost_po::value<size_t>(&(opts.cols))->default_value(300),
     "Number of columns in the noisy image")
    ("sigma",  boost_po::value<double>(&(opts.sigma))->default_value(2),
     "Standard deviation of noise.")
    ("lambda",  boost_po::value<double>(&(opts.lambda))->default_value(10),
     "Smoothness parameter (larger => smoother).")
    ("smoothing",
     boost_po::value<std::string>(&(opts.smoothing))->default_value("laplace"),
     "Options are {square, laplace}")
    ("engine",
     boost_po::value<std::string>(&(opts.engine))->default_value("threaded"),
     "Options are {threaded, sequential}")
    ("scope",
     boost_po::value<std::string>(&(opts.scope))->default_value("edge"),
     "Options are {vertex, edge, full}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("fifo"),
     "Options are {fifo, priority, sampling, splash}")
    ("splashsize",  boost_po::value<size_t>(&(opts.splash_size))->default_value(10000),
     "The desired splash size.")
    ("orig",
     boost_po::value<std::string>(&(opts.orig_fn))->default_value("source_img.pgm"),
     "Original image file name.")
    ("noisy",
     boost_po::value<std::string>(&(opts.noisy_fn))->default_value("noisy_img.pgm"),
     "Noisy image file name.")
    ("pred",
     boost_po::value<std::string>(&(opts.pred_fn))->default_value("pred_img.pgm"),
     "Predicted image file name.")
    ("pred_type",
     boost_po::value<std::string>(&(opts.pred_type))->default_value("map"),
     "Predicted image type {map, exp}")
    ("visualizer",
     boost_po::value<std::string>(&(opts.visualizer))->default_value("false"),
     "Use visualizer server {true, false, console}")
    ("clustersize",
     boost_po::value<size_t>(&(opts.clustersize))->default_value(100),
     "Number of vertices per cluster (only for clustered_priority_scheduler) ")
    ("partmethod",
     boost_po::value<std::string>(&(opts.partmethod))->default_value("metis"),
     "metis / random. Only for clustered_priority_scheduler");

  // Parse the arguments
  boost_po::variables_map vm;
  boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
  boost_po::notify(vm);
  if(vm.count("help")) {
    std::cout << desc << std::endl;
    return false;
  }
  return true;
} // end of parse command line arguments



void display_options(options& opts) {
  std::cout << "ncpus:          " << opts.ncpus << std::endl
            << "bound:          " << opts.bound << std::endl
            << "damping:        " << opts.damping << std::endl
            << "num_rings:      " << opts.num_rings << std::endl
            << "rows:           " << opts.rows << std::endl
            << "cols:           " << opts.cols << std::endl
            << "sigma:          " << opts.sigma << std::endl
            << "lambda:         " << opts.lambda << std::endl
            << "smoothing:      " << opts.smoothing << std::endl
            << "engine:         " << opts.engine << std::endl
            << "scope:          " << opts.scope << std::endl
            << "scheduler:      " << opts.scheduler << std::endl
            << "splash_size:    " << opts.splash_size << std::endl
            << "orig_fn:        " << opts.orig_fn << std::endl
            << "noisy_fn:       " << opts.noisy_fn << std::endl
            << "pred_fn:        " << opts.pred_fn << std::endl
            << "pred_type:      " << opts.pred_type << std::endl
            << "visualizer:     " << opts.visualizer << std::endl
            << "clustersize:    " << opts.clustersize<< std::endl
            << "partmethod:     " << opts.partmethod<< std::endl;
} // end of display options



size_t ctr = 0;
// Implementations
// ============================================================>
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data) {
   ++ctr ;
  //  std::cout << scope.vertex();;
  //  std::getchar();
  assert(shared_data != NULL);

  // Get the shared data
  double bound = shared_data->get_constant(BOUND_ID).as<double>();
  double damping = shared_data->get_constant(DAMPING_ID).as<double>();
  if (ctr == 1) {
    std::cout << "damping: " << damping << "\n";
    std::cout << "bound: " << bound << "\n";
  }
  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  
  // Get the in and out edges by reference
  graphlab::edge_list in_edges = 
    scope.in_edge_ids();
  graphlab::edge_list out_edges = 
    scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size()); // Sanity check

  // Flip the old and new messages to improve safety when using the
  // unsynch scope
  foreach(graphlab::edge_id_t ineid, in_edges) {   
    // Get the in and out edge data
    edge_data& in_edge = scope.edge_data(ineid);
    // Since we are about to receive the current message make it the
    // old message
    in_edge.old_message = in_edge.message;
  }

  // Compute the belief
  // ---------------------------------------------------------------->
  // Initialize the belief as the value of the factor
  v_data.belief = v_data.potential;
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the message
    const edge_data& e_data = scope.const_edge_data(ineid);
    // Notice we now use the old message since neighboring vertices
    // could be changing the new messages
    v_data.belief.times( e_data.old_message );
  }
  v_data.belief.normalize(); // finally normalize the belief
  
  // Compute outbound messages
  // ---------------------------------------------------------------->

  const graphlab::binary_factor edge_factor =
    shared_data->get_constant(EDGE_FACTOR_ID).as<graphlab::binary_factor>();

  
  // Send outbound messages
  graphlab::unary_factor cavity, tmp_msg;
  for(size_t i = 0; i < in_edges.size(); ++i) {
    // Get the edge ids
    graphlab::edge_id_t outeid = out_edges[i];
    graphlab::edge_id_t ineid = in_edges[i];
    // CLEVER HACK: Here we are expoiting the sorting of the edge ids
    // to do fast O(1) time edge reversal
    assert(scope.target(outeid) == scope.source(ineid));
    // Get the in and out edge data
    const edge_data& in_edge = scope.const_edge_data(ineid);
    edge_data& out_edge = scope.edge_data(outeid);
    
    // Compute cavity
    cavity = v_data.belief;
    cavity.divide(in_edge.old_message); // Make the cavity a cavity
    cavity.normalize();


    // convolve cavity with the edge factor storing the result in the
    // temporary message
    tmp_msg.resize(out_edge.message.arity());
    tmp_msg.var() = out_edge.message.var();
    tmp_msg.convolve(edge_factor, cavity);
    tmp_msg.normalize();

    // Damp the message
    tmp_msg.damp(out_edge.message, damping);
    
    // Compute message residual
    double residual = tmp_msg.residual(out_edge.old_message);
    
    // Assign the out message
    out_edge.message = tmp_msg;
    
    if(residual > bound) {
      gl_types::update_task task(scope.target(outeid), bp_update);      
      scheduler.add_task(task, residual);
    }    
  }
} // end of BP_update


void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl_types::graph& graph) {
  // Construct a single blob for the vertex data
  vertex_data vdata;
  vdata.potential.resize(num_rings);
  vdata.belief.resize(num_rings);
  vdata.belief.uniform();
  vdata.potential.uniform();
  vdata.belief.normalize();
  vdata.potential.normalize();
  // Add all the vertices
  double sigmaSq = sigma*sigma;
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      // initialize the potential and belief
      uint32_t pixel_id = img.vertid(i, j);
      vdata.potential.var() = vdata.belief.var() = pixel_id;
      // Set the node potential
      double obs = img.pixel(i, j);
      for(size_t pred = 0; pred < num_rings; ++pred) {
        vdata.potential.logP(pred) = 
          -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
      }
      vdata.potential.normalize();
      // Store the actual data in the graph
      size_t vertid = graph.add_vertex(vdata);
      // Ensure that we are using a consistent numbering
      assert(vertid == img.vertid(i, j));
    } // end of for j in cols
  } // end of for i in rows

  // Add the edges
  edge_data edata;
  edata.message.resize(num_rings);
  edata.message.uniform();
  edata.message.normalize();
  edata.old_message = edata.message;
  
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        edata.message.var() = img.vertid(i-1, j);
        edata.old_message.var() = edata.message.var();
        graph.add_edge(vertid, img.vertid(i-1, j), edata);
      }
      if(i+1 < img.rows()) {
        edata.message.var() = img.vertid(i+1, j);
        edata.old_message.var() = edata.message.var();
        graph.add_edge(vertid, img.vertid(i+1, j), edata);
      }
      if(j-1 < img.cols()) {
        edata.message.var() = img.vertid(i, j-1);
        edata.old_message.var() = edata.message.var();
        graph.add_edge(vertid, img.vertid(i, j-1), edata);
      } if(j+1 < img.cols()) {
        edata.message.var() = img.vertid(i, j+1);
        edata.old_message.var() = edata.message.var();
        graph.add_edge(vertid, img.vertid(i, j+1), edata);
      }
    } // end of for j in cols
  } // end of for i in rows
  graph.finalize();  
} // End of construct graph


