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

#include <graphlab.hpp>

// ============================================================================>
//  Support code for loopy belief propagation
#include "image.hpp"



// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>



const size_t EDGE_POTENTIAL_ID = 0;


typedef graphlab::blob_types gl_types;

// GLOBAL listener handle. FIXME FIXME
gl_types::imonitor* listener(NULL);

// STRUCTS (Edge and Vertex data) =============================================>

/**
 * The data associated with each directed edge in the pairwise markov
 * random field
 */
class edge_data {
  graphlab::unary_factor _message;
  graphlab::unary_factor _old_message;
  size_t* _src_mapassg;
public:
  /** Create an empty edge data */
  edge_data() : 
    _message(NULL),
    _old_message(NULL),
    _src_mapassg(0) { }
  /** Interpret the pointer as edge data using the old_new message
      order specified by old_first */
  edge_data(void* data);
  /** Get the edge factor id */
  size_t& src_mapassg() { return *_src_mapassg; }
  /** Get the message */
  graphlab::unary_factor& message() { return _message; }
  graphlab::unary_factor& old_message() { return _old_message; }
  /** Allocate an edge factor blob */
  static graphlab::blob allocate_edge_data(size_t arity);
}; // End of edge data


/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
class vertex_data {
  graphlab::unary_factor _potential;
  graphlab::unary_factor _belief;
public:
  /** Create an empty vertex data */
  vertex_data() : _potential(NULL), _belief(NULL) { }
  /** Interpret the pointer as vertex data */
  vertex_data(void* data);
  /** Get the potential associated with the vertex */
  graphlab::unary_factor& potential() { return _potential; }
  /** Get the belief */
  graphlab::unary_factor& belief() { return _belief; }
  /** Allocate a vertex blob */
  static graphlab::blob allocate_vertex_data(size_t arity);
}; // End of vertex data



class edge_parameter : public graphlab::ireducible {
public:

  // this is a HACK but the binary_factor is stored in log form
  // and I need a counter. So I will store my counts as exp(counts)
  // zero the acc edge counts
  void zero_acc() {
    for (size_t i = 0; i < acc_edge_counts.arity1() ; ++i) {
      for (size_t j = 0; j < acc_edge_counts.arity2() ; ++j) {
        acc_edge_counts.logP(acc_edge_counts.var1(), i,
                                  acc_edge_counts.var2(), j) = 0;
      }
    }
  }
  
  edge_parameter(const graphlab::ireducible& other_) {
    this->operator=(other_);
  }
  
  edge_parameter(const edge_parameter& other_) {
    this->operator=(other_);
  }
  
  // constructor for creating the initial edge_parameter.
  edge_parameter(double ipfdamp, 
                 graphlab::blob initialedge,
                 graphlab::blob real_edge_c){
    // clone all the blobs we need
    edge_blob = initialedge.copy();
    real_edge_counts_blob = real_edge_c.copy();
    // make an edge accumulator
    acc_edge_counts_blob = real_edge_c.copy();
    
    // initialize the binary factors;
    edge_pot = graphlab::binary_factor(edge_blob.data());
    real_edge_counts = graphlab::binary_factor(real_edge_counts_blob.data());
    acc_edge_counts = graphlab::binary_factor(acc_edge_counts_blob.data());
    
    zero_acc();

    ipf_update_damping = ipfdamp;
  }
  
  ~edge_parameter() {
      acc_edge_counts_blob.clear();
  }
  
  ireducible& operator+=(const ireducible& other_) {
    bool nonzero = false;
    const edge_parameter& other = dynamic_cast<const edge_parameter&>(other_);
    for (size_t i = 0; i < acc_edge_counts.arity1() ; ++i) {
      for (size_t j = 0; j < acc_edge_counts.arity2() ; ++j) {
        double &myd = acc_edge_counts.logP(
                                      acc_edge_counts.var1(), i,
                                      acc_edge_counts.var2(), j);

        const double &otherd = other.acc_edge_counts.logP(
                                      other.acc_edge_counts.var1(), i,
                                      other.acc_edge_counts.var2(), j);
        myd += otherd;
        if (otherd) nonzero++;
      }
    }  
    return *this;
  }
  
  ireducible& operator=(const ireducible& other_) {
    const edge_parameter& other = dynamic_cast<const edge_parameter&>(other_);
    *this = other;
    return *this;
  }
  
  edge_parameter& operator=(const edge_parameter& other) {
    // copy the data
    edge_blob = other.edge_blob;
    real_edge_counts_blob = other.real_edge_counts_blob;
    acc_edge_counts_blob = other.acc_edge_counts_blob.copy();
    // create the binary factors
    acc_edge_counts = graphlab::binary_factor(acc_edge_counts_blob.data());
    edge_pot = graphlab::binary_factor(edge_blob.data());
    real_edge_counts = graphlab::binary_factor(real_edge_counts_blob.data());
    
    ipf_update_damping = other.ipf_update_damping;
    return *this;    
  }
  void add_count(size_t i, size_t j) {
    double &accij = acc_edge_counts.logP(
                                      acc_edge_counts.var1(), i,
                                      acc_edge_counts.var2(), j);
    accij++;
    double &accij2 = acc_edge_counts.logP(
                                  acc_edge_counts.var1(), j,
                                  acc_edge_counts.var2(), i);
    accij2++;
  }
  
  void apply(size_t reductions) { 
    // if counts are too small, don't do anything
    size_t count = 0;
    for (size_t i = 0; i < acc_edge_counts.arity1() ; ++i) {
      for (size_t j = 0; j < acc_edge_counts.arity2() ; ++j) {
          double &accij = acc_edge_counts.logP(
                                      acc_edge_counts.var1(), i,
                                      acc_edge_counts.var2(), j);
          count += accij;
      }
    }
    if (count < 10000) return;
    
    // normalize in the counts
    for (size_t i = 0; i < acc_edge_counts.arity1() ; ++i) {
      for (size_t j = 0; j < acc_edge_counts.arity2() ; ++j) {
          double &accij = acc_edge_counts.logP(
                                      acc_edge_counts.var1(), i,
                                      acc_edge_counts.var2(), j);
          accij /= count;
      }
    }

    std::cout << "estimated counts: " << std::endl;
    std::cout << acc_edge_counts;

    // perform the IPF update.
    for (size_t i = 0; i < acc_edge_counts.arity1() ; ++i) {
      for (size_t j = 0; j < acc_edge_counts.arity2() ; ++j) {
          // recall that the counts are stored in logP itself
          double &accij = acc_edge_counts.logP(
                                      acc_edge_counts.var1(), i,
                                      acc_edge_counts.var2(), j);

          double &realij = real_edge_counts.logP(
                                      acc_edge_counts.var1(), i,
                                      acc_edge_counts.var2(), j);

          double &logpotij = edge_pot.logP(
                                      acc_edge_counts.var1(), i,
                                      acc_edge_counts.var2(), j);
          double potij = exp(logpotij);
          // dampped updates
          potij = ipf_update_damping * potij + 
                  (1 - ipf_update_damping) * ((accij + 0.1) / (realij + 0.1)) * potij;
          
          if (potij <= 1E-5) potij = 1E-5;
          // write it back
          logpotij = log(potij);
      }
    }
    // renormalize the potential
    std::cout << "cur edge pot: " << std::endl;
    std::cout << edge_pot;

    edge_pot.normalize();
  };
  graphlab::blob edge_blob;
  graphlab::blob acc_edge_counts_blob;
  graphlab::blob real_edge_counts_blob;
  
  graphlab::binary_factor edge_pot;
  graphlab::binary_factor acc_edge_counts;
  graphlab::binary_factor real_edge_counts;
  double ipf_update_damping;
};

  
// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     graphlab::blob_graph& graph);

/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler,
               const graphlab::shared_data_manager* shared_data);



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

graphlab::blob get_edge_counts(const image &i, size_t arity);


// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using loopy belief propagation inside " << std::endl
            << "the graphlab framework." << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

  // Parse command line arguments --------------------------------------------->
  options opts; 
  bool success = parse_command_line(argc, argv, opts);
  if(!success) {
    return EXIT_FAILURE;
  }
  display_options(opts);
  

  // Create synthetic images -------------------------------------------------->
  // Creating image for denoising
  std::cout << "Creating a synethic image. " << std::endl;
  image img(opts.rows, opts.cols);
  img.paint_sunset(opts.num_rings);
  
  graphlab::blob realecounts = get_edge_counts(img, opts.num_rings);
  
  std::cout << "Saving image. " << std::endl;
  img.save(opts.orig_fn.c_str());
  std::cout << "Corrupting Image. " << std::endl;
  img.corrupt(opts.sigma);
  std::cout << "Saving corrupted image. " << std::endl;
  img.save(opts.noisy_fn.c_str());

  // Create the graph --------------------------------------------------------->
  graphlab::blob_graph graph;
  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(img, opts.num_rings, opts.sigma, graph);

  
  // Setup global shared variables -------------------------------------------->
  graphlab::shared_data_manager sdm(graph.num_vertices());
  // Initialize the edge agreement factor 
  std::cout << "Initializing shared edge agreement factor. " << std::endl;
  const size_t buffer_size(10000);
  size_t binfactsize = graphlab::binary_factor::get_sizeof(opts.num_rings,
                                                 opts.num_rings);
  assert(binfactsize <= buffer_size);
  char edge_pot_buffer[buffer_size];
  graphlab::blob edge_blob(binfactsize, edge_pot_buffer);
  graphlab::binary_factor edge_potential(edge_pot_buffer);
  edge_potential.var1() = 1;
  edge_potential.var2() = 2;
  edge_potential.arity1() = opts.num_rings;
  edge_potential.arity2() = opts.num_rings;
  if(opts.smoothing == "square") {
    edge_potential.set_as_agreement(opts.lambda);
  } else if (opts.smoothing == "laplace") {
    edge_potential.set_as_laplace(opts.lambda);
  } else {
    std::cout << "Invalid smoothing stype!" << std::endl;
    return EXIT_FAILURE;
  }
  
  edge_parameter ipfupdate(0.8, edge_blob, realecounts);
  ipfupdate.zero_acc();
  sdm.register_sync_data(0,  ipfupdate, 2000);
  sdm.register_shared_data(1, graphlab::blob(sizeof(double),
                                             &(opts.bound)));
  sdm.register_shared_data(2, graphlab::blob(sizeof(double), 
                                             &(opts.damping)));
  
  // Create the engine -------------------------------------------------------->
  gl_types::iengine* engine =
    graphlab::engine_factory::new_engine(opts.engine,
                                         opts.scheduler,
                                         opts.scope,
                                         graph,
                                         opts.ncpus);
  if(engine == NULL) {
    std::cout << "Unable to construct engine!" << std::endl;
    return EXIT_FAILURE;   
  }

  engine->set_shared_data(&sdm);
  // Set the splash size if necessary
  engine->set_option("splash_size", opts.splash_size);
  // Set the update function
  engine->set_option("update_function", bp_update);



  // Running the engine ------------------------------------------------------->
  std::cout << "Running the engine. " << std::endl;
  graphlab::timer timer; timer.start(); 

  // Aapo: made visualizer a dirty global variable FIXME
  if (opts.visualizer == "true") {  
    listener = new graphlab::visualizer_server_listener(&graph); 
  } else if(opts.visualizer == "console") {
    listener = new graphlab::simple_console_listener();
  }
  if(listener != NULL)
    engine->register_listener(listener);

  // engine->get_scheduler().set_sort_comparisor(graphlab::shuffler);
  // engine->get_scheduler().init();
  for (size_t i = 0; i < 5; ++i) {
  // Add an update task for each vertex with initial priority of 2.0
  // if (typeid(*engine) == typeid(graphlab::synchronous_engine)) {
  //   ((graphlab::synchronous_engine*)(engine))->set_update_function(bp_update);
  // } else {
    //for (graphlab::vertex_id_t vertex = 0; vertex < graph.num_vertices(); ++vertex){
    //  engine->get_scheduler().add_task(graphlab::update_task(vertex, bp_update), 100.0+float(rand())/RAND_MAX);
    //}
    // Change by Aapo: need to use add_task_to_all to utilize shuffler
    engine->add_task_to_all(bp_update, 100.0);
  
    engine->start();
  }

  double runtime = timer.current_time();
  size_t update_count = engine->last_update_count();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;  


  // Saving the output -------------------------------------------------------->
  std::cout << "Rendering the cleaned image. " << std::endl;
  if(opts.pred_type == "map") {
    for(size_t v = 0; v < graph.num_vertices(); ++v) {
      vertex_data vdata(graph.vertex_data(v).data());
      img.pixel(v) = vdata.belief().max_asg();    
    }
  } else if(opts.pred_type == "exp") {
    for(size_t v = 0; v < graph.num_vertices(); ++v) {
      vertex_data vdata(graph.vertex_data(v).data());
      img.pixel(v) = vdata.belief().expectation();
    }
  } else {
    std::cout << "Invalid prediction type! : " << opts.pred_type
              << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Saving cleaned image. " << std::endl;
  img.save(opts.pred_fn.c_str());

  std::cout << "Done!" << std::endl;
  
  if(listener != NULL) delete listener;
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
    ("bound",  boost_po::value<double>(&(opts.bound))->default_value(1E-5),
     "Residual termination bound")
    ("damping",  boost_po::value<double>(&(opts.damping))->default_value(0),
     "Residual termination bound")
    ("rings",  boost_po::value<size_t>(&(opts.num_rings))->default_value(5),
     "Number of rings in the noisy image")
    ("rows",  boost_po::value<size_t>(&(opts.rows))->default_value(200),
     "Number of rows in the noisy image")
    ("cols",  boost_po::value<size_t>(&(opts.cols))->default_value(200),
     "Number of columns in the noisy image")
    ("sigma",  boost_po::value<double>(&(opts.sigma))->default_value(1),
     "Standard deviation of noise.")
    ("lambda",  boost_po::value<double>(&(opts.lambda))->default_value(2),
     "Smoothness parameter (larger => smoother).")
    ("smoothing",
     boost_po::value<std::string>(&(opts.smoothing))->default_value("laplace"),
     "Options are {square, laplace}")
    ("engine",
     boost_po::value<std::string>(&(opts.engine))->default_value("threaded"),
     "Options are {threaded, sequential}")
    ("scope",
     boost_po::value<std::string>(&(opts.scope))->default_value("locked"),
     "Options are {unsync, locked}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("fifo"),
     "Options are {fifo, priority, sampling, splash}")
    ("splashsize",  boost_po::value<size_t>(&(opts.splash_size))->default_value(1000),
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
     "Use visualizer server {true, false, console}");
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
            << "visualizer:     " << opts.visualizer << std::endl;
} // end of display options




// Implementations ============================================================>
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler,
               const graphlab::shared_data_manager* shared_data) {
//  std::cout << scope.vertex();;
//  std::getchar();
  assert(shared_data != NULL);
  // Get the shared data
  double bound = shared_data->get_shared(1).as<double>();
  double damping = shared_data->get_shared(2).as<double>();

  const size_t buffer_size(sizeof(double)*100);

  // Create stack allocated cavity factor
  unsigned char cavity_buffer[buffer_size]; 
  graphlab::unary_factor cavity(cavity_buffer);

  // Create stack allocated temporary message
  unsigned char message_buffer[buffer_size]; 
  graphlab::unary_factor tmp_msg(message_buffer);

  // Grab the state from the scope -------------------------------------------->
  // Get the vertex data
  vertex_data v_data( scope.vertex_data().data() );
  // Get the in and out edges by reference
  const std::vector<graphlab::edge_id_t>& in_edges = 
    scope.in_edge_ids();
  const std::vector<graphlab::edge_id_t>& out_edges = 
    scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size());

  // Check that the buffer for the cavity is large enough
  assert(graphlab::unary_factor::get_sizeof(v_data.belief().arity()) 
         <= buffer_size);

  // Flip the old and new messages to improve safety when using the
  // unsynch scope
  foreach(graphlab::edge_id_t ineid, in_edges) {   
    // Get the in and out edge data
    edge_data in_edge( scope.edge_data(ineid).data() );
    // copy the message onto the out edge
    in_edge.old_message().set( in_edge.message() );
  }

  // Compute the belief ------------------------------------------------------->
  // Initialize the belief as the value of the factor 
  v_data.belief().set( v_data.potential() );
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the message
    edge_data e_data( const_cast<void*>(scope.edge_data(ineid).data()) );
    v_data.belief().times( e_data.old_message() );
  }
  v_data.belief().normalize();
  
  // Provide listner current value (by Aapo)
  if (listener != NULL ) 
    listener->app_set_vertex_value(scope.vertex(), 
                                   v_data.belief().expectation());

                                   
  // Get the edge factor

  edge_parameter edgeparam = shared_data->snapshot<edge_parameter>(EDGE_POTENTIAL_ID);
  edgeparam.zero_acc();
  graphlab::binary_factor &edge_factor = edgeparam.edge_pot;
  // build a post to the IPF update

  size_t mymapassg = v_data.belief().max_asg();
  // Compute outbound messages ------------------------------------------------>
  // Initialize cavity factor
  cavity.arity() = v_data.belief().arity();
  // Send outbound messages
  for(size_t i = 0; i < in_edges.size(); ++i) {
    // Get the edge ids
    graphlab::edge_id_t outeid = out_edges[i];
    graphlab::edge_id_t ineid = in_edges[i];
    // CLEVER HACK: Here we are expoiting the sorting of the edge ids
    // to do fast O(1) time edge reversal
    assert(scope.target(outeid) == scope.source(ineid));
    // Get the in and out edge data
    edge_data in_edge( const_cast<void*>(scope.edge_data(ineid).data()) );
    edge_data out_edge( scope.edge_data(outeid).data() );
    
    // Compute cavity
    cavity.set(v_data.belief()); // Copy the data from belief to cavity
    cavity.divide(in_edge.old_message()); // Make the cavity a cavity
    cavity.normalize();


    // convolve cavity with the edge factor storing the result in the
    // temporary message
    tmp_msg.arity() = out_edge.message().arity();
    tmp_msg.convolve(edge_factor, cavity, false);
    tmp_msg.normalize();

    // Damp the message
    tmp_msg.damp(out_edge.message(), damping);
    // Should not be necessary
    tmp_msg.normalize();
    
    // Compute message residual
    double residual = tmp_msg.residual(out_edge.old_message());
    out_edge.src_mapassg() = mymapassg;
    // Assign the out message
    out_edge.message().set(tmp_msg);
    
    // update the IPF counter
    edgeparam.add_count(mymapassg, in_edge.src_mapassg());
    // If The residuals is above threshold
    if(residual > bound) {
      gl_types::update_task task(scope.target(outeid), bp_update);      
      scheduler.add_task(task, residual);
    }
  }
  
  shared_data->post(EDGE_POTENTIAL_ID, scope.vertex(), edgeparam);
} // end of BP_update


// STRUCTS (Edge and Vertex data) =============================================>

edge_data::edge_data(void* data) {
  assert(data != NULL);
  _message = graphlab::unary_factor(data);
  _old_message = graphlab::unary_factor(_message.end());
  _src_mapassg = (size_t*)_old_message.end();
}



graphlab::blob edge_data::allocate_edge_data(size_t arity) {  
  graphlab::blob blob(sizeof(size_t) + 
                      2 * graphlab::unary_factor::get_sizeof(arity));
  graphlab::unary_factor tmp(blob.data());
  tmp.arity() = arity;
  graphlab::unary_factor tmp2(tmp.end());
  tmp2.arity() = arity;
  
  edge_data edge(blob.data());
  edge.src_mapassg() = 0;
  edge.message().arity() = arity;
  edge.message().zero();
  edge.message().normalize();
  edge.old_message().arity() = arity;
  edge.old_message().zero();
  edge.old_message().normalize();
  return blob;
} 

vertex_data::vertex_data(void* data) {
  assert(data != NULL);
  // construct the edge potentials
  _potential = graphlab::unary_factor(data);
  // Construct the belief starting at the end of the potential
  _belief = graphlab::unary_factor(_potential.end());
} 

graphlab::blob vertex_data::allocate_vertex_data(size_t arity) { 
  graphlab::blob blob( 2 * graphlab::unary_factor::get_sizeof(arity) );
  graphlab::unary_factor tmp(blob.data());
  tmp.arity() = arity;
  graphlab::unary_factor tmp2(tmp.end());
  tmp2.arity() = arity;
  vertex_data vdata(blob.data());
  vdata.potential().arity() = arity;
  vdata.potential().zero();
  vdata.potential().normalize();
  vdata.belief().arity() = arity;
  vdata.belief().zero();
  vdata.belief().normalize();
  return blob;  
}

void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     graphlab::blob_graph& graph) {
  // Construct a single blob for the vertex data
  graphlab::blob vblob = 
    vertex_data::allocate_vertex_data(num_rings);
  vertex_data vdata(vblob.data());
  // Add all the vertices
  double sigmaSq = sigma*sigma;
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      // Set the node potential
      double obs = img.pixel(i, j);
      for(size_t pred = 0; pred < num_rings; ++pred) {
        vdata.potential().logP(pred) = 
          -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
      }
      vdata.potential().normalize();
      // Store the actual data in the graph
      size_t vertid = graph.add_vertex(vblob);
      // Ensure that we are using a consistent numbering
      assert(vertid == img.vertid(i, j));
    } // end of for j in cols
  } // end of for i in rows
  // Add the edges
  graphlab::blob eblob = edge_data::allocate_edge_data(num_rings);
  edge_data edata (eblob.data());
  edata.src_mapassg() = 0;
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) 
        graph.add_edge(vertid, img.vertid(i-1, j), eblob);
      if(i+1 < img.rows())
        graph.add_edge(vertid, img.vertid(i+1, j), eblob);
      if(j-1 < img.cols())
        graph.add_edge(vertid, img.vertid(i, j-1), eblob);
      if(j+1 < img.cols())
        graph.add_edge(vertid, img.vertid(i, j+1), eblob);
    } // end of for j in cols
  } // end of for i in rows
  graph.finalize();
} // End of construct graph



graphlab::blob get_edge_counts(const image &img, size_t arity) {
  graphlab::blob b(graphlab::binary_factor::get_sizeof(arity,arity));
  graphlab::binary_factor edgecounts(b.data());
  edgecounts.var1() = 1;
  edgecounts.var2() = 2;
  edgecounts.arity1() = arity;
  edgecounts.arity2() = arity;
  //zero the counts
  for(size_t i = 0; i < arity; ++i) {
    for(size_t j = 0; j < arity; ++j) {
      edgecounts.logP(1, i, 2, j) = 0;
    }
  }
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      if(i-1 < img.rows())
        ++edgecounts.logP(1, img.pixel(i,j), 2, img.pixel(i-1,j));
      if(i+1 < img.rows())
        ++edgecounts.logP(1, img.pixel(i,j), 2, img.pixel(i+1,j));
      if(j-1 < img.cols())
        ++edgecounts.logP(1, img.pixel(i,j), 2, img.pixel(i,j-1));
      if(j+1 < img.cols())
        ++edgecounts.logP(1, img.pixel(i,j), 2, img.pixel(i,j+1));
    }
  }
  
  // normalize in the counts
  size_t c = 0;
  for(size_t i = 0; i < arity; ++i) {
    for(size_t j = 0; j < arity; ++j) {
      c += edgecounts.logP(1, i, 2, j);
    }
  }
  
  for(size_t i = 0; i < arity; ++i) {
    for(size_t j = 0; j < arity; ++j) {
      edgecounts.logP(1, i, 2, j) /= c;
    }
  }
  std::cout << "Real Edge Counts: " << std::endl;
  std::cout << edgecounts;
  return b;
}
