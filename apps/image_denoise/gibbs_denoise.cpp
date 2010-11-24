/**
 * This file contains an example of graphlab used for gibbs sampling
 * in a pairwise markov random field to denoise a synthetic noisy
 * image.
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
#include <boost/bind.hpp>

#include <graphlab.hpp>

// ============================================================================>
//  Support code for loopy belief propagation
#include "image.hpp"

// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

// Determine the shared ata location of the shared edge factor and the
// number of samples
enum constants {EDGE_FACTOR_ID, NUM_SAMPLES_ID};

// STRUCTS (Edge and Vertex data) =============================================>
/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data {
  size_t asg;
  size_t updates;
  unsigned char color;
  graphlab::unary_factor potential;
  //! store the Roa-Blackwell conditional belief estimate
  graphlab::unary_factor belief;
  std::vector<double> counts;
  vertex_data() :
    asg(0), updates(0), color(0) { }
}; // End of vertex data

/**
 * The data associated with each directed edge in the pairwise markov
 * random field
 */
struct edge_data { }; 


typedef graphlab::graph< vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl;


// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl::graph& graph);

/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
template<bool UseCallback>
void gibbs_update(gl::iscope& scope, 
                  gl::icallback& scheduler,
                  gl::ishared_data* shared_data);

/**
 *  The colore schedule is used by the set scheduler to describe an
 *  execution that updates all vertices of the same color (given by a
 *  checkerboard) pattern) in parallel.
 */
void color_schedule(gl::set_scheduler &sched);
// This ugly global variable is needed for the color_schedule
size_t nsamples;



// Command Line Parsing =======================================================>
struct options {
  size_t ncpus;
  size_t samples;
  size_t num_rings;
  size_t rows;
  size_t cols;
  double sigma;
  double lambda;
  std::string smoothing;
  std::string engine;
  std::string scope;
  std::string scheduler;
  std::string orig_fn;
  std::string noisy_fn;
  std::string rb_pred_fn;
  std::string counts_pred_fn;
  std::string pred_type;
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
  std::cout << "This program uses gibbs sampling to denoise a synthetic image."
            << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

  // Parse command line arguments --------------------------------------------->
  options opts; 
  bool success = parse_command_line(argc, argv, opts);
  if(!success)  return EXIT_FAILURE;
  display_options(opts);

  // Create synthetic images -------------------------------------------------->
  // Creating image for denoising
  std::cout << "Creating a synethic image." << std::endl;
  image img(opts.rows, opts.cols);
  img.paint_sunset(opts.num_rings);
  std::cout << "Saving image. " << std::endl;
  img.save(opts.orig_fn.c_str());
  std::cout << "Corrupting Image. " << std::endl;
  img.corrupt(opts.sigma);
  std::cout << "Saving corrupted image. " << std::endl;
  img.save(opts.noisy_fn.c_str());
  img.save_vec("corrupted.tsv");
  
  
  // Create the graph --------------------------------------------------------->
  gl::graph graph;
  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(img, opts.num_rings, opts.sigma, graph);
  image coloring_img(opts.rows, opts.cols);
  for(size_t i = 0; i < graph.num_vertices(); ++i)
    coloring_img.pixel(i) = graph.vertex_data(i).color;
  coloring_img.save("coloring.pgm");
  
 


  

  // Setup global shared variables -------------------------------------------->
  gl::thread_shared_data sdm;
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
  sdm.set_constant(NUM_SAMPLES_ID, graphlab::any(opts.samples));


    
  
  // Create the engine -------------------------------------------------------->
  gl::iengine* engine = NULL;
  if(opts.scheduler == "set") {
    nsamples = opts.samples;
    
    // Here we use the set scheduler which has some additional special setup
    engine =
      graphlab::engine_factory::new_engine(opts.engine,
                                           "set",
                                           "null",
                                           graph,
                                           opts.ncpus);
    if(engine == NULL) {
      std::cout << "Unable to construct engine!" << std::endl;
      return EXIT_FAILURE;   
    }
    // Set the shared data manager for the engine
    engine->set_shared_data_manager(&sdm);
    // Attach the correct schedule
    engine->get_scheduler().set_option(gl::scheduler_options::SCHEDULING_FUNCTION,
                                       (void*)color_schedule);
    
  } else {
    // Here we use the other schedulers which use a more default setup
    engine =
      graphlab::engine_factory::new_engine(opts.engine,
                                           opts.scheduler,
                                           opts.scope,
                                           graph,
                                           opts.ncpus);
    if(engine == NULL) {
      std::cout << "Unable to construct engine!" << std::endl;
      return EXIT_FAILURE;   
    }
    // Set the shared data manager for the engine
    engine->set_shared_data_manager(&sdm);

    // Create a shuffled set of vertices
    std::vector<gl::vertex_id_t> vertex_ids(opts.rows * opts.cols);
    for(size_t i = 0; i < vertex_ids.size(); ++i) vertex_ids[i] = i;
    std::random_shuffle(vertex_ids.begin(), vertex_ids.end());

    // Add the tasks in shuffled order (helps with fifo scheduler)
    const bool use_callback = true;
    double residual = 1.0;
    engine->add_tasks(vertex_ids,
                                      gibbs_update<use_callback>,
                                      residual);
  }

  
  // Running the engine ------------------------------------------------------->
  std::cout << "Running the engine. " << std::endl;
  graphlab::timer timer; timer.start(); 
  engine->start();
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
  image rb_img(opts.rows, opts.cols);
  image counts_img(opts.rows, opts.cols);
  if(opts.pred_type == "map") {
    for(size_t v = 0; v < graph.num_vertices(); ++v) {
      const vertex_data& vdata = graph.vertex_data(v);
      rb_img.pixel(v) = vdata.belief.max_asg();
      counts_img.pixel(v) =
        std::max_element(vdata.counts.begin(), vdata.counts.end()) -
        vdata.counts.begin();
        
    }
  } else if(opts.pred_type == "exp") {
    for(size_t v = 0; v < graph.num_vertices(); ++v) {
      const vertex_data& vdata = graph.vertex_data(v);    
      rb_img.pixel(v) = vdata.belief.expectation();
      double expectation = 0;
      for(size_t i = 0; i < vdata.counts.size(); ++i)
        expectation += i * vdata.counts[i];
      counts_img.pixel(v) = expectation;
    }
  } else {
    std::cout << "Invalid prediction type! : " << opts.pred_type
              << std::endl;
    return EXIT_FAILURE;
  }

  
  std::cout << "Saving cleaned image. " << std::endl;
  rb_img.save(opts.rb_pred_fn.c_str());
  counts_img.save(opts.counts_pred_fn.c_str());
  
  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
} // End of main



// Implementations ============================================================>
template<bool UseCallback>
void gibbs_update(gl::iscope& scope, 
                  gl::icallback& scheduler,
                  gl::ishared_data* shared_data) {

  assert(shared_data != NULL);
  // Get the shared data
  size_t samples = shared_data->get_constant(NUM_SAMPLES_ID).as<size_t>();
  

  // Grab the state from the scope -------------------------------------------->
  // Get the vertex data
  vertex_data& vdata = scope.vertex_data();

  // If this vertex has been updated sufficiently often then we don't
  // update again
  if(vdata.updates >= samples) return;  

  // Construct the conditional------------------------------------------------->
  // Initialize the conditional as the node potential
  graphlab::unary_factor conditional(vdata.potential);

  const graphlab::binary_factor& edge_factor =
    shared_data->get_constant(EDGE_FACTOR_ID).as<graphlab::binary_factor>();

  
  // Condition on all neighbor assignments
  foreach(graphlab::edge_id_t edgeid, scope.in_edge_ids()) {
    graphlab::vertex_id_t source = scope.source(edgeid);
    // Get the neighboring assignment
    const vertex_data& neighbor = scope.neighbor_vertex_data(source);
    // Get the edge factor
    conditional.condition(edge_factor, neighbor.asg);
    conditional.normalize();
  }
  
  // Generate a sample -------------------------------------------------------->
  size_t sample = conditional.sample();
  vdata.asg = sample;  // Record the sample
  vdata.updates++;     // Increment the updates
  vdata.counts[vdata.asg]++; // increment the counts
  
  // Rao-Blackwell belief estimate 
  vdata.belief.plus(conditional);
  
  // Finalize belief on last update ------------------------------------------->
  if(vdata.updates == samples ) { // If this is the last update
    vdata.belief.normalize();
    double Z = 0;
    for(size_t i = 0; i < vdata.counts.size(); ++i)
      Z += vdata.counts[i];
    for(size_t i = 0; i < vdata.counts.size(); ++i)
      vdata.counts[i] /= Z;
  } else if(UseCallback) {
    // Otherwise reschedule self
    double residual = 1.0;
    gl::update_task task(scope.vertex(), gibbs_update<UseCallback> );      
    scheduler.add_task(task, residual);
  }
} // end of BP_update



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
    ("samples",  boost_po::value<size_t>(&(opts.samples))->default_value(100),
     "Number of samples to generate")
    ("rings",  boost_po::value<size_t>(&(opts.num_rings))->default_value(5),
     "Number of rings in the noisy image")
    ("rows",  boost_po::value<size_t>(&(opts.rows))->default_value(200),
     "Number of rows in the noisy image")
    ("cols",  boost_po::value<size_t>(&(opts.cols))->default_value(200),
     "Number of columns in the noisy image")
    ("sigma",  boost_po::value<double>(&(opts.sigma))->default_value(1.2),
     "Standard deviation of noise.")
    ("lambda",  boost_po::value<double>(&(opts.lambda))->default_value(3),
     "Smoothness parameter (larger => smoother).")
    ("smoothing",
     boost_po::value<std::string>(&(opts.smoothing))->default_value("square"),
     "Options are {square, laplace}")
    ("engine",
     boost_po::value<std::string>(&(opts.engine))->default_value("threaded"),
     "Options are {threaded, sequential}")
    ("scope",
     boost_po::value<std::string>(&(opts.scope))->default_value("edge"),
     "Options are {vertex, edge, full}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("fifo"),
     "Options are {fifo, priority, sampling}")
    ("orig",
     boost_po::value<std::string>(&(opts.orig_fn))->default_value("source_img.pgm"),
     "Original image file name.")
    ("noisy",
     boost_po::value<std::string>(&(opts.noisy_fn))->default_value("noisy_img.pgm"),
     "Noisy image file name.")
    ("rb_pred",
     boost_po::value<std::string>(&(opts.rb_pred_fn))->
     default_value("rb_pred_img.pgm"),
     "Predicted image file name for the Rao-Blackwell estimator.")
    ("counts_pred",
     boost_po::value<std::string>(&(opts.counts_pred_fn))->
     default_value("counts_pred_img.pgm"),
     "Predicted image file name for the raw counts estimator.")
    ("pred_type",
     boost_po::value<std::string>(&(opts.pred_type))->default_value("map"),
     "Predicted image type {map, exp}");
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
            << "samples:        " << opts.samples << std::endl
            << "num_rings:      " << opts.num_rings << std::endl
            << "rows:           " << opts.rows << std::endl
            << "cols:           " << opts.cols << std::endl
            << "sigma:          " << opts.sigma << std::endl
            << "lambda:         " << opts.lambda << std::endl
            << "smoothing:      " << opts.smoothing << std::endl
            << "engine:         " << opts.engine << std::endl
            << "scope:          " << opts.scope << std::endl
            << "scheduler:      " << opts.scheduler << std::endl
            << "orig_fn:        " << opts.orig_fn << std::endl
            << "noisy_fn:       " << opts.noisy_fn << std::endl
            << "rb_pred_fn:     " << opts.rb_pred_fn << std::endl
            << "counts_pred_fn: " << opts.counts_pred_fn << std::endl
            << "pred_type:      " << opts.pred_type << std::endl;
}



// STRUCTS (Edge and Vertex data) =============================================>

void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl::graph& graph) {

  // Initialize the vertex data to somethine sensible
  vertex_data vdata;

  vdata.belief.resize(num_rings);
  vdata.belief.uniform(-std::numeric_limits<double>::max());
  vdata.potential.resize(num_rings);
  vdata.potential.uniform();
  vdata.potential.normalize();
  vdata.counts.resize(num_rings, 0);
  for(size_t i = 0; i < num_rings; ++i)
    vdata.counts[i] = 0;

        
  // Add all the vertices
  double sigmaSq = sigma*sigma;
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      // initialize the potential and belief
      uint32_t pixel_id = img.vertid(i, j);
      vdata.potential.var() = vdata.belief.var() = pixel_id;
      // Determine the color (in checkerboard form)
      vdata.color = (i%2 == 0) ^ (j%2 == 0)? 1 : 0;
      // Set the node potential
      double obs = img.pixel(i, j);
      for(size_t pred = 0; pred < num_rings; ++pred) {
        vdata.potential.logP(pred) = 
          -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
      }
      vdata.potential.normalize();
      // Set the initial assignment
      vdata.asg = vdata.potential.sample();
      // Store the actual data in the graph
      size_t vert_id = graph.add_vertex(vdata);
      // Ensure that we are using a consistent numbering
      assert(vert_id == pixel_id);
    } // end of for j in cols
  } // end of for i in rows
  // Construct an edge blob

  // Add all the edges
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) 
        graph.add_edge(vertid, img.vertid(i-1, j));
      if(i+1 < img.rows())
        graph.add_edge(vertid, img.vertid(i+1, j));
      if(j-1 < img.cols())
        graph.add_edge(vertid, img.vertid(i, j-1));
      if(j+1 < img.cols())
        graph.add_edge(vertid, img.vertid(i, j+1));
    } // end of for j in cols
  } // end of for i in rows
} // End of construct graph


/** Color selection function */
template<size_t color>
bool select_color(graphlab::vertex_id_t v, 
                  const vertex_data& vdata) {
  return vdata.color == color;
}

/**
 * The scheduling function used by the set scheduler to execute in
 * colored order
 */
void color_schedule(gl::set_scheduler &sched) {
  // There is a two coloring of this model
  // All sets must be created before scheduling calls
  std::vector<gl::ivertex_set*> colorsets(2);
  // Colors do not change during execution
  bool fixed_selection = true;
  gl::selector_function select_color_0( select_color<0> );
  colorsets[0] = &sched.attach(gl::rvset(select_color_0, fixed_selection),
                               sched.root_set());
  gl::selector_function select_color_1( select_color<1> );
  colorsets[1] = &sched.attach(gl::rvset(select_color_1, fixed_selection),
                               sched.root_set());
  // Actually construct the color sets
  sched.init();

  // Build the execution plan
  gl::execution_plan eplan;
  for(size_t i = 0; i < colorsets.size(); ++i) {
    std::cout << "Set " << i << " has "
              << colorsets[i]->size() << " vertices."
              << std::endl;
    const bool use_callback = false;
    eplan.execute(*colorsets[i], gibbs_update<use_callback>);
  }

  // Compile the execution plan
  graphlab::timer timer; timer.start();
  eplan.generate_plan(sched.get_graph(), sched.num_cpus());
  std::cout << "Execution plan compiled in "
            << timer.current_time() << " seconds" << std::endl;

  // Run the plan nsamples times
  for(size_t iter = 0; iter < nsamples; ++iter)
    sched.execute_plan(eplan);
}

