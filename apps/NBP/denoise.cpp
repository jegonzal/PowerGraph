/**
 * This file contains an example of graphlab used for discrete loopy
 * belief propagation in a pairwise markov random field to denoise a
 * synthetic noisy image.
 *
 *  \author Yucheng
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <iostream>
#include <map>

#include <graphlab.hpp>
#include <limits>
#include "kde.h"
#include <float.h>
#include "image.hpp"

// include itpp
#include <itpp/itstat.h>
#include <itpp/itbase.h>

// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

#define PROPOSAL_STDEV 20

int NSAMP =12;
double EPSILON =1e-5;
int RESAMPLE_FREQUENCY = 0;
int MAX_ITERATIONS = 2;

int ROWS=107;
int COLS=86;


using namespace itpp;
using namespace std;

struct particle {
  float x;
  float weight;
};

float SIGMA;
float LAMBDA;
float damping = 0.8;
const size_t MCMCSTEPS = 30;
const size_t RESAMPLE_FREQUENCY = 5;
size_t MAX_ITERATIONS;
graphlab::atomic<size_t> proposal_count;
graphlab::atomic<size_t> accept_count;

// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
//  vec message;
//  vec old_message;
  kde msg;
  kde edge_pot;
  // vec lfeat_message;
  int update_count;
}; // End of edge data

struct vertex_data: public graphlab::unsupported_serialize {
  kde obs;
  kde bel;
  size_t row, col;
  int rounds;

  vertex_data(){ rounds = 0;}
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;



// Implementations
// ============================================================>
void bp_update(gl_types::iscope& scope,
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data) {


  bool debug = true;

  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  graphlab::vertex_id_t vid = scope.vertex();
  if (debug && vid%100 == 0){
     std::cout<<"Entering node " << (int)vid << " obs: ";
     v_data.obs.matlab_print();
     std::cout << std::endl;
  }

  v_data.rounds++;

  if ((int)vid == 0)
      iiter++;
  //vec prod_message = prod_msg0.get_col(vid);
  //v_data.belief = belief0.get_col(vid);

  // Get the in and out edges by reference
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size()); // Sanity check


  for (size_t j = 0; j < in_edges.size(); ++j){
   
     std::vector<kde> kdes;
     for(size_t i = 0; i < in_edges.size(); ++i) {
    
     // Get the edge ids
      graphlab::edge_id_t ineid = in_edges[i];
      
    // const edge_data& in_edge = scope.edge_data(ineid);
      edge_data& in_edge = scope.edge_data(ineid);
      if (i != j){
         in_edge.msg.verify();
         kdes.push_back(in_edge.msg);
      }
    }

     kdes.push_back(v_data.obs);
     
      graphlab::edge_id_t outeid = out_edges[j];
      edge_data& out_edge = scope.edge_data(outeid);
      kde marg = out_edge.edge_pot.marginal(0);  
      kdes.insert(kdes.begin(), marg);//important: has to be first!     

      prodSampleEpsilon producter; 
      kde m = producter.prodSampleEpsilonRun(kdes.size(), 
                                     NSAMP, 
                                     EPSILON, 
                                     kdes);
     
      m.verify();
      kde mar2 = out_edge.edge_pot.marginal(1);
      mar2.verify();
      imat firstrowind = m.indices(0,0,0,m.indices.cols()-1);
      kde outmsg = mar2.sample(firstrowind,m.weights);
      outmsg.verify(); 
      out_edge.msg = outmsg;

    }


   //compute belief
   if (v_data.rounds == MAX_ITERATIONS){
	if (debug && vid%100 == 0)
	   printf("computing belief node %d\n", vid);

      std::vector<kde> kdes;
    

      for (size_t j = 0; j < in_edges.size(); ++j){
    
         // Get the edge ids
        graphlab::edge_id_t ineid = in_edges[j];
      
        edge_data& in_edge = scope.edge_data(ineid);
        in_edge.msg.verify();
        kdes.push_back(in_edge.msg);
      }

      kdes.push_back(v_data.obs);
      prodSampleEpsilon prod;
      kde m = prod.prodSampleEpsilonRun(kdes.size(), 
                                     NSAMP, 
                                     EPSILON, 
                                     kdes);
     
      m.verify();
      v_data.bel = m;
      if (debug && vid == 0){
	   printf("belief node %d is\n", vid);
           m.matlab_print(); printf("\n");
      }


   }
  /*
  if (v_data.rounds < MAX_ITERATIONS) {
    gl_types::update_task task(scope.vertex(), bp_update);
    scheduler.add_task(task, 1.0);
  }*/


} // end of BP_update





void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     double numparticles,
                     gl_types::graph& graph) {
  // initialize a bunch of particles


  
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      vertex_data vdat;
      vdat.rounds = 0;

      // Set the node potential
      vdat.obs = kde(img.pixel(i, j),30,1);;
      vdat.belief = vdat.obs;
      graph.add_vertex(vdat);
      vdat.obs.verify();
      vdat.bel.verify();

    } // end of for j in cols
  } // end of for i in rows

  // Add the edges
  edge_data edata;
  

  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      
      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        edata.msg = graph.vertex_data(img.vertid(i-1, j)).bel;
        graph.add_edge(vertid, img.vertid(i-1, j), edata);
      }
      if(i+1 < img.rows()) {
        edata.msg = graph.vertex_data(img.vertid(i+1, j)).bel;
        graph.add_edge(vertid, img.vertid(i+1, j), edata);
      }
      if(j-1 < img.cols()) {
        edata.msg = graph.vertex_data(img.vertid(i, j-1)).bel;
        graph.add_edge(vertid, img.vertid(i, j-1), edata);
      }
      if(j+1 < img.cols()) {
        edata.msg = graph.vertex_data(img.vertid(i, j+1)).bel;
        graph.add_edge(vertid, img.vertid(i, j+1), edata);
      }
    } // end of for j in cols
  } // end of for i in rows
  graph.finalize();
} // End of construct graph




















// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using loopy belief propagation inside " << std::endl
            << "the graphlab framework." << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  size_t iterations = 100;
  size_t colors = 5;
  size_t rows = 200;
  size_t cols = 200;
  double sigma = 2;
  double lambda = 10;
  size_t numparticles = 100;
    std::string orig_fn = "source_img.pgm";
  std::string noisy_fn = "noisy_img.pgm";
  std::string pred_fn = "pred_img.pgm";
  std::string pred_type = "map";




  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Loopy BP image denoising");
  clopts.attach_option("iterations",
                       &iterations, iterations,
                       "Number of iterations");
  clopts.attach_option("particles",
                       &numparticles, numparticles,
                       "Number of particlesw");
  clopts.attach_option("epsilon",
                       &EPSILON, EPSILON,
                       "epsilon");
  clopts.attach_option("colors",
                       &colors, colors,
                       "The number of colors in the noisy image");
  clopts.attach_option("rows",
                       &rows, rows,
                       "The number of rows in the noisy image");
  clopts.attach_option("cols",
                       &cols, cols,
                       "The number of columns in the noisy image");
  clopts.attach_option("sigma",
                       &sigma, sigma,
                       "Standard deviation of noise.");
  clopts.attach_option("lambda",
                       &lambda, lambda,
                       "Smoothness parameter (larger => smoother).");
  clopts.attach_option("orig",
                       &orig_fn, orig_fn,
                       "Original image file name.");
  clopts.attach_option("noisy",
                       &noisy_fn, noisy_fn,
                       "Noisy image file name.");
  clopts.attach_option("pred",
                       &pred_fn, pred_fn,
                       "Predicted image file name.");
  clopts.attach_option("pred_type",
                       &pred_type, pred_type,
                       "Predicted image type {map, exp}");
  clopts.attach_option("damping",
                       &damping, damping,
                       "damping");


  clopts.scheduler_type = "round_robin";
  clopts.scope_type = "none";


  bool success = clopts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }

  // fill the global vars
  SIGMA = sigma;
  LAMBDA = lambda;
  MAX_ITERATIONS = iterations;

  std::cout << "ncpus:          " << clopts.ncpus << std::endl
            << "iterations:          " << iterations<< std::endl
            << "particles:          " << numparticles<< std::endl
            << "colors:         " << colors << std::endl
            << "rows:           " << rows << std::endl
            << "cols:           " << cols << std::endl
            << "sigma:          " << sigma << std::endl
            << "lambda:         " << lambda << std::endl
            << "engine:         " << clopts.engine_type << std::endl
            << "scope:          " << clopts.scope_type << std::endl
            << "scheduler:      " << clopts.scheduler_type << std::endl
            << "orig_fn:        " << orig_fn << std::endl
            << "noisy_fn:       " << noisy_fn << std::endl
            << "pred_fn:        " << pred_fn << std::endl
            << "pred_type:      " << pred_type << std::endl;




  // Create synthetic images -------------------------------------------------->
  // Creating image for denoising
  std::cout << "Creating a synethic image. " << std::endl;
  image img(rows, cols);
  img.paint_sunset(colors);
  std::cout << "Saving image. " << std::endl;
  img.save(orig_fn.c_str());
  std::cout << "Corrupting Image. " << std::endl;
  img.corrupt(sigma);
  std::cout << "Saving corrupted image. " << std::endl;
  img.save(noisy_fn.c_str());





  // Create the graph --------------------------------------------------------->
  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);

  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  
  construct_graph(img, colors, sigma, numparticles, core.graph());



  // Running the engine ------------------------------------------------------->
  core.scheduler().set_option(gl_types::scheduler_options::UPDATE_FUNCTION,
                              (void*)bp_update);
  core.scheduler().set_option(gl_types::scheduler_options::MAX_ITERATIONS,
                              (void*)MAX_ITERATIONS);

  std::cout << "Running the engine. " << std::endl;


  // Add the bp update to all vertices
  core.add_task_to_all(bp_update, 100.0);
  // Starte the engine
  double runtime = core.start();

  size_t update_count = core.last_update_count();
  std::cout << "Finished Running engine in " << runtime
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;


  // Saving the output -------------------------------------------------------->
  std::cout << "Rendering the cleaned image. " << std::endl;
   
  for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    float a = vdata.max_asg();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = size_t(a);
  }
  img.save("pred_map.pgm");


  for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    float a = vdata.average();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = (a);
  }
  img.save("pred_exp.pgm");

  std::cout << "Saving cleaned image. " << std::endl;
  img.save(pred_fn.c_str());

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
} // End of main

