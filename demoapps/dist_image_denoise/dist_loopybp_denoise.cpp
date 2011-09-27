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


#include <distributed_graphlab.hpp>

#include "image.hpp"


// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>




// STRUCTS (Edge and Vertex data) =============================================>

/**
 * The data associated with each directed edge in the pairwise markov
 * random field
 */
struct edge_data {
  graphlab::unary_factor message;
  graphlab::unary_factor old_message;
  void save(graphlab::oarchive &oarc) const {
    oarc << message << old_message;
  }
  
  void load(graphlab::iarchive &iarc) {
    iarc >> message >> old_message;
  }
}; // End of edge data


/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data {
  graphlab::unary_factor potential;
  graphlab::unary_factor belief;
  void save(graphlab::oarchive &oarc) const {
    oarc << potential << belief;
  }
  
  void load(graphlab::iarchive &iarc) {
    iarc >> potential >> belief;
  }
}; // End of vertex data


typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::distributed_types<graph_type> gl_types;

gl_types::distributed_glshared<graphlab::binary_factor> EDGE_FACTOR;
gl_types::distributed_glshared<double> BOUND;
gl_types::distributed_glshared<double> DAMPING;


// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl_types::memory_graph& distributed_graph);

/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler);
               

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


// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using loopy belief propagation inside " << std::endl
            << "the graphlab framework." << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);

  bool makegraph = false;
  double bound = 1E-4;
  double damping = 0.1;
  size_t colors = 5;
  size_t rows = 200;
  size_t cols = 200;
  double sigma = 2;
  double lambda = 2;
  std::string smoothing = "laplace";
  std::string orig_fn = "source_img.pgm";
  std::string noisy_fn = "noisy_img.pgm";
  std::string pred_fn = "pred_img.pgm";
  std::string pred_type = "map";




  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Loopy BP image denoising");
  clopts.use_distributed_options();
  clopts.attach_option("makegraph",
                       &makegraph, makegraph,
                       "Creates the disk distributed_graph");
  clopts.attach_option("bound",
                       &bound, bound,
                       "Residual termination bound");
  clopts.attach_option("damping",
                       &damping, damping,
                       "The amount of message damping (higher = more damping)");
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
  clopts.attach_option("smoothing",
                       &smoothing, smoothing,
                       "Options are {square, laplace}");
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
  

  clopts.set_scheduler_type("multiqueue_fifo");
  clopts.set_scope_type("edge");
  

  bool success = clopts.parse(argc, argv);
  if(!success) {    
    return EXIT_FAILURE;
  }


  
  std::cout << "ncpus:          " << clopts.get_ncpus() << std::endl
            << "bound:          " << bound << std::endl
            << "damping:        " << damping << std::endl
            << "colors:         " << colors << std::endl
            << "rows:           " << rows << std::endl
            << "cols:           " << cols << std::endl
            << "sigma:          " << sigma << std::endl
            << "lambda:         " << lambda << std::endl
            << "smoothing:      " << smoothing << std::endl
            << "engine:         " << clopts.get_engine_type() << std::endl
            << "scope:          " << clopts.get_scope_type() << std::endl
            << "scheduler:      " << clopts.get_scheduler_type() << std::endl
            << "orig_fn:        " << orig_fn << std::endl
            << "noisy_fn:       " << noisy_fn << std::endl
            << "pred_fn:        " << pred_fn << std::endl
            << "pred_type:      " << pred_type << std::endl;

  

  std::cout << "Creating a synthetic image. " << std::endl;
  image img(rows, cols);
  img.paint_sunset(colors);
  std::cout << "Saving image. " << std::endl;
  img.save(orig_fn.c_str());
  std::cout << "Corrupting Image. " << std::endl;
  img.corrupt(sigma);

  if (makegraph) {
    // Create synthetic images -------------------------------------------------->
    // Creating image for denoising
    std::cout << "Saving corrupted image. " << std::endl;
    img.save(noisy_fn.c_str());
    
    std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
    gl_types::disk_graph dg("denoise", 32, graphlab::disk_graph_atom_type::WRITE_ONLY_ATOM);
    gl_types::memory_graph g;

    construct_graph(img, colors, sigma, g);
    std::vector<graphlab::graph_partitioner::part_id_type> parts;
    graphlab::graph_partitioner::metis_partition(g, 32, parts);
    dg.create_from_graph(g, parts);
    dg.make_memory_atoms();
    dg.finalize();
    return 0;
  }

  graphlab::mpi_tools::init(argc, argv);
  
  graphlab::dc_init_param param;
  if(graphlab::init_param_from_env(param) == false) {
    graphlab::init_param_from_mpi(param); 
  }
  // create distributed control
  graphlab::distributed_control dc(param);
  // Create the distributed_graph --------------------------------------------------------->
  gl_types::distributed_core core(dc, "denoise.idx", graphlab::disk_graph_atom_type::WRITE_ONLY_ATOM);
  // Set the engine options
  core.set_engine_options(clopts);
  core.build_engine();

  
  // Setup global shared variables -------------------------------------------->
  // Initialize the edge agreement factor 
  std::cout << "Initializing shared edge agreement factor. " << std::endl;

  // dummy variables 0 and 1 and num_rings by num_rings
  graphlab::binary_factor edge_potential(0, colors, 0, colors);
  // Set the smoothing type
  if(smoothing == "square") {
    edge_potential.set_as_agreement(lambda);
  } else if (smoothing == "laplace") {
    edge_potential.set_as_laplace(lambda);
  } else {
    std::cout << "Invalid smoothing stype!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << edge_potential << std::endl;
  
  EDGE_FACTOR.set(edge_potential);
  BOUND.set(bound);
  DAMPING.set(damping);
  


  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function",bp_update);

  std::cout << "Running the engine. " << std::endl;

  
  // Add the bp update to all vertices
  core.add_task_to_all(bp_update, 100.0);
  // Starte the engine
  double runtime = core.start();
  if (dc.procid() == 0) {
    size_t update_count = core.last_update_count();
    std::cout << "Finished Running engine in " << runtime 
              << " seconds." << std::endl
              << "Total updates: " << update_count << std::endl
              << "Efficiency: " << (double(update_count) / runtime)
              << " updates per second "
              << std::endl;  


    // Saving the output -------------------------------------------------------->
    std::cout << "Rendering the cleaned image. " << std::endl;
    if(pred_type == "map") {
      for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
        const vertex_data& vdata = core.graph().get_vertex_data(v);
        img.pixel(v) = vdata.belief.max_asg();    
      }
    } else if(pred_type == "exp") {
      for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
        const vertex_data& vdata = core.graph().get_vertex_data(v);
        img.pixel(v) = vdata.belief.expectation();
      }
    } else {
      std::cout << "Invalid prediction type! : " << pred_type
                << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Saving cleaned image. " << std::endl;
    img.save(pred_fn.c_str());

    std::cout << "Done!" << std::endl;
  }
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main




// Implementations
// ============================================================>
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler) {
  //  std::cout << scope.vertex();;
  //  std::getchar();

  // Get the shared data
  double bound = BOUND.get_val();
  double damping = DAMPING.get_val();

  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  
  // Get the in and out edges by reference
  gl_types::edge_list in_edges = scope.in_edge_ids();
  gl_types::edge_list out_edges = scope.out_edge_ids();
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
    const edge_data& e_data = scope.edge_data(ineid);
    // Notice we now use the old message since neighboring vertices
    // could be changing the new messages
    v_data.belief.times( e_data.old_message );
  }
  v_data.belief.normalize(); // finally normalize the belief
  
  // Compute outbound messages
  // ---------------------------------------------------------------->

  const graphlab::binary_factor edge_factor = EDGE_FACTOR.get_val();
  
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
    const edge_data& in_edge = scope.edge_data(ineid);
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
                     gl_types::memory_graph& distributed_graph) {
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
      // Store the actual data in the distributed_graph
      size_t vertid = distributed_graph.add_vertex(vdata);
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
        distributed_graph.add_edge(vertid, img.vertid(i-1, j), edata);
      }
      if(i+1 < img.rows()) {
        edata.message.var() = img.vertid(i+1, j);
        edata.old_message.var() = edata.message.var();
        distributed_graph.add_edge(vertid, img.vertid(i+1, j), edata);
      }
      if(j-1 < img.cols()) {
        edata.message.var() = img.vertid(i, j-1);
        edata.old_message.var() = edata.message.var();
        distributed_graph.add_edge(vertid, img.vertid(i, j-1), edata);
      } if(j+1 < img.cols()) {
        edata.message.var() = img.vertid(i, j+1);
        edata.old_message.var() = edata.message.var();
        distributed_graph.add_edge(vertid, img.vertid(i, j+1), edata);
      }
    } // end of for j in cols
  } // end of for i in rows
  distributed_graph.finalize();  
} // End of construct distributed_graph


