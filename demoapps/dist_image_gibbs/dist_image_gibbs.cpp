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
 * unused
 */
typedef char edge_data;


/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data {
  size_t sample;  
  graphlab::unary_factor unary;
  void save(graphlab::oarchive &oarc) const {
    oarc << unary << sample;
  }
  
  void load(graphlab::iarchive &iarc) {
    iarc >> unary >> sample;
  }
}; // End of vertex data


typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::distributed_types<graph_type> gl_types;

gl_types::distributed_glshared<graphlab::binary_factor> EDGE_FACTOR;


// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl_types::memory_graph& distributed_graph);

/** 
 * The core gibbs update function.  This update satisfies
 * the graphlab update_function interface.  
 */
void gibbs_update(gl_types::iscope& scope, 
                  gl_types::icallback& scheduler);
               


// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using gibbs sampling inside " << std::endl
            << "the graphlab framework." << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

  bool makegraph = false;
  size_t colors = 5;
  size_t rows = 200;
  size_t cols = 200;
  double sigma = 1;
  double lambda = 2;
  std::string smoothing = "square";
  std::string orig_fn = "source_img.pgm";
  std::string noisy_fn = "noisy_img.pgm";
  std::string pred_fn = "pred_img.pgm";




  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Loopy BP image denoising");
  clopts.use_distributed_options();
  clopts.attach_option("makegraph",
                       &makegraph, makegraph,
                       "Creates the disk distributed_graph");
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
  

  clopts.set_scheduler_type("multiqueue_fifo");
  clopts.set_scope_type("edge");
  

  bool success = clopts.parse(argc, argv);
  if(!success) {    
    return EXIT_FAILURE;
  }


  
  std::cout << "ncpus:          " << clopts.get_ncpus() << std::endl
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
            << "pred_fn:        " << pred_fn << std::endl;

  

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
    gl_types::disk_graph dg("denoise", 64);
    gl_types::memory_graph g;

    construct_graph(img, colors, sigma, g);
    std::vector<graphlab::graph_partitioner::part_id_type> parts;
    graphlab::graph_partitioner::metis_partition(g, 64, parts);
    dg.create_from_graph(g, parts);
    dg.finalize();
    return 0;
  }

  graphlab::mpi_tools::init(argc, argv);
  
  graphlab::dc_init_param param;
  ASSERT_TRUE(graphlab::init_param_from_mpi(param));
  // create distributed control
  graphlab::distributed_control dc(param);
  // Create the distributed_graph --------------------------------------------------------->
  gl_types::distributed_core core(dc, "denoise.idx");
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
  


  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function",gibbs_update);

  std::cout << "Running the engine. " << std::endl;

  
  // Add the bp update to all vertices
  core.add_task_to_all(gibbs_update, 100.0);
  // Starte the engine
  double runtime = core.start();
  std::vector<vertex_data> all_vdata = core.graph().collect_vertices(0);
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
    for(size_t v = 0; v < all_vdata.size(); ++v) {
      const vertex_data& vdata = all_vdata[v];
      img.pixel(v) = vdata.sample;
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
void gibbs_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler) {

  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  
  // Get the in and out edges by reference
  gl_types::edge_list in_edges = scope.in_edge_ids();
 
  graphlab::unary_factor u = v_data.unary;
  graphlab::binary_factor efactor = EDGE_FACTOR.get_val();
  
  foreach(graphlab::edge_id_t ineid, in_edges) {   
    size_t target_val = scope.const_neighbor_vertex_data(scope.source(ineid)).sample;
    for (size_t i = 0;i < u.arity(); ++i) {
      u.logP(i) += efactor.logP(target_val, i);
    }
  }
  u.normalize();
  v_data.sample = u.sample();
} 


void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl_types::memory_graph& distributed_graph) {
  // Construct a single blob for the vertex data
  vertex_data vdata;
  vdata.unary.resize(num_rings);
  // Add all the vertices
  double sigmaSq = sigma*sigma;
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      // Set the node potential
      double obs = img.pixel(i, j);
      for(size_t pred = 0; pred < num_rings; ++pred) {
        vdata.unary.logP(pred) = 
          -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
      }
      vdata.unary.normalize();
      vdata.sample = graphlab::random::uniform<uint32_t>(0, num_rings - 1);
      // Store the actual data in the distributed_graph
      size_t vertid = distributed_graph.add_vertex(vdata);
      // Ensure that we are using a consistent numbering
      assert(vertid == img.vertid(i, j));
    } // end of for j in cols
  } // end of for i in rows

  // Add the edges
 
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        distributed_graph.add_edge(vertid, img.vertid(i-1, j), edge_data());
      }
      if(i+1 < img.rows()) {
        distributed_graph.add_edge(vertid, img.vertid(i+1, j), edge_data());
      }
      if(j-1 < img.cols()) {
        distributed_graph.add_edge(vertid, img.vertid(i, j-1), edge_data());
      } if(j+1 < img.cols()) {
        distributed_graph.add_edge(vertid, img.vertid(i, j+1), edge_data());
      }
    } // end of for j in cols
  } // end of for i in rows
  distributed_graph.compute_coloring();
  distributed_graph.finalize();  
} // End of construct distributed_graph


