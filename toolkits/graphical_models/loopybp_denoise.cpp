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


#include <graphlab.hpp>

// #include "image.hpp"
#include "factors/factor_includes.hpp"


#include <cv.h>
#include <highgui.h>  


// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>


// Global variables
binary_factor EDGE_FACTOR;
size_t NCOLORS;
double SIGMA;
double BOUND;
double DAMPING;



// STRUCTS (Edge and Vertex data) =============================================>

/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data : graphlab::IS_POD_TYPE {
  float obs_color;
  uint16_t true_color, pred_color;
  vertex_data(float obs_color = 0, uint16_t true_color = 0) : 
    obs_color(obs_color), true_color(true_color), pred_color(obs_color) { }
}; // End of vertex data


/**
 * The data associated with each directed edge in the pairwise markov
 * random field
 */
class edge_data {
  unary_factor messages[4];
  size_t message_idx(size_t source_id, size_t target_id, bool is_new) {
    return size_t(source_id < target_id)  + 2 * size_t(is_new);
  }
public:
  edge_data() { }
  edge_data(graphlab::vertex_id_type v1,
            graphlab::vertex_id_type v2,
            size_t ncolors) {
    for(size_t i = 0; i < 4; ++i) {
      messages[i].resize(ncolors);
      messages[i].uniform();
    }
    message(v1, v2).var() = v2;
    old_message(v1, v2).var() = v2;
    message(v2, v1).var() = v1;
    old_message(v2, v1).var() = v1;
  } // end of constructor

  unary_factor& message(size_t source_id, size_t target_id) { 
    return messages[message_idx(source_id, target_id, true)];
  }
  unary_factor& old_message(size_t source_id, size_t target_id) { 
    return messages[message_idx(source_id, target_id, false)];
  }
  void update_old(size_t source_id, size_t target_id) { 
    old_message(source_id, target_id) = message(source_id, target_id);
  }
  void save(graphlab::oarchive& arc) const {
    for(size_t i = 0; i < 4; ++i) arc << messages[i];
  }
  void load(graphlab::iarchive& arc) {
    for(size_t i = 0; i < 4; ++i) arc >> messages[i];
  }
}; // End of edge data



/**
 * The type of the distributed graph representing the MRF.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

/** 
 * The gather_type for the vertex program needs to compute *= in place
 * of += so we create a new type which convertes computes *= for +=.
 */
struct factor_product {
  unary_factor factor;
  factor_product(const unary_factor& factor = unary_factor()) : 
    factor(factor) { }
  factor_product& operator+=(const factor_product& other) {
    ASSERT_EQ(factor.arity(), other.factor.arity());
    factor *= other.factor;
    return *this;
  }
  void save(graphlab::oarchive& arc) const { arc << factor; }
  void load(graphlab::iarchive& arc) { arc >> factor; }
}; // end of struct factor product


/** 
 * Belief Propagation Vertex Program
 *
 */
class bp_vertex_program : 
  public graphlab::ivertex_program< graph_type, factor_product,
                                    graphlab::messages::sum_priority > {
private:
  /**
   * The belief estimate for this vertex program
   */
  unary_factor belief;
public:
  void save(graphlab::oarchive& arc) const { arc << belief; }
  void load(graphlab::iarchive& arc) { arc >> belief; }

  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /**
   * Update the old message to be the new message and collect the
   * message value.
   */
  gather_type gather(icontext_type& context, 
                     const vertex_type& vertex, 
                     edge_type& edge) const {
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    edge_data& edata = edge.data();
    edata.update_old(other_vertex.id(), vertex.id());
    return factor_product(edata.old_message(other_vertex.id(), vertex.id()));
  }; // end of gather function

  /**
   * Multiply message product by node potential and update the belief.
   */
  void apply(icontext_type& context, vertex_type& vertex, 
             const gather_type& total) {
    // construct the node potential
    belief = make_potential(vertex);
    ASSERT_EQ(belief.arity(), total.factor.arity());
    // multiply in the rest of the message product;
    belief *= total.factor;
    belief.normalize();
    // compute the predicted value
    vertex.data().pred_color = belief.max_asg();
  }; // end of apply

  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  /**
   * Compute new message value for each edge.
   */
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {  
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    edge_data& edata = edge.data();
    // construct the cavity
    unary_factor cavity = belief;
    cavity /= edata.old_message(other_vertex.id(), vertex.id());
    cavity.normalize();
    // compute the new message
    unary_factor& new_message = 
      edata.message(vertex.id(), other_vertex.id());
    const unary_factor& old_message = 
      edata.old_message(vertex.id(), other_vertex.id());
    ASSERT_NE(&new_message, &old_message);
    new_message.convolve(EDGE_FACTOR, cavity);
    new_message.normalize();
    new_message.damp(old_message, DAMPING);
    // Compute message residual
    const double residual = new_message.residual(old_message);  
    context.clear_gather_cache(other_vertex);
    // Schedule the adjacent vertex
    if(residual > BOUND) context.signal(other_vertex, residual);
 }; // end of scatter

private:
  /**
   * Construct the unary evidence potential
   */
  unary_factor make_potential(const vertex_type& vertex) const {
    unary_factor potential(vertex.id(), NCOLORS);
    const double obs = vertex.data().obs_color;
    const double sigmaSq = SIGMA*SIGMA;
    for(size_t pred = 0; pred < potential.arity(); ++pred) {
      potential.logP(pred) = 
        -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
    }
    potential.normalize();
    return potential;
  } // end of make_potentail

  /**
   * Return the other vertex
   */
  vertex_type get_other_vertex(edge_type& edge, 
                               const vertex_type& vertex) const {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
  }; // end of other_vertex

}; // end of class bp_vertex_program

/**
 * Define the engine type
 */
//typedef graphlab::synchronous_engine<bp_vertex_program> engine_type;
typedef graphlab::async_consistent_engine<bp_vertex_program> engine_type;


/**
 * construct the synthetic image graph.
 */
void create_synthetic_mrf(graphlab::distributed_control& dc,
                          graph_type& graph, 
                          const size_t rows, const size_t cols); 

               
template<typename T>
struct merge_reduce {
  std::vector<T> values;
  void save(graphlab::oarchive& arc) const { arc << values; }
  void load(graphlab::iarchive& arc) { arc >> values; }
  merge_reduce& operator+=(const merge_reduce& other) {
    values.insert(values.end(), other.values.begin(), 
                  other.values.end());
    return *this;
  }
}; // end of merge_reduce

typedef std::pair<graphlab::vertex_id_type, float> pred_pair_type; 
typedef merge_reduce<pred_pair_type> merge_reduce_type;

merge_reduce_type pred_map_function(graph_type::vertex_type vertex) {
  merge_reduce<pred_pair_type> ret;
  ret.values.push_back(pred_pair_type(vertex.id(), vertex.data().pred_color));
  return ret;
} // end of pred_map_function

merge_reduce_type obs_map_function(graph_type::vertex_type vertex) {
  merge_reduce<pred_pair_type> ret;
  ret.values.push_back(pred_pair_type(vertex.id(), vertex.data().obs_color));
  return ret;
} // end of obs_map_function


/**
 * Save the image data in the vector of pairs to an image file
 */
void save_image(const size_t rows, const size_t cols,
                const std::vector<pred_pair_type>& values,
                const std::string& fname);



// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using loopy belief propagation inside " << std::endl
            << "the graphlab framework." << std::endl;

  // // set the global logger
  // global_logger().set_log_level(LOG_WARNING);
  // global_logger().set_log_to_console(true);

  // Set initial values for members ------------------------------------------->
  NCOLORS = 5;
  SIGMA = 2;
  BOUND = 1E-4;
  DAMPING = 0.1;
 
  size_t nrows = 200;
  size_t ncols = 200;
  double lambda = 2;

  std::string smoothing = "laplace";
  std::string orig_fn = "source_img.jpeg";
  std::string noisy_fn = "noisy_img.jpeg";
  std::string pred_fn = "pred_img.jpeg";

  // std::string orig_fn = "source_img.pgm";
  // std::string noisy_fn = "noisy_img.pgm";
  // std::string pred_fn = "pred_img.pgm";



  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Loopy BP image denoising");
  clopts.attach_option("bound",
                       &BOUND, BOUND,
                       "Residual termination bound");
  clopts.attach_option("damping",
                       &DAMPING, DAMPING,
                       "The amount of message damping (higher = more damping)");
  clopts.attach_option("ncolors",
                       &NCOLORS, NCOLORS,
                       "The number of colors in the noisy image");
  clopts.attach_option("sigma",
                       &SIGMA, SIGMA,
                       "Standard deviation of noise.");
  clopts.attach_option("nrows",
                       &nrows, nrows,
                       "The number of rows in the noisy image");
  clopts.attach_option("ncols",
                       &ncols, ncols,
                       "The number of columns in the noisy image");
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
    

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  const bool success = clopts.parse(argc, argv);
  if(!success) {
    graphlab::mpi_tools::finalize();
    return EXIT_FAILURE;
  }
  ///! Create a distributed control object 
  graphlab::distributed_control dc;
  ///! display settings  
  if(dc.procid() == 0) {
    std::cout << "ncpus:          " << clopts.get_ncpus() << std::endl
              << "bound:          " << BOUND << std::endl
              << "damping:        " << DAMPING << std::endl
              << "colors:         " << NCOLORS << std::endl
              << "nrows:           " << nrows << std::endl
              << "ncols:           " << ncols << std::endl
              << "sigma:          " << SIGMA << std::endl
              << "lambda:         " << lambda << std::endl
              << "smoothing:      " << smoothing << std::endl
              << "scheduler:      " << clopts.get_scheduler_type() << std::endl
              << "orig_fn:        " << orig_fn << std::endl
              << "noisy_fn:       " << noisy_fn << std::endl
              << "pred_fn:        " << pred_fn << std::endl;
  }

  
  

  // Create synthetic images -------------------------------------------------->
  std::cout << "Creating a synthetic noisy image." << std::endl;
  graph_type graph(dc, clopts);
  create_synthetic_mrf(dc, graph, nrows, ncols);
  std::cout << "Finalizing the graph." << std::endl;
  graph.finalize();
  if(dc.procid() == 0) {
    std::cout << "Number of vertices: " << graph.num_vertices() << std::endl
              << "Number of edges:    " << graph.num_edges() << std::endl;
  }

  std::cout << "Collect the noisy image. " << std::endl;
  merge_reduce_type obs_image = 
    graph.map_reduce_vertices<merge_reduce_type>(obs_map_function);
  std::cout << "saving the noisy image." << std::endl;
  if(dc.procid() == 0) {
    save_image(nrows, ncols, obs_image.values, noisy_fn);
  }

  // Initialze the edge factor ----------------------------------------------->
  std::cout << "Initializing shared edge agreement factor. " << std::endl;
  // dummy variables 0 and 1 and num_rings by num_rings
  EDGE_FACTOR = binary_factor(0, NCOLORS, 0, NCOLORS);
  // Set the smoothing type
  if(smoothing == "square") {
    EDGE_FACTOR.set_as_agreement(lambda);
  } else {
    EDGE_FACTOR.set_as_laplace(lambda);
  } 
  if(dc.procid() == 0)
    std::cout << EDGE_FACTOR << std::endl;

  // Create the engine -------------------------------------------------------->
  std::cout << "Creating the engine. " << std::endl;
  engine_type engine(dc, graph, clopts);

  std::cout << "Scheduling all vertices" << std::endl;
  engine.signal_all();
  std::cout << "Starting the engine" << std::endl;
  engine.start();
  const float runtime = engine.elapsed_seconds();
  size_t update_count = engine.num_updates();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;  


  // Saving the output -------------------------------------------------------->
  std::cout << "Saving the predicted image" << std::endl;
  std::cout << "Collect the noisy image. " << std::endl;
  merge_reduce_type pred_image = 
    graph.map_reduce_vertices<merge_reduce_type>(pred_map_function);
  std::cout << "saving the pred image." << std::endl;
  if(dc.procid() == 0) {
    save_image(nrows, ncols, pred_image.values, pred_fn);
  }

  std::cout << "Done!" << std::endl;
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main




graphlab::vertex_id_type sub2ind(size_t rows, size_t cols,
                                 size_t r, size_t c) {
  return r * cols + c;
}; // end of sub2ind

std::pair<int,int> ind2sub(size_t rows, size_t cols,
                           size_t ind) {
  return std::make_pair(ind / cols, ind % cols);
}; // end of sub2ind


void create_synthetic_mrf(graphlab::distributed_control& dc,
                          graph_type& graph,
                          const size_t rows, const size_t cols) {
  dc.barrier();
  const double center_r = rows / 2.0;
  const double center_c = cols / 2.0;
  const double max_radius = std::min(rows, cols) / 2.0;
 
  for(size_t r = dc.procid(); r < rows; r += dc.numprocs()) {
    for(size_t c = 0; c < cols; ++c) {
      // Compute the true pixel value
      const double distance = sqrt((r-center_r)*(r-center_r) + 
                                   (c-center_c)*(c-center_c));
      // Compute ring of sunset
      const uint16_t ring_color =  
        std::floor(std::min(1.0, distance/max_radius) * (NCOLORS - 1) );
      // Compute the true pixel color by masking with the horizon
      const uint16_t true_color = r < rows/2 ? ring_color : 0;
      // compute the predicted color
      const float obs_color = true_color + graphlab::random::normal(0, SIGMA);
      // determine the true pixel id
      const graphlab::vertex_id_type vid = sub2ind(rows,cols,r,c);
      const vertex_data vdata(obs_color, true_color);
      graph.add_vertex(vid, vdata);
      // Add the edges
      if(r + 1 < rows) 
        graph.add_edge(vid, sub2ind(rows,cols,r+1,c),
                       edge_data(vid, sub2ind(rows,cols,r+1,c), NCOLORS));
      if(c + 1 < cols) 
        graph.add_edge(vid, sub2ind(rows,cols,r,c+1),
                       edge_data(vid, sub2ind(rows,cols,r,c+1), NCOLORS));
    } // end of loop over cols
  } // end of loop over rows
  dc.barrier();
}; // end of create synthetic mrf


// void save_image(const size_t rows, const size_t cols,
//                 const std::vector<pred_pair_type>& values,
//                 const std::string& fname) {
//   std::cout << "NPixels: " << values.size() << std::endl;
//   image img(rows, cols);
//   foreach(pred_pair_type pair, values) 
//     img.pixel(pair.first) = pair.second;
//   img.save(fname);
// }


void save_image(const size_t rows, const size_t cols,
                const std::vector<pred_pair_type>& values,
                const std::string& fname) {
  std::cout << "NPixels: " << values.size() << std::endl;
  // determine the max and min colors
  float max_color = -std::numeric_limits<float>::max();
  float min_color =  std::numeric_limits<float>::max();
  foreach(pred_pair_type pair, values) {
    max_color = std::max(max_color, pair.second);
    min_color = std::min(min_color, pair.second);
  }

  cv::Mat img(cols, rows, CV_8UC1);
  foreach(pred_pair_type pair, values) {
    std::pair<int,int> coords = ind2sub(rows,cols, pair.first);
    float value = (pair.second - min_color) / (max_color - min_color);
    int color = 255 * value > 255 ? 255 : 255 * value;
    img.at<unsigned char>(coords.first, coords.second) = color;
  }
  cv::imwrite(fname, img);
}
