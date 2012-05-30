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
 * This file contains an example of graphlab used for MAP inference 
 * in a discrete graphical model (pairwise MRF). The algorithm
 * implemented is the MPLP LP-Relaxation scheme of Globerson & Jaakkola. 
 *
 *  \author Dhruv Batra
 */


#include <vector>
#include <string>
#include <fstream>


#include <Eigen/Dense>

#include <graphlab.hpp>


#include "eigen_serialization.hpp"

#include <graphlab/macros_def.hpp>

typedef Eigen::VectorXd vector;
typedef Eigen::MatrixXd matrix;

// Global variables
size_t NCOLORS;
double SIGMA;
double BOUND;
double DAMPING;

// Shared base edge potential
matrix thetaij; 

// STRUCTS (Edge and Vertex data) =============================================>

/**
 * Each GraphLab vertex is a (pairwise) factor from the MRF
 */
struct vertex_data {
  /** variable ids */
  int i, j; 
  /** observed color for each variable */
  float obs_color_i, obs_color_j;
  /** dual variables being optimized (or messages) */
  vector del_fi, del_fj;  
  // constructor
  vertex_data(): i(-1), j(-1), obs_color_i(-1), obs_color_j(-1) { }
  void save(graphlab::oarchive& arc) const {
    arc << i << j << obs_color_i << obs_color_j << del_fi << del_fj;
  }
  void load(graphlab::iarchive& arc) {
    arc >> i >> j >> obs_color_i >> obs_color_j >> del_fi >> del_fj;
  }
}; // End of vertex data


/**
 * The data associated with a pair of factors in a pairwise MRF
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  // primal labelling; We assume pairwise factors, so intersection has
  // a single node
  int pred_label; 
}; // End of edge data

/**
 * The graph type
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


// GraphLab Vertex Program ====================================================
/**
 * The type passed around during the gather phase
 */
struct gather_type {
  vector del_fi, del_fj;
  gather_type& operator+=(const gather_type& other) {
    if(!other.del_fi.empty()) {
      if(del_fi.empty()) del_fi = other.del_fi;
      else del_fi += other.del_fi;
    }
    if(!other.del_fj.empty()) {
      if(del_fj.empty()) del_fj = other.del_fj;
      else del_fj += other.del_fj;
    }
  } // end of operator +=
  void save(graphlab::oarchive& arc) const {
    arc << del_fi << del_fj;
  }
  void load(graphlab::iarchive& arc) {
    arc >> del_fi >> del_fj;
  }
}; // end of gather type



/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
class mplp_vertex_program : 
  public ivertex_program<graph_type, gather_type, 
                         graphlab::messages::sum_priority> {
public:

  // void save(graphlab::oarchive& arc) const { /** save members */ }
  // void load(graphlab::iarchive& arc) { /** load members */ }

  /**
   * The init function is called once for each vertex before the
   * start of the GraphLab program.  If the vertex program does not
   * implement this function then the default implementation (NOP)
   * is used.
   */
  // void init(icontext_type& context, vertex_type& vertex) { /** NOP */ }



  // I assume that MRF nodes are numbered 1-N. Now if i<k and i<j and
  // (i,j) (j,k) are two factors then I assume that GraphLab vertex
  // (i,j) will have an outgoing edge to (j,k) (edge are small to
  // large number).  This assumption let's me figure out which node i
  // or j to use in vertex data.  
  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /**
   * Since the MRF is undirected we will use all edges for gather and
   * scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 


  /**
   * Recv message is called by the engine to receive a message to this
   * vertex program.  The vertex program can use this to initialize
   * any state before entering the gather phase.  If the vertex
   * program does not implement this function then the default
   * implementation (NOP) is used.
   */
  void recv_message(icontext_type& context, const vertex_type& vertex, 
                    const message_type& msg) { /** NOP */ }

   
  // Run the gather operation over all in edges
  gather_type gather(icontext_type& context, const vertex_type& target_vertex, 
                     edge_type& edge) const {
    const vertex_type source_vertex = get_other_vertex(edge, vertex);
    edge_data& edata = edge.data();
    const vertex_data& source_vdata = source_vertex.data();
    const vertex_data& target_vdata = target_vertex.data();
    // Accumulate message
    gather_type ret_value;
    if (source_vdata.i == target_vdata.i)
      ret_value.del_fi = source_vdata.del_fi;
    else if (source_vdata.i == target_vdata.j)
      ret_value.del_fi = source_vdata.del_fj;
    else if (source_vdata.j == target_vdata.i)
      ret_value.del_fj = source_vdata.del_fi;
    else if (source_vdata.j == target_vdata.j)
      ret_value.del_fj = source_vdata.del_fj;
    else assert(false); // invalid state
    return ret_value;
  } // end of gather
    
  /** Update the dual parameters */
  void apply(icontext_type& context, vertex_type& vertex, 
             const gather_type& sum) {
    vertex_data& vdata = vertex.data();  
    const size_t arity = vdata.thetaij.rows();
    vector thetai = make_unary_potentail(vertex, 'i');
    vector thetaj = make_unary_potentail(vertex, 'j');
    // Update del fi
    vdata.del_fi = -(thetai + sum.del_fi)/2 +
      + (thetaij + sum.del_fj.rowwise().replicate(thetaij.rows())).
      rowwise().maxCoeff()/2;
    // Update del fj
    vdata.del_fj = -(thetaj + sum.del_fj)/2 +
      + ((thetaij + 
          sum.del_fi.rowwise().replicate(thetaij.cols()).transpose()).
         colwise().maxCoeff()).transpose()/2;       
  } // end of apply
  
  /** reschedule neighbors with a given priority and updated
      predictions on each edge*/
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {  
    // Nothing yet. Will hold the LPDG scheduling scheme. 
    const double priority = 1;
    context.signal(get_other_vertex(edge, vertex), priority);
  } // end of scatter

private:

  /**
   * Construct the unary evidence potential
   */
  unary_factor make_unary_potential(const vertex_type& vertex, 
                                    const char varid) const {    
    vector potential(NCOLORS);
    const double obs = varid == 'i'? 
      vertex.data().obs_color_i : vertex.data().obs_color_j;
    const double sigmaSq = SIGMA*SIGMA;
    for(size_t pred = 0; pred < potential.arity(); ++pred) {
      potential(pred) = -(obs - pred)*(obs - pred) / (2.0 * sigmaSq);
    }
    potential /= potential.sum();
    return potential;
  } // end of make_potentail

  /**
   * Return the other vertex
   */
  vertex_type get_other_vertex(edge_type& edge, 
                               const vertex_type& vertex) const {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
  }; // end of other_vertex

}; // end of MPLP vertex program


/**
 * Define the engine type
 */
typedef graphlab::synchronous_engine<bp_vertex_program> engine_type;


/**
 * construct the synthetic image graph.
 */
void create_synthetic_mrf(graphlab::distributed_control& dc,
                          graph_type& graph, 
                          const size_t rows, const size_t cols); 





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
  std::string orig_fn =  "source_img.pgm";
  std::string noisy_fn = "noisy_img.pgm";
  std::string pred_fn = "pred_img.pgm";



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
  engine.initialize();

  std::cout << "Scheduling all vertices" << std::endl;
  engine.signal_all();
  std::cout << "Starting the engine" << std::endl;
  engine.start();
  const float runtime = engine.elapsed_time();
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


void create_synthetic_mrf(graphlab::distributed_control& dc,
                          graph_type& graph,
                          const size_t rows, const size_t cols) {
  dc.barrier();
  // Generate the image on all machines --------------------------------------->
  // Need to ensure that all machines generate the same noisy image
  graphlab::random::generator gen; gen.seed(314);
  std::vector<float>    obs_pixels(rows * cols);
  std::vector<uint16_t> true_pixels(rows * cols);
  const double center_r = rows / 2.0;
  const double center_c = cols / 2.0;
  const double max_radius = std::min(rows, cols) / 2.0;
  for(size_t r = 0; r < rows; ++r) {
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
      const float obs_color = true_color + gen.normal(0, SIGMA);
      // determine the true pixel id
      const graphlab::vertex_id_type vid = sub2ind(rows,cols,r,c);
      true_pixels[vid] = true_color; obs_pixels[vid] = obs_color;
    } // end of loop over cols
  } // end of loop over rows

  // load the graph
  for(size_t r = dc.procid(); r < rows; r += dc.numprocs()) {
    for(size_t c = 0; c < cols; ++c) {
      if(r + 1 < rows) {
        vertex_data vdata;
        vdata.i = sub2ind(rows,cols,r,c);
        vdata.j = sub2ind(rows,cols,r+1,c);
        vdata.obs_color_i = obs_pixels[vdata.i];
        vdata.obs_color_j = obs_pixels[vdata.j];
        //        graph.add_vertex(
        graph.add_edge(vid, sub2ind(rows,cols,r+1,c),
                       edge_data(vid, sub2ind(rows,cols,r+1,c), NCOLORS));
      }
    }
  }



  dc.barrier();
}; // end of create synthetic mrf
