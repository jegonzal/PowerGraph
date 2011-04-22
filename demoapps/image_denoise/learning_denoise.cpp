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
}; // End of edge data


/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data {
  graphlab::unary_factor potential;
  graphlab::unary_factor belief;
}; // End of vertex data



typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;

gl_types::glshared<double> sh_bound;
gl_types::glshared<double> sh_damping;
gl_types::glshared<double> sh_ipfdamping;
gl_types::glshared<graphlab::binary_factor> sh_edgepot;
gl_types::glshared<graphlab::binary_factor> sh_truecounts;

// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     gl_types::graph& graph);

/** Get the counts of the true image */
void get_image_counts(image &trueimg, 
                      graphlab::binary_factor &truebinarycount, 
                      size_t arity);

/** Tests if two binary factors are equal */
bool binary_factor_equal(const graphlab::binary_factor &a, 
                         const graphlab::binary_factor &b);
/** 
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.  
 */
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler);
               
/**
 * The edge potential sync simply counts the belief across the edges
 */
void edgepot_sync(gl_types::iscope &scope,  graphlab::any& acc);

/**
 * Performs an IPF update using the counts and the true counts
 */
void edgepot_apply(graphlab::any& result,  const graphlab::any& acc);

void edgepot_merge(graphlab::any& result,  const graphlab::any& acc);

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

void get_image_counts(image &trueimg, 
                      graphlab::binary_factor &truebinarycount, 
                      size_t arity) {
  truebinarycount.resize(arity,arity);
  for (size_t i = 0; i < trueimg.rows(); ++i) {
    for (size_t j = 0; j < trueimg.cols(); ++j) {
      size_t color = trueimg.pixel(i,j);
      // check neighbors
      if (i > 1) {
        size_t color2 = trueimg.pixel(i-1,j); 
        truebinarycount.logP(color,color2)++;
      }
      if (i < trueimg.rows() - 1) {
        size_t color2 = trueimg.pixel(i+1,j); 
        truebinarycount.logP(color,color2)++;
      }
      if (j > 1) {
        size_t color2 = trueimg.pixel(i,j-1); 
        truebinarycount.logP(color,color2)++;
      }
      if (j < trueimg.cols() - 1) {
        size_t color2 = trueimg.pixel(i,j+1); 
        truebinarycount.logP(color,color2)++;
      }
    }
  }
}

bool binary_factor_equal(const graphlab::binary_factor &a, 
                         const graphlab::binary_factor &b) {
  ASSERT_EQ(a.arity1(), b.arity1());
  ASSERT_EQ(a.arity2(), b.arity2());
  for (size_t i = 0;i < a.arity1(); ++i) {
     for (size_t j = 0;j < a.arity2(); ++j) {
      if ((a.logP(i,j) - b.logP(i,j)) > std::numeric_limits<double>::epsilon()) {
        return false;
      }
    }
  }
  return true;
}
// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using loopy belief propagation inside " << std::endl
            << "the graphlab framework." << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  double bound = 1E-4;
  double damping = 0.1;
  double ipfdamping = 0.1;
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
  clopts.attach_option("bound",
                       &bound, bound,
                       "Residual termination bound");
  clopts.attach_option("damping",
                       &damping, damping,
                       "The amount of message damping (higher = more damping)");
  clopts.attach_option("ipfdamping",
                       &ipfdamping, ipfdamping,
                       "The amount of IPF damping. (lower == more damping)");
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
  
  

  clopts.set_scheduler_type("splash(splash_size=100)");
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

  
  

  // Create synthetic images -------------------------------------------------->
  // Creating image for denoising
  std::cout << "Creating a synethic image. " << std::endl;
  image img(rows, cols);
  img.paint_sunset(colors);
  graphlab::binary_factor truecounts;
  get_image_counts(img, truecounts, colors);
  std::cout << "Saving image. " << std::endl;
  img.save(orig_fn.c_str());
  std::cout << "Corrupting Image. " << std::endl;
  img.corrupt(sigma);
  std::cout << "Saving corrupted image. " << std::endl;
  img.save(noisy_fn.c_str());


  std::cout << "True Counts: " << std::endl;
  std::cout << truecounts;
 
  
  
  // Create the graph --------------------------------------------------------->
  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);
  
  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(img, colors, sigma, core.graph());

  
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
  
  
  // fill the shared variales
  sh_bound.set(bound);
  sh_damping.set(damping);
  sh_truecounts.set(truecounts);
  sh_edgepot.set(edge_potential);
  sh_ipfdamping.set(ipfdamping);
  // make the true counts
  

  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function",bp_update);
  std::cout << "Running the engine. " << std::endl;

  graphlab::binary_factor zero;
  zero.resize(colors,colors);
  zero.set_as_agreement(0);
  core.set_sync(sh_edgepot,
                edgepot_sync,
                edgepot_apply,
                zero,
                rows*cols,
                edgepot_merge);
  graphlab::binary_factor oldedgepot = sh_edgepot.get_val();
  graphlab::timer ti;
  ti.start();
  // loop it a few times
  for (size_t i = 0;i < 10; ++i) {
    std::cout << "restart " << i << "\n";
    // Add the bp update to all vertices
    core.add_task_to_all(bp_update, 100.0);
    // Start the engine
    core.start();
    graphlab::binary_factor newedgepot = sh_edgepot.get_val();
    if (binary_factor_equal(oldedgepot, newedgepot)) break;
    oldedgepot = newedgepot;
  }
  double runtime = ti.current_time();
  size_t update_count = core.last_update_count();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;  


  // Saving the output ----------------------------------------------------->
  std::cout << "Rendering the cleaned image. " << std::endl;
  if(pred_type == "map") {
    for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
      const vertex_data& vdata = core.graph().vertex_data(v);
      img.pixel(v) = vdata.belief.max_asg();    
    }
  } else if(pred_type == "exp") {
    for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
      const vertex_data& vdata = core.graph().vertex_data(v);
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
  return EXIT_SUCCESS;
} // End of main




// Implementations
// ============================================================>
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler) {
  //  std::cout << scope.vertex();;
  //  std::getchar();
  double bound = sh_bound.get_val();
  double damping = sh_damping.get_val();

  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  
  // Get the in and out edges by reference
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();
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
  boost::shared_ptr<const graphlab::binary_factor> edge_factor_ptr = 
    sh_edgepot.get_ptr();
  
  
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
    tmp_msg.convolve(*edge_factor_ptr, cavity);
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


void edgepot_sync(gl_types::iscope &scope,  graphlab::any& acc) {
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size()); // Sanity check
  
  graphlab::binary_factor& counts = acc.as<graphlab::binary_factor>();
  
  // Get the in and out edge data
  // the edge belief of u -- v
  // belief of u / msg_{v->u) * belief of v / msg_{u->v} * edgepot;

  const graphlab::unary_factor& blfu = scope.const_vertex_data().belief;
  
  foreach(graphlab::edge_id_t ineid, in_edges) {   
    
    graphlab::vertex_id_t srcv = scope.source(ineid);
    // message from v->u
    const graphlab::unary_factor &msgvu = 
      scope.const_edge_data(ineid).message;
    // belief at v
    const graphlab::unary_factor& blfv= 
      scope.const_neighbor_vertex_data(srcv).belief;
    // get the message from u->v. requires the reverse edge
    graphlab::edge_id_t outeid = scope.reverse_edge(ineid);
    const graphlab::unary_factor &msguv = 
      scope.const_edge_data(outeid).message;
    
    boost::shared_ptr<const graphlab::binary_factor> edge_factor = 
      sh_edgepot.get_ptr();
    
    graphlab::binary_factor edge_belief;
    edge_belief.resize(blfu.arity(), blfv.arity());
    // loop through my assignments and my neighbor assignments
    
    // using logP to store actual counts
    for (size_t i = 0;i < blfu.arity(); ++i) {
      for (size_t j = 0;j < blfv.arity(); ++j) {
        edge_belief.logP(i,j) = 
          blfu.logP(i) - msgvu.logP(i) + blfv.logP(j) - 
          msguv.logP(j) + edge_factor->logP(i,j);
      }
    }
    edge_belief.normalize();
    for (size_t i = 0;i < edge_belief.arity1(); ++i) {
      for (size_t j = 0;j < edge_belief.arity2(); ++j) {
        counts.logP(i,j) += std::exp(edge_belief.logP(i,j));
      }
    }
    
  }
}

void edgepot_merge(graphlab::any& result,  const graphlab::any& acc) {
  graphlab::binary_factor& res = result.as<graphlab::binary_factor>();
  const graphlab::binary_factor& a = acc.as<graphlab::binary_factor>();

  for (size_t i = 0;i < res.arity1(); ++i) {
    for (size_t j = 0;j < res.arity2(); ++j) {
      res.logP(i,j) += a.logP(i,j);
    }
  }
}

void edgepot_apply(graphlab::any& result,  const graphlab::any& acc) {
  // IPF update
  graphlab::binary_factor& res = result.as<graphlab::binary_factor>();
  graphlab::binary_factor truecounts = sh_truecounts.get_val();
  const graphlab::binary_factor& curcounts = 
    acc.as<graphlab::binary_factor>();
  double ipfdamping = sh_ipfdamping.get_val();
  // perform the IPF update
  // note that BP+IPF can be quite unstable
  // so lets only update the parameter values if they change by > 1E-2
  for (size_t i = 0;i < res.arity1(); ++i) {
    for (size_t j = 0;j < res.arity2(); ++j) {
      // + 100 to avoid divide by 0 problems
      double newval = 
        ipfdamping * log((truecounts.logP(i,j)+100) / 
                         (curcounts.logP(i,j)+100)) + 
        (1 - ipfdamping) * res.logP(i,j);
      if (std::fabs(res.logP(i,j) - newval) >= 1E-2) {
        res.logP(i,j) = newval;
      }
    }
  }
  std::cout << "sync of edge pot!\n";
  std::cout << res << std::endl;
}


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


