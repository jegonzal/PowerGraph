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

double BOUND;
double DAMPING;
double IPF_DAMPING;

// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     graph_type& graph);


class bp_update : 
  public graphlab::iupdate_functor<graph_type, bp_update> {
 
  typedef graphlab::iupdate_functor<graph_type, bp_update> base;
  typedef base::icontext_type      icontext_type;
  typedef base::vertex_id_type   vertex_id_type;
  typedef base::edge_type     edge_type;
  typedef base::edge_list_type   edge_list_type;

private:
  double residual;

public:

  bp_update(double residual = 1000) : residual(residual) { }

  double priority() const { return residual; }
  void operator+=(const bp_update& other) { residual += other.residual; }

  void operator()(icontext_type& context) {

    // Grab the state from the context
    // ---------------------------------------------------------------->
    // Get the vertex data
    vertex_data& vdata = context.vertex_data();
  
    // Get the in and out edges by reference
    const edge_list_type in_edges  = context.in_edges();
    const edge_list_type out_edges = context.out_edges();
    assert(in_edges.size() == out_edges.size()); // Sanity check

    // Compute the belief
    // ---------------------------------------------------------------->
    vdata.belief = vdata.potential;
    foreach(const edge_type& ineid, in_edges) {   
      // Receive the message
      edge_data& edata = context.edge_data(ineid);
      edata.old_message = edata.message;
      vdata.belief.times( edata.old_message );
    }
    vdata.belief.normalize(); // finally normalize the belief
  
    // Compute outbound messages
    // ---------------------------------------------------------------->
    graphlab::binary_factor edge_factor = 
      context.get_global<graphlab::binary_factor>("EDGE_FACTOR");

    // Use cached memory spaces for cavity and tmp msg
    if(!context.is_local("cavity")) {
      context.add_local("cavity", graphlab::unary_factor());
      context.add_local("tmp_msg", graphlab::unary_factor());
    }
    graphlab::unary_factor& cavity = 
      context.get_local<graphlab::unary_factor>("cavity");
    graphlab::unary_factor& tmp_msg = 
      context.get_local<graphlab::unary_factor>("tmp_msg");
    // Send outbound messages
    for(size_t i = 0; i < in_edges.size(); ++i) {
      // Get the edge ids
      const edge_type out = out_edges[i];
      const edge_type in = in_edges[i];
      // CLEVER HACK: Here we are expoiting the sorting of the edge ids
      // to do fast O(1) time edge reversal
      ASSERT_EQ(out.target() , in.source());
      // Get the in and out edge data
      const edge_data& in_edge = context.edge_data(in);
      edge_data& out_edge = context.edge_data(out);
   
      // Compute cavity
      cavity = vdata.belief;
      cavity.divide(in_edge.old_message); // Make the cavity a cavity
      cavity.normalize();


      // convolve cavity with the edge factor storing the result in the
      // temporary message
      tmp_msg.resize(out_edge.message.arity());
      tmp_msg.var() = out_edge.message.var();
      tmp_msg.convolve(edge_factor, cavity);
      tmp_msg.normalize();
      
      // Damp the message
      tmp_msg.damp(out_edge.message, DAMPING);
    
      // Compute message residual
      double residual = tmp_msg.residual(out_edge.old_message);
    
      // Assign the out message
      out_edge.message = tmp_msg;
      
      if(residual > BOUND) {
        context.schedule(out.target(), bp_update(residual));
      }    
    }
  } // end of operator()
}; // end of class bp_update

void get_image_counts(image& trueimg, 
                      graphlab::binary_factor& truebinarycount, 
                      size_t arity) {
  truebinarycount.resize(arity,arity);
  for (size_t i = 0; i < trueimg.rows(); ++i) {
    for (size_t j = 0; j < trueimg.cols(); ++j) {
      size_t color = trueimg.pixel(i,j);
      // check neighbors
      if (i > 1) {
        const size_t color2 = trueimg.pixel(i-1,j); 
        truebinarycount.logP(color,color2)++;
      }
      if (i < trueimg.rows() - 1) {
        const size_t color2 = trueimg.pixel(i+1,j); 
        truebinarycount.logP(color,color2)++;
      }
      if (j > 1) {
        const size_t color2 = trueimg.pixel(i,j-1); 
        truebinarycount.logP(color,color2)++;
      }
      if (j < trueimg.cols() - 1) {
        const size_t color2 = trueimg.pixel(i,j+1); 
        truebinarycount.logP(color,color2)++;
      }
    }
  }
} // end of get image counts

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
} // end of factor equal

class accumulator :
  public graphlab::iaccumulator<graph_type, bp_update, accumulator> {
private:
  graphlab::binary_factor counts;
public:
  typedef graphlab::iaccumulator<graph_type, bp_update, accumulator> base;
  accumulator(size_t arity) : counts(arity, arity) { counts.set_as_agreement(0); }
  void operator()(base::icontext_type& context) {
    // get the current edge factor
    graphlab::binary_factor edge_factor = 
      context.get_global<graphlab::binary_factor>("EDGE_FACTOR");
    const graph_type::edge_list_type in_edges = context.in_edges();
    const graph_type::edge_list_type out_edges = context.out_edges();
    ASSERT_EQ(in_edges.size(), out_edges.size()); // Sanity check
    // Get the in and out edge data
    // the edge belief of u -- v
    // belief of u / msg_{v->u) * belief of v / msg_{u->v} * edgepot;
    const graphlab::unary_factor& blfu = context.const_vertex_data().belief;
    foreach(const edge_type& ineid, in_edges) {   
      const graph_type::vertex_id_type srcv = ineid.source();
      // message from v->u
      const graphlab::unary_factor &msgvu = 
        context.const_edge_data(ineid).message;
      // belief at v
      const graphlab::unary_factor& blfv = 
        context.const_vertex_data(srcv).belief;
      // get the message from u->v. requires the reverse edge
      const graph_type::edge_type outeid = context.reverse_edge(ineid);
      const graphlab::unary_factor &msguv = 
        context.const_edge_data(outeid).message;
      // compute the edge belief
      graphlab::binary_factor edge_belief;
      edge_belief.resize(blfu.arity(), blfv.arity());
      // loop through my assignments and my neighbor assignments
      // using logP to store actual counts
      for (size_t i = 0;i < blfu.arity(); ++i) {
        for (size_t j = 0;j < blfv.arity(); ++j) {
          edge_belief.logP(i,j) = 
            blfu.logP(i) - msgvu.logP(i) + blfv.logP(j) - 
            msguv.logP(j) + edge_factor.logP(i,j);
        }
      }
      edge_belief.normalize();
      for (size_t i = 0;i < edge_belief.arity1(); ++i) {
        for (size_t j = 0;j < edge_belief.arity2(); ++j) {
          counts.logP(i,j) += std::exp(edge_belief.logP(i,j));
        }
      }
    } // end of loop over in edges
  } // end of operator()
  void operator+=(const accumulator& other) { 
    for (size_t i = 0;i < counts.arity1(); ++i) {
      for (size_t j = 0;j < counts.arity2(); ++j) {
        counts.logP(i,j) += other.counts.logP(i,j);
      }
    }
  }
  void finalize(base::iglobal_context_type& context) {
    // get the current edge factor
    graphlab::binary_factor edge_factor =
      context.get_global<graphlab::binary_factor>("EDGE_FACTOR");
    graphlab::binary_factor true_counts =
      context.get_global<graphlab::binary_factor>("TRUE_COUNTS");
    // perform the IPF update
    // note that BP+IPF can be quite unstable.
    // (We do recommend the gradient update in practice)
    // so lets only update the parameter values if they change by > 1E-1
    for (size_t i = 0; i < counts.arity1(); ++i) {
      for (size_t j = 0; j < counts.arity2(); ++j) {
        // + 100 to avoid divide by 0 problems
        double newval = 
          IPF_DAMPING * log((true_counts.logP(i,j)+100) / (counts.logP(i,j)+100)) + 
          (1 - IPF_DAMPING) * edge_factor.logP(i,j);
        if (std::fabs(edge_factor.logP(i,j) - newval) >= 1E-1) {
          edge_factor.logP(i,j) = newval;
        }
      }
    }
    context.set_global("EDGE_FACTOR", edge_factor);
    std::cout << "sync of edge pot!\n";
    std::cout << edge_factor << std::endl;
  } // end of finalize
}; // end of accumulator



// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "This program creates and denoises a synthetic " << std::endl
            << "image using loopy belief propagation inside " << std::endl
            << "the graphlab framework." << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  BOUND = 1E-4;
  DAMPING = 0.1;
  IPF_DAMPING = 0.1;
  size_t colors = 4;
  size_t rows = 100;
  size_t cols = 100;
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
                       &BOUND, BOUND,
                       "Residual termination bound");
  clopts.attach_option("damping",
                       &DAMPING, DAMPING,
                       "The amount of message damping (higher = more damping)");
  clopts.attach_option("ipfdamping",
                       &IPF_DAMPING, IPF_DAMPING,
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
            << "bound:          " << BOUND << std::endl
            << "damping:        " << DAMPING << std::endl
            << "ipfdamping:     " << IPF_DAMPING << std::endl
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
  std::cout << "Creating a synthetic image. " << std::endl;
  image img(rows, cols);
  img.paint_sunset(colors);
  graphlab::binary_factor true_counts;
  get_image_counts(img, true_counts, colors);
  std::cout << "Saving image. " << std::endl;
  img.save(orig_fn.c_str());
  std::cout << "Corrupting Image. " << std::endl;
  img.corrupt(sigma);
  std::cout << "Saving corrupted image. " << std::endl;
  img.save(noisy_fn.c_str());


  std::cout << "True Counts: " << std::endl;
  std::cout << true_counts;
 
  
  
  // Create the graph --------------------------------------------------------->
  graphlab::core<graph_type, bp_update> core;
  // Set the engine options
  core.set_options(clopts);
  
  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(img, colors, sigma, core.graph());

  
  // Setup global shared variables -------------------------------------------->
  // Initialize the edge agreement factor 
  std::cout << "Initializing shared edge agreement factor. " << std::endl;

  // dummy variables 0 and 1 and num_rings by num_rings
  graphlab::binary_factor edge_factor(0, colors, 0, colors);
  // Set the smoothing type
  if(smoothing == "square") {
    edge_factor.set_as_agreement(lambda);
  } else if (smoothing == "laplace") {
    edge_factor.set_as_laplace(lambda);
  } else {
    std::cout << "Invalid smoothing stype!" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << edge_factor << std::endl;
  
  core.add_global("EDGE_FACTOR", edge_factor);
  core.add_global("TURE_COUNTS", true_counts);
  core.add_sync("ipfupdate", accumulator(colors), rows*cols);
  

  // Running the engine ------------------------------------------------------->
  std::cout << "Running the engine. " << std::endl;

 
  graphlab::binary_factor oldedgepot = 
    core.get_global<graphlab::binary_factor>("EDGE_FACTOR");
  graphlab::timer ti;
  ti.start();
  // loop it a few times
  size_t update_count = 0;
  for (size_t i = 0;i < 10; ++i) {
    std::cout << "restart " << i << "\n";
    // Add the bp update to all vertices
    core.schedule_all(bp_update(100.0));
    // Start the engine
    core.start();
    update_count += core.last_update_count();
    graphlab::binary_factor newedgepot =
      core.get_global<graphlab::binary_factor>("EDGE_FACTOR");
    std::cout << newedgepot;
    if (binary_factor_equal(oldedgepot, newedgepot)) break;
    oldedgepot = newedgepot;
  }
  
  double runtime = ti.current_time();
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



void construct_graph(image& img,
                     size_t num_rings,
                     double sigma,
                     graph_type& graph) {
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


