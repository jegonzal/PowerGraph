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
 *  Any changes to the code must include this original license notice in full.
 * Original code by Yucheng Low, CMU
 * Modified to Gaussian Mixture by DAnny Bickson, CMU
 * Based on Matlab code by Alex Ihler, UC Irvine
 * See the paper: Nonparametric Belief Propagation. E. Sudderth, A. Ihler, W. Freeman, and A. Willsky. CVPR, June 2003.
 */


#include <iostream>
#include <map>
#include <graphlab.hpp>
#include <limits>
#include "kde.h"
#include <float.h>
#include "image.hpp"
#include "prodSampleEpsilon.hpp"
#include <set>

#include <itpp/itstat.h>
#include <itpp/itbase.h>

#include <graphlab/macros_def.hpp>

#define PROPOSAL_STDEV 20

using namespace itpp;
using namespace std;


int NSAMP =12; //number of samples
double EPSILON =1e-5; //epsilon (accuracy of product sampling)
int MAX_ITERATIONS = 10;
int iiter = 0;

// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
  kde msg; //the NBP message sent along this edge
  kde edge_pot; //edge potential of this edge
  int update_count;
}; 

struct vertex_data: public graphlab::unsupported_serialize {
  kde obs; //ovservation
  kde bel; //belief
  int rounds;
  vertex_data(){ rounds = 0;}
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;


/**
 * compare MAE (mean average error) of the true image vs. inferred image
 */
double image_compare_mae(image &trueimg, image &infered) {
    assert(trueimg.rows() == infered.rows());
    assert(trueimg.cols() == infered.cols());
    // get the set of colors in the trueimg
    std::set<int> colors;
    for (size_t i = 0; i < trueimg.rows(); ++i) {
      for (size_t j = 0; j < trueimg.cols(); ++j) {
        colors.insert(size_t(trueimg.pixel(i,j)));
      }
    }
    
    // fill a rounding color map
    int colormap[256];
    int previval = -256;
    std::set<int>::iterator curi = colors.begin();
    std::set<int>::iterator nexti = curi;
    nexti++;
    int nextival = (nexti != colors.end())?*nexti:512;
    while (curi != colors.end()) {
      int low = (previval + (*curi)) / 2; if (low < 0) low = 0;
      int high = (nextival + (*curi)) / 2; if (high > 256) high = 256;
      
      for (int i = low; i < high; ++i) {
          colormap[i] = (*curi);
      }
      previval = (*curi);
      curi++;
      nexti++;
      nextival = (nexti != colors.end())?*nexti:512;
    }
    
    // compute absolute difference
    double err = 0;
    for (size_t i = 0; i < infered.rows(); ++i) {
      for (size_t j = 0; j < infered.cols(); ++j) {
        //err  += (infered.pixel(i,j) - trueimg.pixel(i,j)) * (infered.pixel(i,j) - trueimg.pixel(i,j)) ;
        err += fabs(infered.pixel(i,j) - trueimg.pixel(i,j));
      }
    }
    err  /= (infered.rows() * infered.cols());
    return err;
}

/**
 * compare RMSE (root mean square error) of true image vs. inferred image
 */
double image_compare_rmse(image &trueimg, image &infered) {
    assert(trueimg.rows() == infered.rows());
    assert(trueimg.cols() == infered.cols());
    // get the set of colors in the trueimg
    std::set<int> colors;
    for (size_t i = 0; i < trueimg.rows(); ++i) {
      for (size_t j = 0; j < trueimg.cols(); ++j) {
        colors.insert(size_t(trueimg.pixel(i,j)));
      }
    }
    
    // fill a rounding color map
    int colormap[256];
    int previval = -256;
    std::set<int>::iterator curi = colors.begin();
    std::set<int>::iterator nexti = curi;
    nexti++;
    int nextival = (nexti != colors.end())?*nexti:512;
    while (curi != colors.end()) {
      int low = (previval + (*curi)) / 2; if (low < 0) low = 0;
      int high = (nextival + (*curi)) / 2; if (high > 256) high = 256;
      
      for (int i = low; i < high; ++i) {
          colormap[i] = (*curi);
      }
      previval = (*curi);
      curi++;
      nexti++;
      nextival = (nexti != colors.end())?*nexti:512;
    }
    
    // compute absolute difference
    double err = 0;
    for (size_t i = 0; i < infered.rows(); ++i) {
      for (size_t j = 0; j < infered.cols(); ++j) {
        err  += (infered.pixel(i,j) - trueimg.pixel(i,j)) * (infered.pixel(i,j) - trueimg.pixel(i,j)) ;
//        err += fabs(infered.pixel(i,j) - trueimg.pixel(i,j));
      }
    }
    err  /= (infered.rows() * infered.cols());
    return err;
}


/**
 * Non-parametric BP update function
 */
void nbp_update(gl_types::iscope& scope,
               gl_types::icallback& scheduler) {


  bool debug = true;

  vertex_data& v_data = scope.vertex_data();
  graphlab::vertex_id_t vid = scope.vertex();
  if (debug && vid%1000000 == 0){
     std::cout<<"Entering node " << (int)vid << " obs: ";
     v_data.obs.matlab_print();
     std::cout << std::endl;
  }

  v_data.rounds++;

  if ((int)vid == 0)
      iiter++;

  gl_types::edge_list in_edges = scope.in_edge_ids();
  gl_types::edge_list out_edges = scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size()); // Sanity check

  //for each incoming message
  for (size_t j = 0; j < in_edges.size(); ++j){
   
     std::vector<kde> kdes;
     for(size_t i = 0; i < in_edges.size(); ++i) {
    
      graphlab::edge_id_t ineid = in_edges[i];
      edge_data& in_edge = scope.edge_data(ineid);
      if (i != j){
         in_edge.msg.verify();
         //add the message into the mixture list
         kdes.push_back(in_edge.msg);
      }
     }

     kdes.push_back(v_data.obs);
     
     graphlab::edge_id_t outeid = out_edges[j];
     edge_data& out_edge = scope.edge_data(outeid);
     kde marg = out_edge.edge_pot.marginal(0);  
      //insert the marginal of this dimension as the first item of the mixture
     kdes.insert(kdes.begin(), marg);//important: has to be first!     

      //compute the mixtures product
     prodSampleEpsilon producter; 
     kde m = producter.prodSampleEpsilonRun(kdes.size(), NSAMP, EPSILON, kdes);
     
     m.verify();
     kde mar2 = out_edge.edge_pot.marginal(1);
     mar2.verify();
     imat firstrowind = m.indices(0,0,0,m.indices.cols()-1);
      //sample from marginal, using indices taken from product
     kde outmsg = mar2.sample(firstrowind,m.weights);
     outmsg.verify(); 
     out_edge.msg = outmsg;

  }

   //compute belief by multiplying self potential with all incoming message
   if (v_data.rounds == MAX_ITERATIONS){
	if (debug && vid%100000 == 0)
	   printf("computing belief node %d\n", vid);

      std::vector<kde> kdes;
      for (size_t j = 0; j < in_edges.size(); ++j){
        graphlab::edge_id_t ineid = in_edges[j];
        edge_data& in_edge = scope.edge_data(ineid);
        in_edge.msg.verify();
        //add all incoming message to mixture list
        kdes.push_back(in_edge.msg);
      }

      //add self potential
      kdes.push_back(v_data.obs);
      //compute the product
      prodSampleEpsilon prod;
      kde m = prod.prodSampleEpsilonRun(kdes.size(), NSAMP, EPSILON, kdes);
     
      m.verify();
      //store the result
      v_data.bel = m;
      if (debug && vid == 0){
	   printf("belief node %d is\n", vid);
           m.matlab_print(); printf("\n");
      }
   }

} // end of nbp_update





void construct_graph(image& img,
                     kde & edge_pot,
                     gl_types::graph& graph) {

  
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      vertex_data vdat;
      vdat.rounds = 0;

      // Set the node potentials
      vec cent = zeros(2);
      //center of mixture component is around pixel color
      cent[0] = img.pixel(i,j);
      cent[1] = img.pixel(i,j);
      mat cent2 = cent; cent2 = transpose(cent2);
      vec bw = "30 30";
      mat bw2 = bw; bw2 = transpose(bw2);
      vec wght = "1 1";
      //create a mixture
      vdat.obs = kde(cent2, bw2, wght);
      vdat.bel = vdat.obs;
      if (i == 0 && j == 0)
       vdat.obs.matlab_print();
      graph.add_vertex(vdat);
      vdat.obs.verify();
      vdat.bel.verify();

    } // end of for j in cols
  } // end of for i in rows

  edge_data edata;
  edata.edge_pot = edge_pot;
  edata.edge_pot.matlab_print();

  //add the edges to the grid graph
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
  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

  std::string gmmfile= "";
  std::string inputfile = "";

    // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("NBP image denoising");
  clopts.attach_option("epsilon",
                       &EPSILON, EPSILON,
                       "epsilon - product accuracy");
  clopts.attach_option("gmmfile",
                       &gmmfile, std::string(""),
                       "true image + self and edge potential file");
  clopts.attach_option("inputfile",
                       &inputfile, std::string(""),
                       "the input noisy image");
  clopts.attach_option("max_iter", &MAX_ITERATIONS, MAX_ITERATIONS, "maximum number of iterations. In this round the belief is compted");

  clopts.set_scheduler_type("round_robin");

  bool success = clopts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }

  // load the potentials mixture components
  it_ifile f(gmmfile.c_str());

  mat edgecenter, edgesigma, edgeweight;
  mat nodecenter, nodesigma, nodeweight;
  ivec truedata;
  ivec imgsize;
  imat integermat;
  vec doublevec;
  //read edge potentials
  f >> Name("edge_ce") >> integermat;   edgecenter = to_mat(integermat);
  f >> Name("edge_alpha") >> doublevec; edgeweight = doublevec;
  f >> Name("edge_sigma") >> doublevec; edgesigma = doublevec;

  //read self potential
  f >> Name("like_ce") >> nodecenter;
  f >> Name("like_alpha") >> doublevec; nodeweight = doublevec;
  f >> Name("like_sigma") >> doublevec; nodesigma = doublevec;
  
  //read true image
  f >> Name("img1") >> truedata;
  //read image size
  f >> Name("isize") >> imgsize;

  size_t rows = imgsize(0);
  size_t cols = imgsize(1);
  std::cout << "Image size is " << rows << " x " << cols << std::endl;
  mat edgesigma2 = edgesigma;
  edgesigma2 = transpose(edgesigma2);
  mat edgeweight2 = edgeweight;
  edgeweight2 = transpose(edgeweight2); 
  if (edgesigma2.cols() > edgecenter.cols())
	edgesigma2 = edgesigma2(0,0,0,edgecenter.cols()-1);
  kde edge_pot = kde(edgecenter, edgesigma2, edgeweight2);

// convert the true image to an image
  image trueimg(rows, cols);
  for (size_t i = 0; i < size_t(truedata.size()); ++i) {
    trueimg.pixel(i) = truedata(i);
  }

  //read noisy image
  it_ifile imgfile(inputfile.c_str());
  vec observations;
  imgfile >> Name("obs2") >> observations;
  // convert observations to an image
  image img(rows, cols);
  for (size_t i = 0;i < size_t(observations.size()); ++i) {
    img.pixel(i) = observations(i);
  }
  img.save("noisy.pgm");
  trueimg.save("source_img.pgm");

  // Create the graph --------------------------------------------------------->
  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);

  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(img, edge_pot, core.graph());

  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function", nbp_update);
  std::cout << "Running the engine. " << std::endl;

  // Add the bp update to all vertices
  core.add_task_to_all(nbp_update, 100.0);
  // Starte the engine
  const double runtime = core.start();

  // Saving the output -------------------------------------------------------->
  std::cout << "Rendering the cleaned image. " << std::endl;
  
  //parse belief to find the reconstructed image 
  for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    float a = vdata.bel.max();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = size_t(a);
  }
  double err = sqrt(image_compare_rmse(trueimg, img));
  double err2 = image_compare_mae(trueimg, img);
  img.save("inferred.pgm");
  std::cout << "RMSE: " << err << " MAE: "<< err2<<std::endl;
  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
} // End of main

