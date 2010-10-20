/*
 * kbp.cpp
 *
 * This file contains codes for kernel
 * belief propagation in a pairwise markov random field to
 * estimate depth from features from still images;
 *
 *  Created on: Oct 10, 2010
 *      Author: lesong
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <iostream>
#include <map>

#include <graphlab.hpp>
#include "kde.h"

// include itpp
#include <itpp/itstat.h>
#include <itpp/itbase.h>
#include "../kernelbp/image_compare.hpp"
#include "prodSampleEpsilon.hpp"

#include <fstream>
#include <cmath>
#include <cstdio>
#include <cfloat>
//#include <graphlab/graph/graph.hpp>

#include <graphlab/schedulers/round_robin_scheduler.hpp>
#include <graphlab/macros_def.hpp>

#define NSAMP 24
#define EPSILON 1e-5
int RESAMPLE_FREQUENCY = 0;
int MAX_ITERATIONS = 6;

using namespace itpp;
using namespace std;

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

mat* pUud;
mat* pUdu;
mat* pUlr;
mat* pUrl;
mat ipUud;
mat ipUdu;
mat ipUlr;
mat ipUrl;

ivec isize;
vec m0;
vec predy_k;
mat predfy_k;
mat prod_msg0;
mat belief0;
vec testy;
vec truey;
mat cfybasis;
mat ftesty;

size_t msgdim;
size_t basisno;
size_t testno;

int iiter = 0;

double damping = 0.90;

// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */

/**
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.
 */
void bp_update(gl_types::iscope& scope,
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data);



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
  if (debug && vid == 0){
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
      kdes.push_back(marg);      

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
   if (v_data.rounds == MAX_ITERATIONS - 1){
	if (debug && vid == 0)
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
	   printf("computing belief node %d\n", vid);
           m.matlab_print(); printf("\n");
      }


   }
  /*
  if (v_data.rounds < MAX_ITERATIONS) {
    gl_types::update_task task(scope.vertex(), bp_update);
    scheduler.add_task(task, 1.0);
  }*/


} // end of BP_update



void construct_graph(std::string gmmfile,
                    double numparticles,
                     gl_types::graph& graph,
                     size_t rows,
                     size_t cols) {


  image img(rows, cols);
  // initialize a bunch of particles
  for(size_t i = 0; i < rows; ++i) {
    
    it_ifile gmms((gmmfile+"_part"+boost::lexical_cast<std::string>(i+1)+".it").c_str());
    
    mat nodeweight;
    vec doublevec;
    gmms >> Name("like_ce") >> doublevec; mat nodecenter(1,doublevec.size()); nodecenter = doublevec; 
    gmms >> Name("like_alpha") >> nodeweight;
    gmms >> Name("like_sigma") >> doublevec; 
    mat nodesigma(1,doublevec.size()); nodesigma = doublevec;
    
    ASSERT_EQ(nodeweight.rows(), cols);
    
    for(size_t j = 0; j < cols; ++j) {
      vertex_data vdat;
      //vdat.rounds = 0;

      // Set the node potential
      itpp::vec weights = nodeweight.get_row(j);
      if (nodecenter.rows() > nodecenter.cols())
      	nodecenter = itpp::transpose(nodecenter); 
      if (nodesigma.rows() > nodesigma.cols())
         nodesigma = itpp::transpose(nodesigma);
      vdat.obs = kde(nodecenter, nodesigma, weights);
      vdat.obs.verify();
      //vdat.p.simplify();
    
      //vdat.bel = vdat.obs.sample();//TOOD
      vdat.bel = vdat.obs;
      vdat.bel.verify();
      /*for (size_t n = 0;n < numparticles; ++n) {
        particle p;
        p.x = vdat.p.sample();
        p.weight = 1.0/numparticles;
        vdat.belief.push_back(p);
      }

      if(i == 0 && j == 0) {
        std::cout << "vertex 0 particles: ";
        for (size_t i = 0;i < vdat.belief.size(); ++i) {
          std::cout  << vdat.belief[i].x << " ";
        }
        std::cout << "\n";
      }*/
      graph.add_vertex(vdat);

    } // end of for j in cols
  } // end of for i in rows



  // Add the edges
  
  edge_data edata;

  //GaussianMixture<2> lrpot;
  //GaussianMixture<2> udpot;
  //GaussianMixture<2> lastudpot;
  kde lrpot, udpot, lastudpot;
  //
  for(size_t i = 0; i < rows; ++i) {
    
    it_ifile gmms(gmmfile+"_part"+boost::lexical_cast<std::string>(i+1)+".it");
    
    mat lrcenter;
    vec doublevec;
    gmms >> Name("lr_edge_ce") >> lrcenter;
    gmms >> Name("lr_edge_alpha") >> doublevec; 
    mat lrweight(1,doublevec.size());  
    lrweight = doublevec;
    gmms >> Name("lr_edge_sigma") >> doublevec; 
    mat lrsigma(1,doublevec.size());
    lrsigma = doublevec;

    mat udcenter;
    gmms >> Name("ud_edge_ce") >> udcenter;
    gmms >> Name("ud_edge_alpha") >> doublevec; 
    mat udweight(1,doublevec.size());
    udweight = doublevec;
    gmms >> Name("ud_edge_sigma") >> doublevec; 
    mat udsigma(1,doublevec.size());
    udsigma = doublevec;
    lastudpot = udpot;
    lrsigma = transpose(lrsigma); 
    udsigma = transpose(udsigma); 
    lrweight = transpose(lrweight);
    udweight = transpose(udweight);

    //lrpot = GaussianMixture<2>(lrcenter, lrsigma, lrweight);
    lrpot = kde(lrcenter, lrsigma, lrweight);
    lrpot.verify();
    //udpot = GaussianMixture<2>(udcenter, udsigma, udweight);
    udpot = kde(udcenter, udsigma, udweight);   
    udpot.verify();

       for(size_t j = 0; j < cols; ++j) {

      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        edata.msg = graph.vertex_data(img.vertid(i-1, j)).bel;
        edata.edge_pot = lastudpot;
        graph.add_edge(vertid, img.vertid(i-1, j), edata);
      }
      if(i+1 < img.rows()) {
        edata.msg = graph.vertex_data(img.vertid(i+1, j)).bel;
        edata.edge_pot = udpot;
        graph.add_edge(vertid, img.vertid(i+1, j), edata);
      }
      if(j-1 < img.cols()) {
        edata.msg = graph.vertex_data(img.vertid(i, j-1)).bel;
        edata.edge_pot = lrpot;
        graph.add_edge(vertid, img.vertid(i, j-1), edata);
      }
      if(j+1 < img.cols()) {
        edata.msg = graph.vertex_data(img.vertid(i, j+1)).bel;
        edata.edge_pot = lrpot;
        graph.add_edge(vertid, img.vertid(i, j+1), edata);
      }
    } // end of for j in cols
  } // end of for i in rows
  graph.finalize();
} // End of construct graph



void test(){
  test_marginal();
  test_max();
  test_sample();
  test_sample2();
  test_ROT();
  test_product();
  exit(0);
}


// MAIN =======================================================================>
int main(int argc, char** argv) {
  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

  //test();

  size_t iterations = 100;
  size_t numparticles = 100;
  std::string pred_type = "map";

  std::string gmmfile= "";
  std::string kbpfile = "";
  std::string logfile = "";

  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Loopy BP image denoising");
  clopts.attach_option("iterations",
                       &iterations, iterations,
                       "Number of iterations");
  clopts.attach_option("particles",
                       &numparticles, numparticles,
                       "Number of particlesw");
  clopts.attach_option("gmmfile",
                       &gmmfile, std::string(""),
                       "gmm mixture file");
  clopts.attach_option("kbpfile",
                       &kbpfile, std::string(""),
                       "kbpfile");
  clopts.attach_option("damping",
                       &damping, damping,
                       "damping");
  clopts.attach_option("resample",
                       &RESAMPLE_FREQUENCY, NSAMP,
                       "resampling frequency");
  //clopts.attach_option("mcmc",
  //                     &MCMCSTEPS, MCMCSTEPS,
  //                     "mcmc steps per resample");
  clopts.attach_option("logfile",
                       &logfile, std::string(""),
                       "log file");

  // set default scheduler type
  clopts.scheduler_type = "round_robin";
  clopts.scope_type = "none";

  bool success = clopts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }

  // fill the global vars
  MAX_ITERATIONS = iterations;

  // load ground truth
  vec testy;
  vec truey;
  it_ifile kbppart0((kbpfile+"_part"+boost::lexical_cast<std::string>(0)+".it").c_str());
  kbppart0 >> Name("testy") >> testy;
  kbppart0 >> Name("truey") >> truey;

  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);
  
  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(gmmfile, numparticles, core.graph(),107,86);
  


  // Running the engine ------------------------------------------------------->
  core.scheduler().set_option(gl_types::scheduler_options::UPDATE_FUNCTION,
                              (void*)bp_update);
  core.scheduler().set_option(gl_types::scheduler_options::MAX_ITERATIONS,
                              (void*)&MAX_ITERATIONS);

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

  vec pred(107*86);
  image img(107, 86);
  image trueimg(107, 86);
  image transposedimg(86, 107);
  for (size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    // ok... this is unbelievably annoying but my traversal order is opposite
    // of the data's traversal order
    // real physical location
    std::pair<size_t, size_t> loc = img.loc(v);
    // index in transposed image
    size_t v2 = transposedimg.vertid(loc.second, loc.first);
    pred(v2) = vdata.bel.max();
    img.pixel(v) = pred(v2);
    trueimg.pixel(v) = truey(v2);
  }
  img.save("pred.pgm", false, 0, 2);
  trueimg.save("true.pgm", false, 0, 2);
  double mae = mean(abs(pred - truey));
  std::cout << "Mean absolute error: " << mae << std::endl;
  
  if (logfile.length() != 0) {
    ofstream fout;
    fout.open(logfile.c_str(), ios::app);
    gmmfile = gmmfile.rfind("/") == std::string::npos ?
                                       gmmfile:
                                       gmmfile.substr(gmmfile.rfind("/")+1);
    
    fout << "pbp\t" << gmmfile << "\t" << numparticles << "\t" << RESAMPLE_FREQUENCY << "\t" << mae << "\t"
    << runtime << std::endl;
    fout.close();
  }

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
} // End of main

