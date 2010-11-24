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
#include <float.h>
#include "image.hpp"

// include itpp
#include <itpp/itstat.h>
#include <itpp/itbase.h>

#include "gaussian_mixture.hpp"
#include "image_compare.hpp"
// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

#define PROPOSAL_STDEV 20

using namespace itpp;
using namespace std;



struct particle {
  float x;
  float weight;
};

float damping = 0.2;
size_t MCMCSTEPS = 10;
size_t RESAMPLE_FREQUENCY = 1000;
size_t MAX_ITERATIONS;
graphlab::atomic<size_t> proposal_count;
graphlab::atomic<size_t> accept_count;


// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
  graphlab::unary_factor message;
  graphlab::unary_factor old_message;
  graphlab::binary_factor edgepot;
}; // End of edge data

struct vertex_data: public graphlab::unsupported_serialize {
  // the belief at this vertex. Also provides the shape of all the incoming messages
  graphlab::unary_factor potential;
  graphlab::unary_factor belief;
  float obs;
  uint rounds;
};



typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;



// Implementations
// ============================================================>
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data) {
  //  std::cout << scope.vertex();;
  //  std::getchar();

  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  if (v_data.rounds >= MAX_ITERATIONS) return;
  if (scope.vertex() == 0) {
    std::cout << "Vertex 0 Round " << v_data.rounds << std::endl;
  }
  v_data.rounds++;
  // Get the in and out edges by reference
  graphlab::edge_list in_edges = 
    scope.in_edge_ids();
  graphlab::edge_list out_edges = 
    scope.out_edge_ids();
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
    tmp_msg.convolve(out_edge.edgepot, cavity);
    tmp_msg.normalize();

    // Damp the message
    tmp_msg.damp(out_edge.message, damping);
    
    
    // Assign the out message
    out_edge.message = tmp_msg;
  }
  if (v_data.rounds < MAX_ITERATIONS) {
    gl_types::update_task task(scope.vertex(), bp_update);
    scheduler.add_task(task, 1.0);
  }

} // end of BP_update
















void construct_graph(std::string gmmfile,
                    double numparticles,
                     gl_types::graph& graph,
                     size_t rows,
                     size_t cols,
                     std::vector<double> discretizationpts) {
  size_t numasgs = discretizationpts.size();
  image img(rows, cols);
  // initialize a bunch of particles
  for(size_t i = 0; i < rows; ++i) {
    
    it_ifile gmms((gmmfile+"_part"+boost::lexical_cast<std::string>(i+1)+".it").c_str());
    
    mat nodecenter, nodesigma, nodeweight;
    vec doublevec;
    gmms >> Name("like_ce") >> doublevec; nodecenter = doublevec;
    gmms >> Name("like_alpha") >> nodeweight;
    gmms >> Name("like_sigma") >> doublevec; nodesigma = doublevec;
    
    ASSERT_EQ(nodeweight.rows(), cols);
    
    for(size_t j = 0; j < cols; ++j) {
      vertex_data vdat;
      vdat.rounds = 0;

      // Set the node potential
      GaussianMixture<1> vpot(nodecenter, nodesigma, mat(nodeweight.get_row(j)));
      vdat.potential.resize(numasgs);
      for(size_t pred = 0; pred < numasgs; ++pred) {
        vdat.potential.logP(pred) = log(vpot.likelihood(discretizationpts[pred]));
      }
      vdat.potential.normalize();
      vdat.belief = vdat.potential;
      graph.add_vertex(vdat);

    } // end of for j in cols
  } // end of for i in rows



  // Add the edges
  
  edge_data edata;

  GaussianMixture<2> lrpot;
  GaussianMixture<2> udpot;
  GaussianMixture<2> lastudpot;
  graphlab::binary_factor lrpotfactor, udpotfactor, lastudpotfactor;
  for(size_t i = 0; i < rows; ++i) {
    
    it_ifile gmms(gmmfile+"_part"+boost::lexical_cast<std::string>(i+1)+".it");
    
    mat lrcenter, lrsigma, lrweight;
    vec doublevec;
    gmms >> Name("lr_edge_ce") >> lrcenter;
    gmms >> Name("lr_edge_alpha") >> doublevec; lrweight = doublevec;
    gmms >> Name("lr_edge_sigma") >> doublevec; lrsigma = doublevec;

    mat udcenter, udsigma, udweight;
    gmms >> Name("ud_edge_ce") >> udcenter;
    gmms >> Name("ud_edge_alpha") >> doublevec; udweight = doublevec;
    gmms >> Name("ud_edge_sigma") >> doublevec; udsigma = doublevec;
    lastudpot = udpot;
    lrpot = GaussianMixture<2>(lrcenter, lrsigma, lrweight);
    udpot = GaussianMixture<2>(udcenter, udsigma, udweight);
    
    lastudpotfactor = udpotfactor;
    lrpotfactor.resize(numasgs, numasgs);
    for (size_t p = 0; p < numasgs; ++p) {
      for (size_t q = 0; q < numasgs; ++q) {
        double d[2]; d[0] = discretizationpts[p]; d[1] = discretizationpts[q];
        lrpotfactor.logP(p,q) = log(lrpot.likelihood(d));
      }
    }
    
    udpotfactor.resize(numasgs, numasgs);
    for (size_t p = 0; p < numasgs; ++p) {
      for (size_t q = 0; q < numasgs; ++q) {
        double d[2]; d[0] = discretizationpts[p]; d[1] = discretizationpts[q];
        udpotfactor.logP(p,q) = log(udpot.likelihood(d));
      }
    }
    for(size_t j = 0; j < cols; ++j) {

      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        edata.message = graph.vertex_data(img.vertid(i-1, j)).belief;
        edata.edgepot = lastudpotfactor;
        graph.add_edge(vertid, img.vertid(i-1, j), edata);
      }
      if(i+1 < img.rows()) {
        edata.message = graph.vertex_data(img.vertid(i+1, j)).belief;
        edata.edgepot = udpotfactor;
        graph.add_edge(vertid, img.vertid(i+1, j), edata);
      }
      if(j-1 < img.cols()) {
        edata.message = graph.vertex_data(img.vertid(i, j-1)).belief;
        edata.edgepot = lrpotfactor;
        graph.add_edge(vertid, img.vertid(i, j-1), edata);
      }
      if(j+1 < img.cols()) {
        edata.message = graph.vertex_data(img.vertid(i, j+1)).belief;
        edata.edgepot = lrpotfactor;
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
                       &RESAMPLE_FREQUENCY, RESAMPLE_FREQUENCY,
                       "resampling frequency");
  clopts.attach_option("mcmc",
                       &MCMCSTEPS, MCMCSTEPS,
                       "mcmc steps per resample");
  clopts.attach_option("logfile",
                       &logfile, std::string(""),
                       "log file");


  // set default scheduler type
  clopts.scheduler_type = "splash(100)";
  clopts.scope_type = "edge";

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



  // Create the graph --------------------------------------------------------->
  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);

  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  // collect the discretization points
  std::vector<double> discretizationpts;
  double d = -0.5;
  while (d <= 2.5) {
    d += 0.1;
    discretizationpts.push_back(d);
  }


  std::cout << std::endl;
  std::cout << "Arity: " << discretizationpts.size() << std::endl;


  construct_graph(gmmfile, numparticles, core.graph(),107,86,discretizationpts);
  
  
  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function", bp_update);

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
    pred(v2) = discretizationpts[vdata.belief.max_asg()];
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
    fout << "dbp\t" << gmmfile << "\t" << discretizationpts.size() << "\t" << 0 << "\t" << mae << "\t"
          << runtime << std::endl;
    fout.close();
  }
  
  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
} // End of main

