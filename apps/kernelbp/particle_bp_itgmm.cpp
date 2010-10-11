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

float damping = 0.8;
size_t MCMCSTEPS = 10;
size_t RESAMPLE_FREQUENCY = 1000;
size_t MAX_ITERATIONS;
graphlab::atomic<size_t> proposal_count;
graphlab::atomic<size_t> accept_count;

// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
  std::vector<particle> message;
  GaussianMixture<2> p;
  
  float edgepot(float xsrc, float xdest) const {
    double arr[2];
    arr[0] = xsrc;
    arr[1] = xdest;
    return p.likelihood(arr);
  }

  float log_edgepot(float xsrc, float xdest) const {
    return std::log(edgepot(xsrc, xdest));
  }


}; // End of edge data

struct vertex_data: public graphlab::unsupported_serialize {
  // the belief at this vertex. Also provides the shape of all the incoming messages
  std::vector<particle> belief;
  float obs;
  uint rounds;
  GaussianMixture<1> p;
  float vertexpot(float x) const {
    return p.likelihood(x);
  }

  float log_vertexpot(float x) const {
    return std::log(p.likelihood(x));
  }

  float average() const {
    float v = 0;
    for (size_t i = 0;i < belief.size(); ++i) {
      v = v  + belief[i].weight * belief[i].x;
    }
    return v;
  }

  void print() const{
    for (size_t i = 0;i < belief.size(); ++i) {
      std::cout << belief[i].x << " " << belief[i].weight << "\t";
    }
    std::cout << std::endl;
  }
  float max_asg()  const{
    float ret = 0;
    float best = 0;
    for (size_t i = 0;i < belief.size(); ++i) {
      if (belief[i].weight > best) {
        ret = belief[i].x;
        best = belief[i].weight;
      }
    }
    return ret;
  }
};



typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;


/**
  Generates a spherical guassian of a particular stddev
  Returns log likelihood
*/
double gaussianrng(float &v,  double stdev) {
  // box muller generator
  double U0 = graphlab::random::rand01() + DBL_MIN;
  double U1 = graphlab::random::rand01() + DBL_MIN;
  double z0 = sqrt(-2 * log(U0)) * cos(2 * M_PI * U1);

  v = z0;
  return (v * v) / (2 * stdev * stdev);;
}

/**
  Makes the largest particle have a size of 1
*/
void normalize_particles(std::vector<particle> &particles) {
  float totalweight = 0;
  for (size_t i = 0;i < particles.size(); ++i) {
    totalweight += particles[i].weight;
  }
  for (size_t i = 0;i < particles.size(); ++i) {
    particles[i].weight /= totalweight;
    // set a minimum threshold
    if (particles[i].weight < 1E-12) particles[i].weight = 1E-12;
  }
}

/**
  Computes the message along the edge.
  The message requires the shape of the target vertex
*/
void sample_message(gl_types::iscope& scope,
                    graphlab::edge_id_t edge) {
  // source and dest vertices
  const vertex_data &srcvdata = scope.const_vertex_data();
  const vertex_data &destvdata = scope.const_neighbor_vertex_data(scope.target(edge));
  // in and out edges
  edge_data &outedgedata = scope.edge_data(edge);
  const edge_data &inedgedata = scope.const_edge_data(scope.reverse_edge(edge));

  // get the shape of the particles
  outedgedata.message = destvdata.belief;
  // update the weights
  for (size_t i = 0; i < outedgedata.message.size(); ++i) {
    double oldweight = outedgedata.message[i].weight;
    outedgedata.message[i].weight = 0;
    for (size_t j = 0;j < inedgedata.message.size(); ++j) {
      double scaleweight = srcvdata.belief[j].weight / inedgedata.message[j].weight;
      if (scaleweight > 1E-7) {
        // outgoing potential * src belief / by incoming message
        outedgedata.message[i].weight  +=
                  outedgedata.edgepot(srcvdata.belief[j].x, outedgedata.message[i].x) * scaleweight;
      }
    }
    outedgedata.message[i].weight = outedgedata.message[i].weight * damping
                                    + oldweight * (1.0 - damping);
  }
 normalize_particles(outedgedata.message);
}


/**
Updates the belief on the current vertex.
Requires that the current vertex has the same particles
as all the incoming messages
*/
void update_belief(gl_types::iscope& scope) {
  // source and dest vertices
  vertex_data &vdata = scope.vertex_data();
  for (size_t i = 0;i < vdata.belief.size(); ++i) {
    vdata.belief[i].weight = vdata.vertexpot(vdata.belief[i].x);
  }
  normalize_particles(vdata.belief);
  foreach(graphlab::edge_id_t ineid, scope.in_edge_ids()) {
    const edge_data &inedgedata = scope.const_edge_data(ineid);
    for (size_t i = 0;i < vdata.belief.size(); ++i) {
      vdata.belief[i].weight *= inedgedata.message[i].weight;
    }
    normalize_particles(vdata.belief);
  }
}



/**
Estimates  log P(v) by reevaluating all the incoming messages
*/
double second_order_belief_ll(gl_types::iscope& scope, float v) {
  graphlab::edge_list out_edges = scope.out_edge_ids();
  double ll = 0;
  foreach(graphlab::edge_id_t outeid, out_edges) {
    const edge_data &outedgedata = scope.const_edge_data(outeid);
    const vertex_data &othervdata = scope.const_neighbor_vertex_data(scope.target(outeid));
    double w = 0;
    for (size_t j = 0;j < othervdata.belief.size(); ++j) {
      w += outedgedata.edgepot(v, othervdata.belief[j].x) *
                        othervdata.belief[j].weight /
                        outedgedata.message[j].weight;
    }
    ll += log(w);
  }
  ll += (scope.vertex_data().log_vertexpot(v));
  return ll;
}


double eff_particles(gl_types::iscope &scope) {
// lets move all the particles around using MCMC
  vertex_data &srcvdata = scope.vertex_data();

  normalize_particles(srcvdata.belief);
  // get effective number of particle weight
  double wsqr = 0;
  for (size_t i = 0;i < srcvdata.belief.size(); ++i) {
    wsqr += srcvdata.belief[i].weight * srcvdata.belief[i].weight;
  }
  return 1.0 / wsqr;
}

/**
  Resamples all the partciles at the current scope
*/
float resample_particles(gl_types::iscope& scope,
                        double gaussian_proposal_stdev) {
  int acceptcount = 0;
  int proposecount = 0;

  // lets move all the particles around using MCMC
  vertex_data &srcvdata = scope.vertex_data();

  // if the number of particles is fine. don't resample
  if (eff_particles(scope) > srcvdata.belief.size() / 2 + 1) return 0;
  update_belief(scope);
  std::vector<float> newparticles;
  // sample from the old belief
  for (size_t i = 0;i < srcvdata.belief.size(); ++i) {
    // pick a particle to use
    double r = graphlab::random::rand01();
    // scan the CDF
    float cursample = srcvdata.belief[srcvdata.belief.size()-1].x;
    double accweight = 0;
    for (size_t j = 0;j < srcvdata.belief.size(); ++j) {
      accweight += srcvdata.belief[j].weight;
      if (accweight >= r) {
        cursample = srcvdata.belief[j].x;
        break;
      }
    }
    double curll = second_order_belief_ll(scope, cursample);
    for (size_t i = 0; i < MCMCSTEPS; ++i) {
      float gauss = graphlab::random::gaussian_rand() * gaussian_proposal_stdev;
      proposal_count.inc();
      float proposal = cursample + gauss;
      // compute forward probability
      double proposell = second_order_belief_ll(scope, proposal);
      // since the gaussian is symmetric...
      double acceptll = proposell - curll;
      proposecount++;
      bool accept = false;
      if (acceptll >= 0)  {
        accept = true;
      }
      else {
        // acceptreject
        acceptll = std::exp(acceptll);
        if (graphlab::random::rand01() <= acceptll)  accept = true;
      }
      if (accept) {
        acceptcount++;
        accept_count.inc();
        curll = proposell;
        cursample = proposal;
      }
    }

    newparticles.push_back(cursample);
  }

  // fill the new belief
  for (size_t i = 0;i < srcvdata.belief.size(); ++i) {
    srcvdata.belief[i].x = newparticles[i];
    srcvdata.belief[i].weight = 1.0 / srcvdata.belief.size();
  }

  
  if(scope.vertex() == 0) {
    std::cout << "vertex 0 particles: ";
    for (size_t i = 0;i < srcvdata.belief.size(); ++i) {
      std::cout  << srcvdata.belief[i].x << ":" << srcvdata.belief[i].weight << " ";
    }
    std::cout << "\n";
  }
  graphlab::edge_list in_edges = scope.in_edge_ids();
  foreach(graphlab::edge_id_t ineid, in_edges) {
    sample_message(scope, ineid);
  }
  return float(acceptcount) / proposecount;
}




/**
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.
 */
void bp_update(gl_types::iscope& scope,
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data) {
  // resample the particles here.

  vertex_data &v = scope.vertex_data();
  if (v.rounds >= MAX_ITERATIONS) return;
  if (scope.vertex() == 0) {
    std::cout << "vertex " << scope.vertex()
              << " round " << v.rounds << std::endl;
  }
  if ((v.rounds+1) % RESAMPLE_FREQUENCY == 0) {
    double resampleround = (v.rounds+1) / RESAMPLE_FREQUENCY;
    double stddev = 1.0 / (resampleround) * 0.8 + (PROPOSAL_STDEV - 0.8);
    resample_particles(scope, stddev);
  }
  update_belief(scope);
  if (scope.vertex() % 1000 == 0){
    std::cout << accept_count.value << " / " << proposal_count.value << " = " << double(accept_count.value) / proposal_count.value << ": "
    << eff_particles(scope) << std::endl;
  }

  // compute the new outgoing messages
  graphlab::edge_list out_edges = scope.out_edge_ids();
  foreach(graphlab::edge_id_t outeid, out_edges) {
    sample_message(scope, outeid);
  }

  v.rounds++;
 /* if (v.rounds < MAX_ITERATIONS && accprop <= 0.3) {
    foreach(graphlab::edge_id_t outeid, out_edges) {
      gl_types::update_task task(scope.target(outeid), bp_update);
      scheduler.add_task(task, 1.0);
    }
  }*/
  if(scope.vertex() == 0) {
    std::cout << "vertex 0 particles: ";
    for (size_t i = 0;i < v.belief.size(); ++i) {
      std::cout  << v.belief[i].x << ":" << v.belief[i].weight << " ";
    }
    std::cout << "\n";
  }
  if (v.rounds < MAX_ITERATIONS) {
    gl_types::update_task task(scope.vertex(), bp_update);
    scheduler.add_task(task, 1.0);
  }

} // end of BP_update








void construct_graph(image& img,
                    const GaussianMixture<2> &edgepot,
                    const GaussianMixture<2> &nodepot,
                     double numparticles,
                     gl_types::graph& graph) {
                     
  // initialize a bunch of particles
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      vertex_data vdat;
      vdat.rounds = 0;

      // Set the node potential
      vdat.obs = img.pixel(i,j);
      vdat.p = gmm_conditional<2>()(nodepot, 1, vdat.obs);
      vdat.p.simplify();
      for (size_t n = 0;n < numparticles; ++n) {
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
      }
      graph.add_vertex(vdat);

    } // end of for j in cols
  } // end of for i in rows



  // Add the edges
  
  edge_data edata;
  edata.p = edgepot;

  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {

      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        edata.message = graph.vertex_data(img.vertid(i-1, j)).belief;
        graph.add_edge(vertid, img.vertid(i-1, j), edata);
      }
      if(i+1 < img.rows()) {
        edata.message = graph.vertex_data(img.vertid(i+1, j)).belief;
        graph.add_edge(vertid, img.vertid(i+1, j), edata);
      }
      if(j-1 < img.cols()) {
        edata.message = graph.vertex_data(img.vertid(i, j-1)).belief;
        graph.add_edge(vertid, img.vertid(i, j-1), edata);
      }
      if(j+1 < img.cols()) {
        edata.message = graph.vertex_data(img.vertid(i, j+1)).belief;
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
  std::string inputfile = "";
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
  clopts.attach_option("inputfile",
                       &inputfile, std::string(""),
                       "input file");
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

  // load the potentials mixture components
  it_ifile f(gmmfile.c_str());

  // weigghts
  mat edgecenter, edgesigma, edgeweight;
  mat nodecenter, nodesigma, nodeweight;
  ivec truedata;
  ivec imgsize;
  // intermediate types to use...
  imat integermat;
  vec doublevec;
  f >> Name("edge_ce") >> integermat;   edgecenter = to_mat(integermat);
  f >> Name("edge_alpha") >> doublevec; edgeweight = doublevec;
  f >> Name("edge_sigma") >> doublevec; edgesigma = doublevec;

  f >> Name("like_ce") >> nodecenter;
  f >> Name("like_alpha") >> doublevec; nodeweight = doublevec;
  f >> Name("like_sigma") >> doublevec; nodesigma = doublevec;
  f >> Name("img1") >> truedata;
  f >> Name("isize") >> imgsize;

  size_t rows = imgsize(0);
  size_t cols = imgsize(1);
  std::cout << "Image size is "
            << rows << " x " << cols << std::endl;
  // make the GMMs
  GaussianMixture<2> edgepot(edgecenter, edgesigma, edgeweight);
  GaussianMixture<2> nodepot(nodecenter, nodesigma, nodeweight);
  edgepot.simplify();
  edgepot.print();
  std::cout << "node pot: " << std::endl;
  nodepot.print();
  // convert the true image to an image
  image trueimg(rows, cols);
  for (size_t i = 0;i < truedata.size(); ++i) {
    trueimg.pixel(i) = truedata(i);
  }

  // load the observations
  it_ifile imgfile(inputfile.c_str());
  vec observations;
  imgfile >> Name("obs2") >> observations;
  // convert observations to an image
  image img(rows, cols);
  for (size_t i = 0;i < observations.size(); ++i) {
    img.pixel(i) = observations(i);
  }
  img.save("noisy.pgm");
  trueimg.save("source_img.pgm");

  // Create the graph --------------------------------------------------------->
  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);

  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;

  construct_graph(img, edgepot, nodepot, numparticles, core.graph());



  // Running the engine ------------------------------------------------------->
  core.scheduler().set_option(gl_types::scheduler_options::UPDATE_FUNCTION,
                              (void*)bp_update);

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
    img.pixel(v) = a;
  }
  
  double rmse = image_compare(trueimg, img);
  std::cout << "RMSE: " << rmse << std::endl;
  if (logfile.length() != 0) {
    ofstream fout;
    fout.open(logfile.c_str());
    fout << gmmfile << "\t" << inputfile << "\t" << numparticles << "\t" << RESAMPLE_FREQUENCY << "\t" << rmse  << std::endl;
    fout.close();
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


  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
} // End of main

