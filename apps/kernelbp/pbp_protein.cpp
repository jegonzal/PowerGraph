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

#define PROPOSAL_STDEV 1

using namespace itpp;
using namespace std;



struct particle {
  std::vector<float> x;
  float weight;
};

float damping = 0.8;
size_t MCMCSTEPS = 10;
size_t RESAMPLE_FREQUENCY = 1000;

graphlab::atomic<size_t> proposal_count;
graphlab::atomic<size_t> accept_count;

// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
  std::vector<particle> message;
  GaussianMixture<6> p;
  
  float edgepot(const std::vector<float> &xsrc, const std::vector<float> &xdest) const {
    float f[6];
    for (size_t i = 0;i < 3; ++i) f[i] = xsrc[i];
    for (size_t i = 0;i < 3; ++i) f[i+3] = xdest[i];
    return p.likelihood(f);
  }

  float log_edgepot(const std::vector<float> &xsrc, const std::vector<float> &xdest) const {
    return std::log(edgepot(xsrc, xdest));
  }


}; // End of edge data

struct vertex_data: public graphlab::unsupported_serialize {
  // the belief at this vertex. Also provides the shape of all the incoming messages
  std::vector<particle> belief;
  uint rounds;
  GaussianMixture<3> p;
  float vertexpot(const std::vector<float> &x) const {
    return p.likelihood(x);
  }

  float log_vertexpot(const std::vector<float> &x) const {
    return log(p.likelihood(x));
  }

  void print() const{
    for (size_t i = 0;i < belief.size(); ++i) {
      std::cout << belief[i].weight << "\t";
    }
    std::cout << std::endl;
  }
  std::vector<float> max_asg() const{
    std::vector<float> ret;
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
double second_order_belief_ll(gl_types::iscope& scope, std::vector<float> v) {
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
  std::vector<std::vector<float> > newparticles;
  // sample from the old belief
  for (size_t i = 0;i < srcvdata.belief.size(); ++i) {
    // pick a particle to use
    double r = graphlab::random::rand01();
    // scan the CDF
    std::vector<float> cursample = srcvdata.belief[srcvdata.belief.size()-1].x;
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

      std::vector<float> proposal = cursample;
      for (size_t j = 0;j < proposal.size(); ++j) {
        proposal[j] += graphlab::random::gaussian_rand() * gaussian_proposal_stdev;;
      }
      proposal_count.inc();
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
      std::cout  << srcvdata.belief[i].weight << " ";
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
      std::cout  << v.belief[i].weight << " ";
    }
    std::cout << "\n";
  }


} // end of BP_update






void get_edge(std::string foldbase, size_t srcatype, size_t destatype, edge_data &ed) {
  std::string edgefname = foldbase + "_edge" + boost::lexical_cast<std::string>(srcatype) + 
        "_" + boost::lexical_cast<std::string>(destatype) + ".it";
  it_ifile f(edgefname.c_str());
  mat edgecenter, edgesigma, edgeweight;  
  vec doublevec;
  f >> Name("edge_ce") >> edgecenter;
  f >> Name("edge_alpha") >> doublevec; edgeweight = doublevec;
  f >> Name("edge_sigma") >> doublevec; edgesigma = doublevec;

  ed.p = GaussianMixture<6>(edgecenter, edgesigma, edgeweight);
}


mat construct_graph(std::string foldbase,
                    std::string testfile,
                    double numparticles,
                    gl_types::graph& graph) {
  // open the test file 
  it_ifile f(testfile.c_str());
  mat coord;
  mat feat;
  ivec atype;
  ivec chainlen;
  vec doublevec;
  mat nodecenter, nodesigma, nodeweight;

  f >> Name("coord") >> coord;
  f >> Name("atype") >> atype;
  f >> Name("chainlen") >> chainlen;
  f >> Name("like_ce") >> nodecenter;
  f >> Name("like_alpha") >> nodeweight;
  f >> Name("like_sigma") >> doublevec; nodesigma = doublevec;
  
  // build the chain
  for (size_t i = 0; i < chainlen[0]; ++i) {
    vertex_data vdat;
    vdat.rounds = 0;
    vdat.p = GaussianMixture<3>(nodecenter, nodesigma, mat(nodeweight.get_col(i)));
    for (size_t n = 0;n < numparticles; ++n) {
      particle p;
      p.x = vdat.p.sample();
      p.weight = 1.0/numparticles;
      vdat.belief.push_back(p);
    }
    graph.add_vertex(vdat);
    // add the edge to and from the node behind me
    if (i > 0) {
      // forward edge
      {
        edge_data ed;
        get_edge(foldbase, atype[i-1], atype[i], ed);
        ed.message = graph.vertex_data(i).belief;
        graph.add_edge(i-1, i, ed);
      }
      // backward edge
      {
        edge_data ed;
        get_edge(foldbase, atype[i], atype[i-1], ed);
        ed.message = graph.vertex_data(i-1).belief;
        graph.add_edge(i, i-1, ed);
      }
    }
  }
  return coord;
} // End of construct graph




















// MAIN =======================================================================>
int main(int argc, char** argv) {
  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  size_t iterations = 100;
  size_t numparticles = 100;
  std::string pred_type = "map";

  std::string testfile= "";
  std::string foldbase = "";
  std::string logfile = "";

  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Loopy BP image denoising");
  clopts.attach_option("iterations",
                       &iterations, iterations,
                       "Number of iterations");
  clopts.attach_option("particles",
                       &numparticles, numparticles,
                       "Number of particlesw");
  clopts.attach_option("foldbase",
                       &foldbase, std::string(""),
                       "gmm mixture file");
  clopts.attach_option("testfile",
                       &testfile, std::string(""),
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
  gl_types::core core;
  // Set the engine options
  core.set_engine_options(clopts);
  
  mat coord = construct_graph(foldbase, testfile, numparticles, core.graph());



  // Running the engine ------------------------------------------------------->
  core.scheduler().set_option(gl_types::scheduler_options::UPDATE_FUNCTION,
                              (void*)bp_update);

  std::cout << "Running the engine. " << std::endl;


  // Add the bp update to all vertices
  std::vector<graphlab::vertex_id_t> vforward, vbackward;
  for (size_t i = 0;i < core.graph().num_vertices(); ++i) {
    vforward.push_back(i);
    vbackward.push_back(core.graph().num_vertices() - i - 1);
  }
  
  graphlab::timer ti;
  ti.start();
  for (size_t iter = 0; iter < iterations; ++iter) {
    if (iter == iterations - 1) RESAMPLE_FREQUENCY = 100000;
    std::cout << "iteration " << iter << " forward" << std::endl;
    core.add_tasks(vforward, bp_update, 100.0);
    core.start();
    std::cout << "iteration " << iter << " backward" << std::endl;
    core.add_tasks(vbackward, bp_update, 100.0);
    core.start();
  }
  
  double runtime = ti.current_time();
  size_t update_count = core.last_update_count();
  std::cout << "Finished Running engine in " << runtime
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;

  // ok extract the predictions
  mat predcoord;
  predcoord.set_size(3, core.graph().num_vertices());
  for (size_t i = 0;i < core.graph().num_vertices(); ++i) {
    vertex_data &vdata = core.graph().vertex_data(i);
    std::vector<float> vasg = vdata.max_asg();
    for (size_t j = 0;j < 3; ++j) {
      predcoord(j,i) = vasg[j];
    }
  }
  vec err = sum(elem_mult(coord, predcoord), 1);
  double meanalign = mean(err);
  std::cout << "average alignment " << meanalign << std::endl;
  
  if (logfile.length() != 0) {
    ofstream fout;
    fout.open(logfile.c_str(), ios::app);
    foldbase = foldbase.rfind("/") == std::string::npos ?
                                       foldbase:
                                       foldbase.substr(foldbase.rfind("/")+1);
    
    fout << "pbp\t" << foldbase << "\t" << numparticles << "\t" << RESAMPLE_FREQUENCY << "\t" << meanalign << "\t"
    << runtime << std::endl;
    fout.close();
  }
  
  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
} // End of main

