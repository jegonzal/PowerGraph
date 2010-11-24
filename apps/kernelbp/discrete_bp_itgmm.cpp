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

graphlab::binary_factor global_edge_pot;

// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
  graphlab::unary_factor message;
  graphlab::unary_factor old_message;
  
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
    tmp_msg.convolve(global_edge_pot, cavity);
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







void construct_graph(image& img,
                    const GaussianMixture<2> &edgepot,
                    const GaussianMixture<2> &nodepot,
                     double numparticles,
                     gl_types::graph& graph,
                    std::vector<double> discretizationpts,
                    bool usetruenodepot) {
  size_t numasgs = discretizationpts.size();
  // initialize a bunch of particles
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      vertex_data vdat;
      vdat.rounds = 0;

      // Set the node potential
      vdat.obs = img.pixel(i,j);
      vdat.potential.resize(numasgs);
      if (usetruenodepot == false) {
        GaussianMixture<1> np = gmm_conditional<2>()(nodepot, 01, vdat.obs);
        for (size_t p = 0;p < numasgs; ++p) {
          vdat.potential.logP(p) = log(np.likelihood(discretizationpts[p]));
        }
      }
      else {
        for(size_t pred = 0; pred < numasgs; ++pred) {
          vdat.potential.logP(pred) =
            -(vdat.obs - discretizationpts[pred])*(vdat.obs - discretizationpts[pred]) / (2.0 * 30 * 30);
        }
      }
      vdat.potential.normalize();

      vdat.belief = vdat.potential;
      vdat.belief.normalize();
      graph.add_vertex(vdat);
    } // end of for j in cols
  } // end of for i in rows


  // make the edge potential
  edge_data edata;
  global_edge_pot.resize(numasgs, numasgs);
  for (size_t p = 0; p < numasgs; ++p) {
    for (size_t q = 0; q < numasgs; ++q) {
      double d[2]; d[0] = discretizationpts[p]; d[1] = discretizationpts[q];
      global_edge_pot.logP(p,q) = log(edgepot.likelihood(d));
    }
  }

  //std::cout << edataedge_pot << std::endl;

  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {

      size_t vertid = img.vertid(i,j);
      if(i-1 < img.rows()) {
        edata.message = graph.vertex_data(img.vertid(i-1, j)).belief;
        edata.old_message = edata.message;
        graph.add_edge(vertid, img.vertid(i-1, j), edata);
      }
      if(i+1 < img.rows()) {
        edata.message = graph.vertex_data(img.vertid(i+1, j)).belief;
        edata.old_message = edata.message;
        graph.add_edge(vertid, img.vertid(i+1, j), edata);
      }
      if(j-1 < img.cols()) {
        edata.message = graph.vertex_data(img.vertid(i, j-1)).belief;
        edata.old_message = edata.message;
        graph.add_edge(vertid, img.vertid(i, j-1), edata);
      }
      if(j+1 < img.cols()) {
        edata.message = graph.vertex_data(img.vertid(i, j+1)).belief;
        edata.old_message = edata.message;
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
  bool potfromtest = false;
  std::string gmmfile= "";
  std::string inputfile = "";
  std::string logfile = "";
  bool learnpot = false;
  double ipf_update_damping;
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
  clopts.attach_option("potfromtest",
                       &potfromtest, false,
                       "Get edge potentials from test image");
clopts.attach_option("learnpot",
                       &learnpot, false,
                       "Learn an edge potential");
clopts.attach_option("ipfdamp",
                       &ipf_update_damping, 0.5,
                       "IPF damping");
std::string potfromfile;
clopts.attach_option("potfromfile",
                       &potfromfile, std::string(""),
                       "potfromfile");

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
  // collect the discretization points
  std::set<int> discretizationpts_set;
  for (size_t i = 0; i < trueimg.rows(); ++i) {
    for (size_t j = 0; j < trueimg.cols(); ++j) {
      discretizationpts_set.insert(size_t(trueimg.pixel(i,j)));
    }
  }
  std::vector<double> discretizationpts;
  std::map<size_t, size_t> discreterevmap;
  foreach(int i, discretizationpts_set) {
    discreterevmap[i] = discretizationpts.size();
    discretizationpts.push_back(i);
    std::cout << *(discretizationpts.end()-1) << " ";
  }
  std::cout << std::endl;
  std::cout << "Arity: " << discretizationpts.size() << std::endl;


  construct_graph(img, edgepot, nodepot, numparticles, core.graph(), discretizationpts, potfromtest);
  // if I am supposed to construct the potential from the test image....
  graphlab::binary_factor truebinarycount;
  double disagreecount = 0;
  if (potfromtest) {
    graphlab::binary_factor altedgepot;
    graphlab::unary_factor unarycount;
    size_t numasgs = discretizationpts.size();
    altedgepot.resize(numasgs, numasgs);
    truebinarycount.resize(numasgs, numasgs);
    unarycount.resize(numasgs);
    // we use the logP to store the actual values
    for (size_t i = 0;i < numasgs; ++i) {
      for (size_t j = 0;j < numasgs; ++j) {
        altedgepot.logP(i,j) = 1;  // initialize with something to avoid divide by 0s
        truebinarycount.logP(i,j) = 0.0;
      }
      unarycount.logP(i) = 1;
    }
    // count the pixels!
    for (size_t i = 0; i < trueimg.rows(); ++i) {
      for (size_t j = 0; j < trueimg.cols(); ++j) {
       assert(discreterevmap.find(trueimg.pixel(i,j)) != discreterevmap.end());
      }
    } 
    for (size_t i = 0; i < trueimg.rows(); ++i) {
      for (size_t j = 0; j < trueimg.cols(); ++j) {
        size_t color = discreterevmap[trueimg.pixel(i,j)]; 
        unarycount.logP(color)++;
        // check neighbors
        if (i > 1) {
         size_t color2 = discreterevmap[trueimg.pixel(i-1,j)]; 
         altedgepot.logP(color, color2)++; 
         truebinarycount.logP(color,color2)++;
         disagreecount += std::fabs((double)color - color2);
        }
        if (i < trueimg.rows() - 1) {
         size_t color2 = discreterevmap[trueimg.pixel(i+1,j)]; 
         altedgepot.logP(color, color2)++; 
         truebinarycount.logP(color,color2)++;
         disagreecount += std::fabs((double)color - color2);
        }
        if (j > 1) {
         size_t color2 = discreterevmap[trueimg.pixel(i,j-1)]; 
         altedgepot.logP(color, color2)++; 
         truebinarycount.logP(color,color2)++;
         disagreecount += std::fabs((double)color - color2);
        }
        if (j < trueimg.cols() - 1) {
         size_t color2 = discreterevmap[trueimg.pixel(i,j+1)]; 
         altedgepot.logP(color, color2)++; 
         truebinarycount.logP(color,color2)++;
         disagreecount += std::fabs((double)color - color2);
        }
      }
    }
    std::cout << unarycount;
    for (size_t i = 0;i < numasgs; ++i) {
      for (size_t j = 0;j < numasgs; ++j) {
        altedgepot.logP(i,j) = std::log(altedgepot.logP(i,j));
      }
      unarycount.logP(i) = std::log(unarycount.logP(i));
    }
    unarycount.normalize();
    altedgepot.normalize();
    for (size_t i = 0;i < numasgs; ++i) {
      for (size_t j = 0;j < numasgs; ++j) {
        altedgepot.logP(i,j) = altedgepot.logP(i,j) - unarycount.logP(i) - unarycount.logP(j);
      }
    }
    altedgepot.normalize();
    
    global_edge_pot = altedgepot;
    if (!(learnpot || potfromfile.length() > 0)) {
    for (size_t v = 0;v < core.graph().num_vertices(); ++v) {
      graphlab::unary_factor& pot = core.graph().vertex_data(v).potential;
      pot.times(unarycount);
      pot.normalize();
    }
    }
  }


  if (potfromfile.length() > 0) {
    std::ifstream fin;
    fin.open(potfromfile.c_str());
    graphlab::iarchive iarc(fin);
    iarc >> global_edge_pot;
    fin.close();
  }

  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function", bp_update);
  if (learnpot) {
    
    
    
    double bestpot = 0.5;
    double bestgap = 100000000000000;
    double potstrength = 0.5;
    for (size_t numlearniter = 0;numlearniter < 10; ++numlearniter) {
      global_edge_pot.set_as_laplace(potstrength);

      //  std::cout << global_edge_pot;
      std::cout << "Running the engine. " << std::endl;
      std::cout << "potential strength: " << potstrength << std::endl;
      size_t numasgs = discretizationpts.size();
      core.add_task_to_all(bp_update, 100.0);
      double runtime = core.start();
      //compute the edge counts
      double counts = 0;
      for (size_t v = 0;v < core.graph().num_vertices(); ++v) {
        vertex_data& v_data = core.graph().vertex_data(v);
        graphlab::edge_list ed = core.graph().out_edge_ids(v);
        foreach(graphlab::edge_id_t outeid, ed) {   
          vertex_data& w_data = core.graph().vertex_data(core.graph().target(outeid));
          for (size_t i = 0;i < numasgs; ++i) {
            for (size_t j = 0;j < numasgs; ++j) if (j!= i) counts += std::fabs((double)i-j) * std::exp(v_data.belief.logP(i) + w_data.belief.logP(j));
          }
        }
      }
      //std::cout << truebinarycount;
      //std::cout << counts;
      // compute the an IPF update    
      double gap = disagreecount - counts;
      //double step = 10 * gap / core.graph().num_edges();
      //potstrength -= (1 - ipf_update_damping) * step;
      potstrength = counts/disagreecount * potstrength;
      global_edge_pot.set_as_laplace(potstrength);
      gap = std::fabs(gap);
      std::cout << "learniter " << numlearniter << ": " << gap << std::endl;
      if (gap < bestgap) {
        bestpot = potstrength;
        bestgap = gap;
        std::cout << "better!" << std::endl;
        MAX_ITERATIONS += 20;
      }
      else if (gap / bestgap > 1.2) {
        ipf_update_damping = 1 - ((1 - ipf_update_damping) / 2);
        std::cout << "worse! Cutting learning rate to " << ipf_update_damping << std::endl;
        potstrength = bestpot;
        MAX_ITERATIONS += 20;
      }
      else {
        MAX_ITERATIONS += 20;
      }
    }
    std::string gmmfile2 = gmmfile.rfind("/") == std::string::npos ?
                                       gmmfile:
                                       gmmfile.substr(gmmfile.rfind("/")+1);

    std::ofstream fout(("pot"+gmmfile2+".bin").c_str());
    graphlab::oarchive oarc(fout);
    global_edge_pot.set_as_laplace(bestpot);
    oarc << global_edge_pot;
    fout.close();
    return 0;
  }
  
  std::cout << "Running the engine. " << std::endl;

  std::vector<graphlab::vertex_id_t> v;
  for (size_t i =0;i < core.graph().num_vertices(); ++i) v.push_back(i);
  std::random_shuffle(v.begin(), v.end());
  // Add the bp update to all vertices
  core.add_tasks(v, bp_update, 100.0);
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
    float a = discretizationpts[vdata.belief.max_asg()];
    img.pixel(v) = a;
  }

  double err = image_compare(trueimg, img);
  std::cout << "RMSE: " << err << std::endl;
  if (logfile.length() != 0) {
    ofstream fout;
    fout.open(logfile.c_str(), ios::app);
    gmmfile = gmmfile.rfind("/") == std::string::npos ?
                                       gmmfile:
                                       gmmfile.substr(gmmfile.rfind("/")+1);
    inputfile= inputfile.rfind("/") == std::string::npos ?
                                       inputfile:
                                       inputfile.substr(inputfile.rfind("/")+1);
    if (potfromtest == false) {
      fout << "dbp\t" << gmmfile << "\t" << inputfile << "\t" << discretizationpts.size() << "\t" << 0 << "\t" << err << "\t"
           << runtime << std::endl;
    }
    else {
      fout << "dbp2\t" << gmmfile << "\t" << inputfile << "\t" << discretizationpts.size() << "\t" << 0 << "\t" << err << "\t"
           << runtime << std::endl;
    }
    fout.close();
  }
  
  img.save("pred_map.pgm");


  for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    float a = vdata.belief.expectation();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = (a);
  }
  img.save("pred_exp.pgm");


  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
} // End of main

