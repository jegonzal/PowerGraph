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
#include "kde.h"
#include <float.h>
#include "../kernelbp/image.hpp"
#include "prodSampleEpsilon.hpp"
#include <set>

// include itpp
#include <itpp/itstat.h>
#include <itpp/itbase.h>

// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>
#define NDEBUG

#define PROPOSAL_STDEV 20

int NSAMP =12;
double EPSILON =1e-5;

int ROWS=107;
int COLS=86;


using namespace itpp;
using namespace std;

struct particle {
  float x;
  float weight;
};

float SIGMA;
float LAMBDA;
float damping = 0.8;
size_t MAX_ITERATIONS = 2;
graphlab::atomic<size_t> proposal_count;
graphlab::atomic<size_t> accept_count;
int iiter = 0;
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
    //for (size_t i = 0; i < 256; ++i) std::cout << colormap[i] << " ";
    //std::cout << std::endl;
    // round the infered image
   /* for (size_t i = 0; i < infered.rows(); ++i) {
      for (size_t j = 0; j < infered.cols(); ++j) {
        if (infered.pixel(i,j) >= 255) infered.pixel(i,j) = 255;
        if (infered.pixel(i,j) < 0) infered.pixel(i,j) = 0;
        infered.pixel(i,j) = colormap[(size_t)(infered.pixel(i,j) + 0.5)];
      }
    }*/
    
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
    //for (size_t i = 0; i < 256; ++i) std::cout << colormap[i] << " ";
    //std::cout << std::endl;
    // round the infered image
   /* for (size_t i = 0; i < infered.rows(); ++i) {
      for (size_t j = 0; j < infered.cols(); ++j) {
        if (infered.pixel(i,j) >= 255) infered.pixel(i,j) = 255;
        if (infered.pixel(i,j) < 0) infered.pixel(i,j) = 0;
        infered.pixel(i,j) = colormap[(size_t)(infered.pixel(i,j) + 0.5)];
      }
    }*/
    
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
  if (debug && vid%1000000 == 0){
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
      kdes.insert(kdes.begin(), marg);//important: has to be first!     

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
   if (v_data.rounds == MAX_ITERATIONS){
	if (debug && vid%100000 == 0)
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
	   printf("belief node %d is\n", vid);
           m.matlab_print(); printf("\n");
      }


   }
  /*
  if (v_data.rounds < MAX_ITERATIONS) {
    gl_types::update_task task(scope.vertex(), bp_update);
    scheduler.add_task(task, 1.0);
  }*/


} // end of BP_update





void construct_graph(image& img,
                     kde & edge_pot,
                     gl_types::graph& graph) {
  // initialize a bunch of particles


  
  for(size_t i = 0; i < img.rows(); ++i) {
    for(size_t j = 0; j < img.cols(); ++j) {
      vertex_data vdat;
      vdat.rounds = 0;

      // Set the node potentiala
      vec cent = zeros(2);
      cent[0] = img.pixel(i,j);
      cent[1] = img.pixel(i,j);
      mat cent2 = cent; cent2 = transpose(cent2);
      vec bw = "30 30";
      mat bw2 = bw; bw2 = transpose(bw2);
      vec wght = "1 1";
      vdat.obs = kde(cent2, bw2, wght);
      vdat.bel = vdat.obs;
      if (i == 0 && j == 0)
       vdat.obs.matlab_print();
      graph.add_vertex(vdat);
      vdat.obs.verify();
      vdat.bel.verify();

    } // end of for j in cols
  } // end of for i in rows

  // Add the edges
  edge_data edata;
  kde _edge_pot = kde("255    0  226    0  141  113    0  198   85    0  170   56; \
-                     255    0    0  226  141    0  113  198   85  170    0   28",
"70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001;  \
 70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001      70.7001",
  "0.24635     0.17462           1e-5           1e-5     0.18841           1e-5           1e-5     0.15803     0.14531           1e-5           1e-5    0.087277"); 
  edata.edge_pot = edge_pot;
  edata.edge_pot.matlab_print();
  _edge_pot.verify(); 

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
  clopts.attach_option("epsilon",
                       &EPSILON, EPSILON,
                       "epsilon");
  clopts.attach_option("gmmfile",
                       &gmmfile, std::string(""),
                       "gmm mixture file");
  clopts.attach_option("inputfile",
                       &inputfile, std::string(""),
                       "input file");


  clopts.scheduler_type = "round_robin";
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
  mat edgesigma2 = edgesigma;
  edgesigma2 = transpose(edgesigma2);
  mat edgeweight2 = edgeweight;
  edgeweight2 = transpose(edgeweight2); 
  if (edgesigma2.cols() > edgecenter.cols())
	edgesigma2 = edgesigma2(0,0,0,edgecenter.cols()-1);
  kde edge_pot = kde(edgecenter, edgesigma2, edgeweight2);





// convert the true image to an image
  image trueimg(rows, cols);
  for (size_t i = 0;i < truedata.size(); ++i) {
    trueimg.pixel(i) = truedata(i);
  }

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

   construct_graph(img, edge_pot, core.graph());



  // Running the engine ------------------------------------------------------->
  core.sched_options().add_option("update_function", bp_update);
  core.sched_options().add_option("max_iterations", MAX_ITERATIONS);

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
    float a = vdata.bel.max();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = size_t(a);
  }
  double err = sqrt(image_compare_rmse(trueimg, img));
  double err2 = image_compare_mae(trueimg, img);
  std::cout << "RMSE: " << err << " MAE: "<< err2<<std::endl;

  //img.save("pred_map.pgm");


 /* for(size_t v = 0; v < core.graph().num_vertices(); ++v) {
    const vertex_data& vdata = core.graph().vertex_data(v);
    float a = vdata.average();
    if (a < 0) a = 0;
    if (a > 255) a = 255;
    img.pixel(v) = (a);
  }
  img.save("pred_exp.pgm");
*/ //TODO?

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
} // End of main

