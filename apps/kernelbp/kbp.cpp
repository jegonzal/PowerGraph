/*
 * kbp.cpp
 *
 * This file contains codes for kernel
 * belief propagation in a pairwise markov random field to denoise a
 * synthetic noisy image.
 *
 *  Created on: Oct 10, 2010
 *      Author: lesong
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <iostream>
#include <map>

#include <graphlab.hpp>
#include "image.hpp"

// include itpp
#include <itpp/itstat.h>
#include <itpp/itbase.h>
#include "image_compare.hpp"
// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

using namespace itpp;
using namespace std;
const size_t MAX_ITERATIONS = 30;
// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
  vec message;
  vec old_message;
  // vec lfeat_message;
}; // End of edge data

struct vertex_data: public graphlab::unsupported_serialize {
  vec belief;
  size_t rounds;
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;



mat pU;
mat pL;
mat lfeat;
ivec test;
mat ftest;
ivec isize;
vec m0;
vec fmu;
mat tmp_prod_msg;
mat tmp_ftest_pL;
size_t msgdim;
size_t basisno;
size_t testno;

mat fobs2;
ivec img1;
mat tmp_prod_msg_fobs2;
mat tmp_ftest_pL_fobs2;

double damping = 0.8;

// GraphLab Update Function ===================================================>

/** Construct denoising ising model based on the image */
void construct_graph(gl_types::graph& graph);

/**
 * The core belief propagation update function.  This update satisfies
 * the graphlab update_function interface.
 */
void bp_update(gl_types::iscope& scope,
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data);


// MAIN =======================================================================>
int main(int argc, char** argv) {
  std::cout << "kernel bp with denoising data" << std::endl;

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);

//  it_ifile f("../denoising_data2.it");
  it_ifile fmodel(argv[1]);
  fmodel >> Name("pU") >> pU;
  fmodel >> Name("pL") >> pL;
  fmodel >> Name("test") >> test;
  fmodel >> Name("fimgbasis") >> lfeat;
  fmodel >> Name("ftest") >> ftest;
  fmodel >> Name("isize") >> isize;
  fmodel >> Name("m0") >> m0;
  fmodel >> Name("fmu") >> fmu;
  fmodel >> Name("tmp_prod_msg") >> tmp_prod_msg;
  fmodel >> Name("tmp_ftest_pL") >> tmp_ftest_pL;

  msgdim = lfeat.rows();
  basisno = lfeat.cols();
  testno = ftest.cols();

  it_ifile fnoisy(argv[2]);
  fnoisy >> Name("fobs2") >> fobs2;
  fnoisy >> Name("img1") >> img1;

  tmp_prod_msg_fobs2 = tmp_prod_msg * fobs2;
  tmp_ftest_pL_fobs2 = tmp_ftest_pL * fobs2;

  // Create the graph --------------------------------------------------------->
  gl_types::graph graph;
  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(graph);
  size_t ncpus = 1;
  if (argc > 4) ncpus = atoi(argv[4]);
  // Create the engine -------------------------------------------------------->
  gl_types::iengine* engine =
    graphlab::engine_factory::new_engine("async",
                                         "sweep",
                                         "edge",
                                         graph,
                                         ncpus);
  if(engine == NULL) {
    std::cout << "Unable to construct engine!" << std::endl;
    return EXIT_FAILURE;
  }


  // Running the engine ------------------------------------------------------->

  // Set options that may be useful to some of the schedulers
  //size_t splash_size = 10000;
  //engine->get_scheduler().set_option(gl_types::scheduler_options::SPLASH_SIZE,
  //                                   (void*)&splash_size);
  //engine->get_scheduler().set_option(gl_types::scheduler_options::UPDATE_FUNCTION,
  //                                    (void*)bp_update);

  std::cout << "Running the engine. " << std::endl;
  graphlab::timer timer; timer.start();

  // Add the bp update to all vertices
  engine->add_task_to_all(bp_update, 100);
  // Starte the engine
  engine->start();
  double runtime = timer.current_time();
  size_t update_count = engine->last_update_count();
  std::cout << "Finished Running engine in " << runtime
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;


  std::cout << "Computing estimation error" << std::endl;
  ivec pred(isize[0]*isize[1]);
  image predimg(isize[0], isize[1]);
  for (size_t v = 0; v < graph.num_vertices(); ++v) {
	  const vertex_data& vdata = graph.vertex_data(v);
	  pred(v) = test(max_index(vdata.belief));
    predimg.pixel(v) = pred(v);
  }
  std::cout << "Mean absolute error: " << mean(abs(pred - img1)) << std::endl;

  std::cout << "Creating a synethic image. " << std::endl;
  
  image trueimg(isize[0], isize[1]);
  for (size_t i = 0;i < img1.size(); ++i) {
    trueimg.pixel(i) = img1(i);
  }
  trueimg.save("source_img.pgm");
  double err = image_compare(trueimg, predimg);
  predimg.save("result.pgm");
  std::cout << "RMSE: " << err << std::endl;
  std::string logfile;
  if (argc > 3) logfile = argv[3];
  if (logfile.length() != 0) {
    ofstream fout;
    fout.open(logfile.c_str(), ios::app);
    std::string gmmfile = argv[1];
    std::string inputfile = argv[2];
    gmmfile = gmmfile.rfind("/") == std::string::npos ?
                                       gmmfile:
                                       gmmfile.substr(gmmfile.rfind("/")+1);
    inputfile= inputfile.rfind("/") == std::string::npos ?
                                       inputfile:
                                       inputfile.substr(inputfile.rfind("/")+1);

    fout << "kbp\t" << gmmfile << "\t" << inputfile << "\t" << 0 << "\t" << 0 << "\t" << err << "\t" << runtime << std::endl;
    fout.close();
  }
  

  std::cout << "Done!" << std::endl;

  if(engine != NULL) delete engine;

  return EXIT_SUCCESS;
} // End of main

// Implementations
// ============================================================>
void bp_update(gl_types::iscope& scope,
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data) {
  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  if (v_data.rounds >= MAX_ITERATIONS) return;
  v_data.rounds++; 
 
  graphlab::vertex_id_t vid = scope.vertex();

  vec prod_message = tmp_prod_msg_fobs2.get_col(vid);
  v_data.belief = tmp_ftest_pL_fobs2.get_col(vid);

  // Get the in and out edges by reference
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size()); // Sanity check

  // Flip the old and new messages to improve safety when using the
  // unsynch scope
  std::map<graphlab::edge_id_t, vec> lfeat_message;
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the in and out edge data
    edge_data& e_data = scope.edge_data(ineid);
    // Since we are about to receive the current message make it the
    // old message
    e_data.old_message = e_data.message;
    lfeat_message[ineid] = lfeat.transpose() * e_data.old_message;

    elem_mult_inplace(
      ftest.transpose() * e_data.old_message,
      v_data.belief
    );

    elem_mult_inplace(
      lfeat_message[ineid],
      prod_message
    );
  }

  elem_mult_inplace(
    ftest.transpose() * fmu,
    v_data.belief
  );

  /*
  if (vid == 0) {
    it_file f2("outfile2.it");
    f2<<Name("vdata_belief")<<v_data.belief;
    f2<<Name("prod_message")<<prod_message<<flush;
    exit(0);
  }
  */

  vec cavity;
  cavity.set_size(basisno);
  vec tmp_msg;
  tmp_msg.set_size(msgdim);
  double residual = 0;
  for(size_t i = 0; i < in_edges.size(); ++i) {
    // Get the edge ids
    graphlab::edge_id_t outeid = out_edges[i];
    graphlab::edge_id_t ineid = in_edges[i];
    // CLEVER HACK: Here we are expoiting the sorting of the edge ids
    // to do fast O(1) time edge reversal
    assert(scope.target(outeid) == scope.source(ineid));
    // const edge_data& in_edge = scope.edge_data(ineid);
    edge_data& out_edge = scope.edge_data(outeid);

    // Compute cavity
    elem_div_out(prod_message, lfeat_message[ineid], cavity);

    // convolve cavity with the edge factor storing the result in the
    // temporary message
    tmp_msg = pU * cavity;
    // tmp_msg /= sqrt(sum_sqr(tmp_msg));

    // Damp the message
    tmp_msg = damping * out_edge.old_message
                     + (1-damping) * tmp_msg;
    tmp_msg /= sqrt(sum_sqr(tmp_msg));

    // Compute message residual
    double r = sqrt(sum_sqr(tmp_msg - out_edge.old_message));

    // Assign the out message
    out_edge.message = tmp_msg;
    residual += r;
    if (v_data.rounds < MAX_ITERATIONS) {
      gl_types::update_task task(scope.vertex(), bp_update);
      scheduler.add_task(task, 1.0);
    }
  }
} // end of BP_update

void construct_graph(gl_types::graph& graph) {
  // Construct a single blob for the vertex data
  vertex_data vdata;
  vdata.rounds = 0;
  vdata.belief.set_size(testno);
  // Add all the vertices
  for (size_t i = 0; i < (size_t) isize[0] * isize[1]; i++) {
    size_t vertid = graph.add_vertex(vdata);
    assert(vertid == i);
  }

  // Add the edges
  edge_data edata;
  edata.message = m0;
  edata.old_message = edata.message;
  // edata.lfeat_message.set_size(basisno);
  // edata.lfeat_message.ones();

  for(size_t j = 0; j < (size_t) isize[1]; ++j) {
    for(size_t i = 0; i < (size_t) isize[0]; ++i) {
      size_t vertid = j*isize[0] + i;
      if(i-1 < (size_t)isize[0]) {
        graph.add_edge(vertid, j*isize[0] + (i-1), edata);
      }
      if(i+1 < (size_t)isize[0]) {
        graph.add_edge(vertid, j*isize[0] + (i+1), edata);
      }
      if(j-1 < (size_t)isize[1]) {
        graph.add_edge(vertid, (j-1)*isize[0] + i, edata);
      } if(j+1 < (size_t)isize[1]) {
        graph.add_edge(vertid, (j+1)*isize[0] + i, edata);
      }
    } // end of for j in cols
  } // end of for i in rows

  graph.finalize();

} // End of construct graph



