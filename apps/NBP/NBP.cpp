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
#include "image.hpp"
#include "kde.h"

// include itpp
#include <itpp/itstat.h>
#include <itpp/itbase.h>


#include <fstream>
#include <cmath>
#include <cstdio>
#include <cfloat>
//#include <graphlab/graph/graph.hpp>

#include <graphlab/macros_def.hpp>


using namespace itpp;
using namespace std;

// STRUCTS (Edge and Vertex data) =============================================>
struct edge_data: public graphlab::unsupported_serialize {
//  vec message;
//  vec old_message;
  kde msg;
  // vec lfeat_message;
  int update_count;
}; // End of edge data

struct vertex_data: public graphlab::unsupported_serialize {
  kde obs;
  //vec belief;
  size_t row, col;
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

double damping = 0.90;

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

  std::string prefix(argv[1]);
  std::string fname = prefix + "0.it";

  std::cout << "file to load: " << fname << std::endl;

  it_ifile part0(fname.c_str());
  part0 >> Name("isize") >> isize;
  part0 >> Name("m0") >> m0;
  part0 >> Name("predy_k") >> predy_k;
  part0 >> Name("predfy_k") >> predfy_k;
  part0 >> Name("prod_msg0") >> prod_msg0;
  part0 >> Name("belief0") >> belief0;
  part0 >> Name("testy") >> testy;
  part0 >> Name("truey") >> truey;
  part0 >> Name("cfybasis") >> cfybasis;
  part0 >> Name("ftesty") >> ftesty;

  pUud = new mat[isize[0]];
  pUdu = new mat[isize[0]];
  pUlr = new mat[isize[0]];
  pUrl = new mat[isize[0]];

  for (size_t i = 0; i < (size_t) isize[0]; i++) {
	  std::stringstream ss;
	  ss << i + 1;
	  fname = prefix + ss.str() + ".it";
	  it_ifile part(fname.c_str());

	  std::cout << "file to load: " << fname << std::endl;

	  if (i > 0) part >> Name("ipUud") >> pUud[i];
	  if (i < (size_t) (isize[0]-1)) part >> Name("ipUdu") >> pUdu[i];
	  part >> Name("ipUlr") >> pUlr[i];
	  part >> Name("ipUrl") >> pUrl[i];

	  part.close();
  }

  msgdim = cfybasis.rows();
  basisno = cfybasis.cols();
  testno = ftesty.cols();

  // Create the graph --------------------------------------------------------->
  gl_types::graph graph;
  std::cout << "Constructing pairwise Markov Random Field. " << std::endl;
  construct_graph(graph);

  // Create the engine -------------------------------------------------------->
  gl_types::iengine* engine =
    graphlab::engine_factory::new_engine("threaded",
                                         "fifo",
                                         "edge",
                                         graph,
                                         1);
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
  engine->get_scheduler().add_task_to_all(bp_update, 100);
  // Starte the engine
  engine->start();

  double runtime = timer.current_time();
 /* size_t update_count = engine->get_profiling_info("update_count");
  std::cout << "Finished Running engine in " << runtime
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;
*/
  std::cout << "Computing estimation error" << std::endl;
  vec pred(isize[0]*isize[1]);
/*  for (size_t v = 0; v < graph.num_vertices(); ++v) {
	  const vertex_data& vdata = graph.vertex_data(v);
	  pred(v) = testy(max_index(vdata.belief));
  }
  std::cout << "Mean absolute error: " << mean(abs(pred - truey)) << std::endl;
*/
  delete [] pUud;
  delete [] pUdu;
  delete [] pUlr;
  delete [] pUrl;

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
  graphlab::vertex_id_t vid = scope.vertex();

  //vec prod_message = prod_msg0.get_col(vid);
  //v_data.belief = belief0.get_col(vid);

  // Get the in and out edges by reference
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size()); // Sanity check

/*  // Flip the old and new messages to improve safety when using the
  // unsynch scope
  std::map<graphlab::edge_id_t, vec> lfeat_message;
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the in and out edge data
    edge_data& e_data = scope.edge_data(ineid);
    // Since we are about to receive the current message make it the
    // old message
    e_data.old_message = e_data.message;
    lfeat_message[ineid] = cfybasis.transpose() * e_data.old_message;

    elem_mult_inplace(
      ftesty.transpose() * e_data.old_message,
      v_data.belief
    );

    elem_mult_inplace(
      lfeat_message[ineid],
      prod_message
    );
  }*/

//  std::cout << "vid: " << vid << std::endl;

//  if (vid == 0) {
//    it_file f2("outfile2.it");
//    f2<<Name("vdata_belief")<<v_data.belief;
//    f2<<Name("prod_message")<<prod_message<<flush;
//    exit(0);
//  }

  vec cavity;
  cavity.set_size(basisno);
  vec tmp_msg;
  tmp_msg.set_size(msgdim);

//  if (v_data.row == 0 && v_data.col == 0) {
//	  std::cout << "in edge size: " << in_edges.size() << std::endl;
//  }

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
    //elem_div_out(prod_message, lfeat_message[ineid], cavity);

    vertex_data& tmp_v_data = scope.neighbor_vertex_data(scope.target(outeid));
//      scheduler.add_task(task, -vid);
    }

} // end of BP_update

void construct_graph(gl_types::graph& graph) {
  // Construct a single blob for the vertex data
  vertex_data vdata;
  //vdata.belief.set_size(testno);
  // Add all the vertices
  for(size_t j = 0; j < (size_t) isize[1]; ++j) {
    for(size_t i = 0; i < (size_t) isize[0]; ++i) {
    	size_t tmpvertid = j*isize[0] + i;
    	vdata.row = i;
    	vdata.col = j;
        size_t vertid = graph.add_vertex(vdata);
        assert(vertid == tmpvertid);
    }
  }

  // Add the edges
  edge_data edata;
  //edata.message = m0;
  //edata.old_message = edata.message;
  edata.update_count = 0;

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



#include <graphlab/macros_undef.hpp>
