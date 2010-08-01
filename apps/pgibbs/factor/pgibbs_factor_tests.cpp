/**
 *
 * Parallel blocked gibbs using graphlab
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <set> 
#include <algorithm>
#include <limits>
#include <cmath>




#include <graphlab.hpp>



// Image reading/writing code
#include "data_structures.hpp"
#include "update_functions.hpp"



// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>




void test_toy_graph() {

  
  variable_t v0(0, 2);
  variable_t v1(1, 3);
  variable_t v2(2, 4);
  variable_t v3(3, 2);
  
  factor_t fact1(domain_t(v0, v1));
  factor_t fact2(domain_t(v1, v2));
  factor_t fact3(domain_t(v2, v3));

  fact1.set_as_agreement(3);
  fact2.set_as_laplace(2);
  fact3.set_as_laplace(10);


  // Consturct the model
  factorized_model model;
  model.add_factor(fact1);
  model.add_factor(fact2);
  model.add_factor(fact3);


  graph_type graph;
  
  // Make the graph
  construct_clique_graph(model, graph);

  // Make a shared data
  gl::thread_shared_data sdm;
  add_factors_to_sdm(sdm, model.factors());
  set_tree_sampler_constants(sdm,  1, 10, false, 0);

   
  
  // Make a simple engine
  gl::iengine* engine = make_engine(graph, sdm, 1, "fifo");
  assert(engine != NULL);
  
  // Add the root
  vertex_id_t root = 1;
  // Create a task to grow a tree from the root
  gl::update_task task(root, grow_root_update);
  engine->get_scheduler().add_task(task, grow_root_residual);
  
  engine->set_timeout(2);
  
  engine->start();





  delete engine;
}



int main(int argc, char** argv) {

  std::cout << "Testing factor graph inference code" << std::endl;
  
  test_toy_graph();

  return EXIT_SUCCESS;
}















