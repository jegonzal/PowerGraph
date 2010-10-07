#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <graphlab.hpp>

#include "data_structures.hpp"
#include "sequential_jt_gibbs.hpp"


#include <graphlab/macros_def.hpp>


// void test_jt_building(int argc, char** argv) {
//   std::cout << "This is simple test code" << std::endl;


//   jt_sampler sampler;
//   size_t rows = 5;
//   size_t cols = 7;


    
//   // Make all variables
//   std::vector<variable_t> vars(rows * cols);
//   for(size_t i = 0; i < rows * cols; ++i) {
//     vars[i] = variable_t(i, 3);
//   }

//   // add all edges
//   for(size_t i = 0; i < rows; ++i) {
//     for(size_t j = 0; j < cols; ++j) {
//       size_t vertid = i * cols + j;
//       variable_t var = vars[vertid];
//       if(i + 1 < rows)
//         sampler.all_neighbors[var].insert(vars[(i+1)*cols+j]);
//       if(i - 1 < rows)
//         sampler.all_neighbors[var].insert(vars[(i-1)*cols+j]);
//       if(j + 1 < cols)
//         sampler.all_neighbors[var].insert(vars[i*cols+j+1]);
//       if(j - 1 < cols)
//         sampler.all_neighbors[var].insert(vars[i*cols+j-1]);      
//     }
//   }

//   // Print the tree
//   typedef std::set<variable_t> clique_type;
//   typedef std::pair<variable_t, clique_type> hood_type;

//   foreach(const hood_type& hood, sampler.all_neighbors) {
//     std::cout << hood.first << ": ";
//     foreach(variable_t var, hood.second) {
//       std::cout << var << " ";
//     }
//     std::cout << std::endl;
//   }

//   // Construct the factors 
//   sampler.in_tree.insert(vars.begin(), vars.end());


//   std::cout << std::endl << std::endl;

//   // Run the clique building
//   sampler.rebuild_cliques();

//   // foreach(const clique_type& clique, sampler.cliques) {    
//   //   foreach(variable_t var, clique) {
//   //     std::cout << var << " ";
//   //   }
//   //   std::cout << std::endl;
//   // }

// }


void test_alchemy(int argc, char** argv) {
  factorized_model model;
  model.load_alchemy(argv[1]);
  model.save_alchemy(argv[2]);
}


void test_compute_tree_width(int argc, char** argv) {
  factorized_model model;
  model.load_alchemy(argv[1]);
  
  // Build the maps
  vset_map var2factors;
  vset_map factor2vars;
  vertex_id_t fid = 0;
  foreach(const factor_t& factor, model.factors()){
    domain_t dom = factor.args();    
    std::cout << "Factor: " << fid << ": ";
    for(size_t i = 0; i < dom.num_vars(); ++i) {
      vertex_id_t varid = dom.var(i).id;
      var2factors[varid].insert(fid);
      factor2vars[fid].insert(varid);
      std::cout << varid << " ";
    }
    std::cout << std::endl;
    fid++;
  }

  size_t tree_width = compute_tree_width(var2factors, factor2vars);
  std::cout << "Tree Width: " << tree_width << std::endl;
}


void test_compute_tree_width2(int argc, char** argv) {
  // Build the maps
  vset_map var2factors;
  vset_map factor2vars;


  factor2vars[0].insert(0);
  factor2vars[0].insert(1);
  factor2vars[0].insert(2);
  factor2vars[0].insert(3);

  factor2vars[1].insert(1);
  factor2vars[1].insert(3);
  factor2vars[1].insert(4);
  factor2vars[1].insert(5);

  foreach(vset_map::value_type pair, factor2vars){
    vertex_id_t fid = pair.first;
    std::cout << "Factor: " << fid << ": ";
    foreach(vertex_id_t varid, pair.second) {
      std::cout << varid << " ";
      var2factors[varid].insert(fid);
    }
    std::cout << std::endl;   
  }

  size_t tree_width = compute_tree_width(var2factors, factor2vars);
  std::cout << "Tree Width: " << tree_width << std::endl;
}

void test_fast_set(int argc, char** argv) {
  const size_t SET_SIZE(10);
  typedef graphlab::fast_set<SET_SIZE, size_t> set_t;


  set_t set;
  set += 1;
  set += 7;
  set += 5;
  set += 2;
  set -= 4;
  set -= 2;
  set += 3;
  std::cout << set << std::endl;

  set_t set2 = set_t(3) +  9 + 4 + 5;
  std::cout << set2 << std::endl;

  std::cout << set - set2 << std::endl;
  std::cout << set * set2 << std::endl;

  foreach(size_t elem, (set + set2)) {
    std::cout << elem << ", ";
  }
  std::cout << std::endl;

}




int main(int argc, char** argv) {
  // test_jt_building(argc, argv);
  //  test_alchemy(argc, argv);
  // test_fast_set(argc, argv);
  test_compute_tree_width2(argc, argv);


  return EXIT_SUCCESS;
} // end of main













#include <graphlab/macros_undef.hpp>

