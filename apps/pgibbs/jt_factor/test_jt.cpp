#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <graphlab.hpp>

#include "data_structures.hpp"



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



int main(int argc, char** argv) {
  // test_jt_building(argc, argv);
  test_alchemy(argc, argv);

  return EXIT_SUCCESS;
} // end of main













#include <graphlab/macros_undef.hpp>

