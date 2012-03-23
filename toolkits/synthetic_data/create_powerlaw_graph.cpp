/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#include <iostream>
#include <algorithm>
#include <fstream>
#include <set>
#include <boost/unordered_set.hpp>
#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>




int main(int argc, char** argv) {

  double alpha = 2.1;
  size_t nverts = 10;
  size_t fanout = 2;
  std::string fname = "graph.tsv";
  graphlab::command_line_options 
    clopts("Generate synthetic graph data.", true);
  clopts.attach_option("graph", &fname, fname,
                       "The name of the graph file."); 
  clopts.attach_option("alpha", &alpha, alpha, 
                       "Power law constant deg^{-a}.");
  clopts.attach_option("nverts", &nverts, nverts,
                       "Number of vertices.");
  clopts.attach_option("fanout", &fanout, fanout,
                       "The fanout of each page");

  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }



  std::vector<double> prob(nverts, 0);
  double Z = 0;
  for(size_t i = 0; i < nverts; ++i) Z += (prob[i] = std::pow(i+1,-alpha));
  // Normalize and convert to CDF
  for(size_t i = 0; i < nverts; ++i) {
    //    std::cout << prob[i] << '\t';
    prob[i] = prob[i]/Z + ((i>0)? prob[i-1] : 0);
    // std::cout << prob[i] << std::endl;
  }
 
  std::ofstream fout(fname.c_str());
  
  boost::unordered_set<size_t> targets;
  for(size_t source = 0; source < nverts; ++source) {
    targets.clear();
    while(targets.size() < fanout) {
      const size_t target =
        std::upper_bound(prob.begin(), prob.end(), 
                         graphlab::random::rand01()) - prob.begin();
      if(source != target) targets.insert(target);
    }
    foreach(size_t target, targets)
      fout << source << '\t' << target << '\n';
  }
  fout.close();
  
 
  return EXIT_SUCCESS;

} // end of main
