/*  
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


#include <graphlab.hpp>
#include <fstream>

int main(int argc, char** argv) {
  std::ofstream fout("synthetic.txt");
  
  size_t num_clusters = 2;
  size_t dim = 2;
  size_t ndata = 10000;
  if (argc >= 2) num_clusters = atoi(argv[1]);
  if (argc >= 3) dim = atoi(argv[2]);
  if (argc >= 4) ndata = atoi(argv[3]);

  std::cout << "Usage: generate_synthetic [NClusters] [Dimensions] [Ndata]\n";

  std::vector< std::vector<double> > centers(num_clusters);
  for (size_t i = 0;i < centers.size(); ++i) {
    std::cout << "Center " << i << " at: " ; 
    for (size_t j = 0; j < dim; ++j) {
      double r = graphlab::random::fast_uniform<double>(-10,10);
      std::cout << r << "\t";
      centers[i].push_back(r);
    }
    std::cout << "\n";
  }

  // now generate data points
  // 
  for (size_t i = 0;i < ndata; ++i) {
    size_t c = graphlab::random::fast_uniform<size_t>(0, centers.size() - 1);
    for (size_t j = 0; j < dim; ++j) {
      double d = graphlab::random::gaussian() + centers[c][j];
      fout << d << "\t";
    }
    fout << "\n";
  }
}
