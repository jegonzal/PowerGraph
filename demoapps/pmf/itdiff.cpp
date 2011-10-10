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

#include "mathlayer.hpp"
#include "graphlab.hpp"
#include <cassert>
#include <iostream>

int main(int argc, char** argv) {
  assert(argc == 3);
  
  mat U1, U2, V1, V2;
  it_file input(argv[1]);
  input >> Name("User");
  input  >> U1;
  input >> Name("Movie");
  input  >> V1;
  it_file input2(argv[2]);
  input2 >> Name("User");
  input2  >> U2;
  input2 >> Name("Movie");
  input2 >> V2;
  
  assert(U1.cols() == U2.cols());
  assert(U1.rows() == U2.rows());
  assert(V1.cols() == V2.cols());
  assert(V1.rows() == V2.rows());

  mat diff = U1-U2;  
  if (abs_sum(diff) > U1.rows() * U1.cols() * 1E-10)
    logstream(LOG_WARNING) << "Itdiff got U diference of: " << abs_sum(diff) << std::endl;
  diff = V1-V2;
  if (abs_sum(diff) > V1.rows() * V1.cols() * 1E-10)
    logstream(LOG_WARNING) << "Itdiff got V d0iference of: " << abs_sum(diff) << std::endl;
  return 0;
}
