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


#include <graphlab/factors/binary_factor.hpp>

std::ostream& operator<<(std::ostream& out, 
                         const graphlab::binary_factor& fact) {
  out << "Binary Factor(v_" << fact.var1() << " in {1..."
      << fact.arity1() << "}, " 
      << ", v_ " << fact.var2() << " in {1..." 
      << fact.arity2() << "})" << std::endl;
  for(uint16_t i = 0; i < fact.arity1(); ++i) {
    for(uint16_t j = 0; j < fact.arity2(); ++j) {
      out << fact.logP(i,j) << " ";
    }
    out << std::endl;
  }
  return out;
} // end of operator<<

