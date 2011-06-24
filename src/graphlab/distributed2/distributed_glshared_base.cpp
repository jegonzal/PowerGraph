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


#include <vector>
#include <graphlab/distributed2/distributed_glshared.hpp>

namespace graphlab {
namespace distgl_impl {
	

std::vector<distributed_glshared_base*>& get_global_dist_glshared_registry() {
  static std::vector<distributed_glshared_base*> global_registry;
  return global_registry;
}	

void register_dist_glshared(distributed_glshared_base* glsharedobj) {
  get_global_dist_glshared_registry().push_back(glsharedobj);
}



} // namespace distgl_impl
} // namespace graphlab

