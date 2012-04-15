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
 
/**
 * @file dispatcher.hpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#ifndef GRAPHLAB_DISPATCHER_HPP
#define GRAPHLAB_DISPATCHER_HPP

#include <graphlab.hpp>
#include "process.hpp"

struct vertex_data {};
struct edge_data {};
typedef graphlab::graph<vertex_data, edge_data> graph_type;

//////////////////////////////// UPDATE FUNCTOR ////////////////////////////////

/** dispatcher update function */
class dispatcher_update : 
  public graphlab::iupdate_functor<graph_type, dispatcher_update> {
private:
public:
  dispatcher_update();
  dispatcher_update(const dispatcher_update& other);
  inline void operator+=(const dispatcher_update& other);
  void operator()(icontext_type& context);
}; // end of dispatcher update functor

////////////////////////////////////////////////////////////////////////////////

#endif /* #ifndef GRAPHLAB_DISPATCHER_HPP */
