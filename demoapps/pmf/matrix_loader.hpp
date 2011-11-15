/* Copyright (c) 2009 Carnegie Mellon University. 
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
 *  Code written by Danny Bickson, CMU
 *  Any changes to the code must include this original license notice in full.
 */


#ifndef MATRIX_LOADER_HPP
#define MATRIX_LOADER_HPP

#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>

template<typename Graph>
struct matrix_entry {
  typedef Graph graph_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  vertex_id_type source, target;
  edge_data_type edata;
}; // end of matrix_entry



template<typename Graph>
bool load_graph(const std::string& fname,
                const std::string& format,
                const double test_prop,
                Graph& graph,
                std::vector< matrix_entry<Graph> >& test_set) {

  return false;
} // end of load graph










#include <graphlab/macros_undef.hpp>
#endif // end of matrix_loader_hpp
