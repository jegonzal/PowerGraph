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
#ifndef GRAPHLAB_GRAPH_BUILTIN_PARSERS_HPP
#define GRAPHLAB_GRAPH_BUILTIN_PARSERS_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/logger/logger.hpp>
namespace graphlab {

  namespace builtin_parsers {
  
    template <typename Graph>
    bool snap_parser(Graph& graph, const std::string& srcfilename,
                     const std::string& str) {
      if (str.empty()) return true;
      else if (str[0] == '#') {
        std::cout << str << std::endl;
      }
      else {
        std::stringstream strm(str);
        size_t source, target;
        strm >> source >> target;
        if(source != target) graph.add_edge(source, target);
      }
      return true;
    }

    template <typename Graph>
    bool tsv_parser(Graph& graph, const std::string& srcfilename,
                    const std::string& str) {
      if (str.empty()) return true;
      std::stringstream strm(str);
      size_t source, target;
      strm >> source >> target;
      if(source != target) graph.add_edge(source, target);
      return true;
    }

    template <typename Graph>
    bool adj_parser(Graph& graph, const std::string& srcfilename,
                    const std::string& str) {
      if (str.empty()) return true;
      std::stringstream strm(str);
  
      size_t source, nneighbors;
      strm >> source >> nneighbors;
      if(!strm.good()) {
        logstream(LOG_ERROR) << "Adj format error on line: " << str << std::endl;
        return false; // failed to read the line
      }
      graph.add_vertex(source);
      for(size_t i = 0; i < nneighbors; ++i) {
        size_t target;
        strm >> target;
        if (!strm.good()) {
          logstream(LOG_ERROR) << "Adj format error on source vertex: " 
                               << source << std::endl;
          return false;
        }
        if(source != target) graph.add_edge(source, target);
      }
      return true;
    }


    template <typename Graph>
    struct tsv_writer{
      typedef typename Graph::vertex_type vertex_type;
      typedef typename Graph::edge_type edge_type;
      std::string save_vertex(vertex_type) { return ""; }
      std::string save_edge(edge_type e) {
        return tostr(e.source().id()) + "\t" + tostr(e.target().id()) + "\n";
      }
    };


  } // namespace builtin_parsers
} // namespace graphlab

#endif
