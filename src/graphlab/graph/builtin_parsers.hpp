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
#include <graphlab/serialization/serialization_includes.hpp>
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
        size_t source, target;
        char* targetptr;
        source = strtoul(str.c_str(), &targetptr, 10);
        if (targetptr == NULL) return false;
        target = strtoul(targetptr, NULL, 10);
        if(source != target) graph.add_edge(source, target);
      }
      return true;
    }

    template <typename Graph>
    bool tsv_parser(Graph& graph, const std::string& srcfilename,
                    const std::string& str) {
      if (str.empty()) return true;
      size_t source, target;
      char* targetptr;
      source = strtoul(str.c_str(), &targetptr, 10);
      if (targetptr == NULL) return false;
      target = strtoul(targetptr, NULL, 10);
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



    
    template <typename Graph>
    struct graphjrl_writer{
      typedef typename Graph::vertex_type vertex_type;
      typedef typename Graph::edge_type edge_type;

      /**
       * \internal
       * Replaces \255 with \255\1
       * Replaces \\n with \255\0
       */
      static std::string escape_newline(charstream& strm) {
        size_t ctr = 0;
        char *ptr = strm->c_str();
        size_t strmlen = strm->size();
        for (size_t i = 0;i < strmlen; ++i) {
          ctr += (ptr[i] == (char)255 || ptr[i] == '\n');
        }

        std::string ret(ctr + strmlen, 0);

        size_t target = 0;
        for (size_t i = 0;i < strmlen; ++i, ++ptr) {
          if ((*ptr) == (char)255) {
            ret[target] = (char)255;
            ret[target + 1] = 1;
            target += 2;
          }
          else if ((*ptr) == '\n') {
            ret[target] = (char)255;
            ret[target + 1] = 0;
            target += 2;
          }
          else {
            ret[target] = (*ptr);
            ++target;
          }
        }
        assert(target == ctr + strmlen);
        return ret;
      }

      /**
       * \internal
       * Replaces \255\1 with \255
       * Replaces \255\0 with \\n
       */
      static std::string unescape_newline(const std::string& str) {
        size_t ctr = 0;
        // count the number of escapes
        for (size_t i = 0;i < str.length(); ++i) {
          ctr += (str[i] == (char)255);
        }
        // real length is string length - escapes
        std::string ret(str.size() - ctr, 0);

        const char escapemap[2] = {'\n', (char)255};
        
        size_t target = 0;
        for (size_t i = 0;i < str.length(); ++i, ++target) {
          if (str[i] == (char)255) {
            // escape character
            ++i;
            ASSERT_MSG(str[i] == 0 || str[i] == 1,
                       "Malformed escape sequence in graphjrl file.");
            ret[target] = escapemap[(int)str[i]];
          }
          else {
            ret[target] = str[i];
          }
        }
        return ret;
      }
      
      std::string save_vertex(vertex_type v) {
        charstream strm(128);
        oarchive oarc(strm);
        oarc << char(0) << v.id() << v.data();
        strm.flush();
        return escape_newline(strm) + "\n";
      }
      
      std::string save_edge(edge_type e) {
        charstream strm(128);
        oarchive oarc(strm);
        oarc << char(1) << e.source().id() << e.target().id() << e.data();
        strm.flush();
        return escape_newline(strm) + "\n";
      }
    };



    template <typename Graph>
    bool graphjrl_parser(Graph& graph, const std::string& srcfilename,
                    const std::string& str) {
      std::string unescapedstr = graphjrl_writer<Graph>::unescape_newline(str);
      boost::iostreams::stream<boost::iostreams::array_source>
                      istrm(unescapedstr.c_str(), unescapedstr.length());
      iarchive iarc(istrm);
      
      char entrytype;
      iarc >> entrytype;
      if (entrytype == 0) {
        typename Graph::vertex_id_type vid;
        typename Graph::vertex_data_type vdata;
        iarc >> vid >> vdata;
        graph.add_vertex(vid, vdata);
      }
      else if (entrytype == 1) {
        typename Graph::vertex_id_type srcvid, destvid;
        typename Graph::edge_data_type edata;
        iarc >> srcvid >> destvid >> edata;
        graph.add_edge(srcvid, destvid, edata);
      }
      else {
        return false;
      }
      return true;
    }
    
  } // namespace builtin_parsers
} // namespace graphlab

#endif
