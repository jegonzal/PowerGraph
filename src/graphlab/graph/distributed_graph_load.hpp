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

#ifndef GRAPHLAB_DISTRIBUTED_GRAPH_LOAD_HPP
#define GRAPHLAB_DISTRIBUTED_GRAPH_LOAD_HPP
#include <boost/function.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/builtin_parsers.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/fs_util.hpp>
#include <graphlab/util/hdfs.hpp>
namespace graphlab {

  namespace graph_ops {
   

    template<typename VertexDataType, typename EdgeDataType, typename Fstream>
    bool load_from_stream(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                          const std::string& srcfilename,
                          Fstream& fin, 
                          boost::function<bool(graphlab::distributed_graph<VertexDataType, EdgeDataType>&,
                                               const std::string&,
                                               const std::string&)> callback) {
      size_t linecount = 0;
      timer ti; ti.start();
      while(fin.good() && !fin.eof()) {
        std::string line;
        std::getline(fin, line);
        if (!callback(graph, srcfilename, line)) {
          logstream(LOG_WARNING) 
            << "Error parsing line " << linecount << " in "
            << srcfilename << ": " << std::endl
            << "\t\"" << line << "\"" << std::endl;  
          return false;
        }
        ++linecount;      
        if (ti.current_time() > 5.0) {
          logstream(LOG_INFO) << linecount << " Lines read" << std::endl;
          ti.start();
        }
      }
      return true;
    } // end of load from stream



    template <typename VertexDataType, typename EdgeDataType>
    void load_from_posixfs(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                           std::string path, 
                           boost::function<bool(graphlab::distributed_graph<VertexDataType, EdgeDataType>&,
                                                const std::string&,
                                                const std::string&)> callback) {
      // force a "/" at the end of the path
      // make sure to check that the path is non-empty. (you do not
      // want to make the empty path "" the root path "/" )
      if (path.length() > 0 && path[path.length() - 1] != '/') path = path + "/";
    
      std::vector<std::string> graph_files;
      graphlab::fs_util::list_files_with_prefix(path, "", graph_files);
      for(size_t i = 0; i < graph_files.size(); ++i) {
        graph_files[i] = path + graph_files[i];
      }
    
      for(size_t i = 0; i < graph_files.size(); ++i) {
        if (i % graph.numprocs() == graph.procid()) {
          std::cout << "Loading graph from file: " << graph_files[i] << std::endl;
          // is it a gzip file ?
          const bool gzip = boost::ends_with(graph_files[i], ".gz");
          // open the stream
          std::ifstream in_file(graph_files[i].c_str(), 
                                std::ios_base::in | std::ios_base::binary);
          // attach gzip if the file is gzip
          boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
          // Using gzip filter
          if (gzip) fin.push(boost::iostreams::gzip_decompressor());
          fin.push(in_file);
          const bool success = 
            graphlab::graph_ops::load_from_stream(graph, graph_files[i], fin, callback);
          ASSERT_TRUE(success);
          fin.pop();
          if (gzip) fin.pop();
        }
      }
    } // end of load from posixfs



    template <typename VertexDataType, typename EdgeDataType>
    void load_from_hdfs(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                        std::string path, 
                        boost::function<bool(graphlab::distributed_graph<VertexDataType, EdgeDataType>&,
                                             const std::string&,
                                             const std::string&)> callback) {
      // force a "/" at the end of the path
      // make sure to check that the path is non-empty. (you do not
      // want to make the empty path "" the root path "/" )
      if (path.length() > 0 && path[path.length() - 1] != '/') path = path + "/";
    
      ASSERT_TRUE(hdfs::has_hadoop());
      hdfs& hdfs = hdfs::get_hdfs();
    
      std::vector<std::string> graph_files;
      graph_files = hdfs.list_files(path);
    
      for(size_t i = 0; i < graph_files.size(); ++i) {
        if (i % graph.numprocs() == graph.procid()) {
          std::cout << "Loading graph from file: " << graph_files[i] << std::endl;
          // is it a gzip file ?
          const bool gzip = boost::ends_with(graph_files[i], ".gz");
          // open the stream
          graphlab::hdfs::fstream in_file(hdfs, graph_files[i]);
          boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
          fin.set_auto_close(false);
          if(gzip) fin.push(boost::iostreams::gzip_decompressor());
          fin.push(in_file);      
          const bool success = 
            graphlab::graph_ops::load_from_stream(graph, graph_files[i], fin, callback);
          ASSERT_TRUE(success);
          fin.pop();
          if (gzip) fin.pop();
        }
      }
    } // end of load from hdfs


  
    template <typename VertexDataType, typename EdgeDataType>
    void load(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
              std::string path, 
              boost::function<bool(graphlab::distributed_graph<VertexDataType, EdgeDataType>&,
                                   const std::string&,
                                   const std::string&)> callback) {
      graph.dc().full_barrier();
      if(boost::starts_with(path, "hdfs://")) {
        load_from_hdfs(graph, path, callback);
      } else {
        load_from_posixfs(graph, path, callback);
      }
      graph.dc().full_barrier();
    } // end of load
  

    // Synthetic Generators ===================================================>
    
    template <typename VertexDataType, typename EdgeDataType>
    void load_synthetic_powerlaw(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                                 size_t nverts, bool in_degree = false,
                                 double alpha = 2.1, size_t truncate = (size_t)(-1)) {
      graph.dc().full_barrier(); 
      std::vector<double> prob(std::min(nverts, truncate), 0);
      std::cout << "constructing pdf" << std::endl;
      for(size_t i = 0; i < prob.size(); ++i)
        prob[i] = std::pow(double(i+1), -alpha);
      std::cout << "constructing cdf" << std::endl;
      random::pdf2cdf(prob);
      std::cout << "Building graph" << std::endl;
      size_t target_index = graph.procid();
      size_t addedvtx = 0;

      for(size_t source = graph.procid(); source < nverts;
          source += graph.numprocs()) {
        const size_t out_degree = random::sample(prob) + 1;
        for(size_t i = 0; i < out_degree; ++i) {
          target_index = (target_index + 2654435761)  % nverts;
          if(source == target_index) {
            target_index = (target_index + 2654435761)  % nverts;
          }
          if(in_degree) graph.add_edge(target_index, source);
          else graph.add_edge(source, target_index);
        }
        ++addedvtx;
        if (addedvtx % 10000000 == 0) {
          std::cout << addedvtx << " inserted\n";
        }
      }
      graph.dc().full_barrier(); 
    } // end of load random powerlaw

    template <typename VertexDataType, typename EdgeDataType>
    void load(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
              std::string path, 
              std::string format) {
      boost::function<bool(graphlab::distributed_graph<VertexDataType, EdgeDataType>&,
                           const std::string&,
                           const std::string&)> callback;
      if (format == "snap") {
        callback = builtin_parsers::snap_parser<VertexDataType, EdgeDataType>;
      } else if (format == "adj") {
        callback = builtin_parsers::adj_parser<VertexDataType, EdgeDataType>;
      } else if (format == "tsv") {
        callback = builtin_parsers::tsv_parser<VertexDataType, EdgeDataType>;
      } else {
        logstream(LOG_ERROR)
          << "Unrecognized Format \"" << format << "\"!" << std::endl;
        return;
      }
      load(graph, path, callback);
    } // end of load



  } // namespace graph_ops

} // namespace graphlab
#endif
