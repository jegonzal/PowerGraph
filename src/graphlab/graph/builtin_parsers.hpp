#ifndef GRAPHLAB_GRAPH_BUILTIN_PARSERS_HPP
#define GRAPHLAB_GRAPH_BUILTIN_PARSERS_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/logger/logger.hpp>
namespace graphlab {

namespace builtin_parsers {
  
template <typename VertexDataType, typename EdgeDataType>
bool snap_parser(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                 const std::string& srcfilename,
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


template <typename VertexDataType, typename EdgeDataType>
bool tsv_parser(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                const std::string& srcfilename,
                const std::string& str) {
  if (str.empty()) return true;
  std::stringstream strm(str);
  size_t source, target;
  strm >> source >> target;
  if(source != target) graph.add_edge(source, target);
  return true;
}

template <typename VertexDataType, typename EdgeDataType>
bool adj_parser(graphlab::distributed_graph<VertexDataType, EdgeDataType>& graph,
                const std::string& srcfilename,
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


template <typename VertexDataType, typename EdgeDataType>
struct tsv_writer{
  typedef typename graphlab::distributed_graph<VertexDataType, EdgeDataType>::vertex_type vertex_type;
  typedef typename graphlab::distributed_graph<VertexDataType, EdgeDataType>::edge_type edge_type;

  std::string save_vertex(vertex_type) { return ""; }
  std::string save_edge(edge_type e) {
    return tostr(e.source().id()) + "\t" + tostr(e.target().id()) + "\n";
  }
};


} // namespace builtin_parsers
} // namespace graphlab

#endif