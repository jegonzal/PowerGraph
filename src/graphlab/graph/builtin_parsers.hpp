#include <string>
#include <sstream>
#include <iostream>
#include <graphlab/logger/logger.hpp>
namespace graphlab {

namespace builtin_parsers {
  
template <typename VertexType, typename EdgeType>
bool snap_parser(graphlab::distributed_graph<VertexType, EdgeType>& graph,
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


template <typename VertexType, typename EdgeType>
bool tsv_parser(graphlab::distributed_graph<VertexType, EdgeType>& graph,
                const std::string& srcfilename,
                const std::string& str) {
  if (str.empty()) return true;
  std::stringstream strm(str);
  size_t source, target;
  strm >> source >> target;
  if(source != target) graph.add_edge(source, target);
  return true;
}

template <typename VertexType, typename EdgeType>
bool adj_parser(graphlab::distributed_graph<VertexType, EdgeType>& graph,
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

} // namespace builtin_callbacks
} // namespace graphlab