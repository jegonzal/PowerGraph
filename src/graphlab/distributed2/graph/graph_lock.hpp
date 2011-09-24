#ifndef GRAPHLAB_GRAPH_LOCK_INTERFACE_HPP
#define GRAPHLAB_GRAPH_LOCK_INTERFACE_HPP

#include <boost/function.hpp>
#include <graphlab/scope/iscope.hpp>

namespace graphlab {

template <typename GraphType>
class graph_lock {
 public:
  typedef typename GraphType::vertex_id_type vertex_id_type;  
  virtual void scope_request(vertex_id_type globalvid,
                           boost::function<void(vertex_id_type)> handler,
                           scope_range::scope_range_enum scopetype) = 0;

  virtual void scope_unlock(vertex_id_type globalvid,
                           scope_range::scope_range_enum scopetype) = 0;
                           
  virtual void print_state() = 0; // for debugging                           
};

}
#endif
