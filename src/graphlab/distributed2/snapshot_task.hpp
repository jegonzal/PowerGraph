#ifndef GRAPHLAB_SNAPSHOT_TASK_HPP
#define GRAPHLAB_SNAPSHOT_TASK_HPP

#include <graphlab/logger/assertions.hpp>
#include <graphlab/scope/iscope.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {
namespace gl_impl {

template <typename Graph, typename CallbackType>
void snapshot_update(iscope<Graph>& scope, 
                     icallback<Graph>& _callback) {
  CallbackType& callback = *dynamic_cast<CallbackType*>(&_callback);
  
  if (callback.get_snapshot_token(scope.vertex())) return;
  
  if (callback.vertex_modified_since_last_snapshot(scope.vertex())) {
    callback.save_vertex(scope.vertex(), scope.const_vertex_data());
  }
                       
  foreach(edge_id_t ineid, scope.in_edge_ids()) {
    if (callback.get_snapshot_token(scope.source(ineid)) == false) {
      if(callback.edge_modified_since_last_snapshot(ineid)) {
        // serialize this edge!
        callback.save_edge(ineid,
                           scope.source(ineid),
                           scope.vertex(), 
                           scope.const_edge_data(ineid));
      }
      update_task<Graph> task(scope.source(ineid), snapshot_update<Graph, CallbackType>);      
      callback.add_task(task, 100.0); 
    }
  }
                      
  foreach(edge_id_t outeid, scope.out_edge_ids()) {
    if (callback.get_snapshot_token(scope.target(outeid)) == false) {
      if(callback.edge_modified_since_last_snapshot(outeid)) {
        // serialize this edge!
        callback.save_edge(outeid, 
                           scope.vertex(), 
                           scope.target(outeid),
                           scope.const_edge_data(outeid));
      }
      update_task<Graph> task(scope.target(outeid), snapshot_update<Graph, CallbackType>);      
      callback.add_task(task, 100.0); 
    }
  }
  callback.set_and_synchronize_token(scope.vertex());
}

} // namespace gl_impl

} // namespace graphlab

#include <graphlab/macros_undef.hpp>
#endif
