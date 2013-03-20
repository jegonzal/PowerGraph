#ifndef GRAPHLAB_GL3_CONTEXT_HPP
#define GRAPHLAB_GL3_CONTEXT_HPP
#include <boost/shared_ptr.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/engine/gl3task.hpp>
#include <graphlab/parallel/qthread_future.hpp>
namespace graphlab {


template <typename EngineType>
struct gl3context {
  typedef typename EngineType::graph_type graph_type;
  typedef typename graph_type::vertex_data_type vertex_data_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef typename EngineType::message_type message_type;

  EngineType* engine;

  lvid_type lvid;

  any map_reduce(size_t taskid,
                 edge_dir_type edir) {
    map_reduce_neighbors_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    return engine->spawn_task(lvid, taskid, any(task_param));
  }


  void broadcast_signal(edge_dir_type edir,
                        const message_type& msg = message_type()) {
    broadcast_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    task_param.message = msg;
    engine->spawn_task(lvid, GL3_BROADCAST_TASK_ID, any(task_param), true);
  }

};


}; // namespace graphlab

#undef FRESULT
#undef REMOVE_CONST_REF

#endif

