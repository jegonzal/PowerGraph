#ifndef GRAPHLAB_GL3_CONTEXT_HPP
#define GRAPHLAB_GL3_CONTEXT_HPP
#include <boost/shared_ptr.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/engine/gl3task.hpp>
namespace graphlab {


/**
 * The context object which provides access to the registered tasks in
 * the gl2engine.
 */
template <typename EngineType>
struct gl3context {
  typedef typename EngineType::graph_type graph_type;

  typedef typename graph_type::local_vertex_type local_vertex_type;
  typedef typename graph_type::vertex_type vertex_type;
  typedef typename graph_type::vertex_data_type vertex_data_type;
  typedef typename graph_type::edge_type edge_type;
  typedef typename graph_type::local_edge_type local_edge_type;
  typedef typename graph_type::edge_data_type edge_data_type;
  typedef typename EngineType::message_type message_type;

  EngineType* engine;
  lvid_type lvid;
  size_t thread_id;

  /// returns the process ID of the current process
  procid_t procid() const {
    return engine->procid();
  }

  /// Returns a reference to the engine
  EngineType& get_engine() const {
    return (*engine);
  }

  /**
   * Perform a map reduce operation over neighbors.
   * The taskid must have been registered with the engine before utilization.
   * See the gl3engine::register_map_reduce() function.
   * edir can be either IN_EDGES, OUT_EDGES or ALL_EDGES.
   *
   * Can be called from within a vertex program or a parfor_all_local_vertices()
   */
  template <typename T>
  T map_reduce(size_t taskid,
               edge_dir_type edir) {
    ASSERT_NE(lvid, (lvid_type)(-1));
    // can only be called from within an vertex program.
    // thus spawn_subtask, which requires the lvid to be locked, is safe
    map_reduce_neighbors_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    any ret = engine->spawn_subtask(lvid, taskid, any(task_param));
    std::pair<T, bool>& retp = ret.as<std::pair<T, bool> >();
    return retp.first;
  }

  /**
   * Broadcasts a signal (schedules) to a set of neighbors.
   * edir can be either IN_EDGES, OUT_EDGES or ALL_EDGES.
   *
   * Can be called from within a vertex program or a parfor_all_local_vertices().
   */
  void broadcast_signal(edge_dir_type edir,
                        const message_type& msg = message_type()) {
    ASSERT_NE(lvid, (lvid_type)(-1));
    // can only be called from within an vertex program.
    // thus spawn_subtask, which requires the lvid to be locked, is safe
    broadcast_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    task_param.message = msg;
    engine->spawn_subtask(lvid, GL3_BROADCAST_TASK_ID, any(task_param), true);
  }

  /**
   * Issues an edge transform operation on a set of neighboring edges.
   * The taskid must have been registered with the engine before utilization.
   * See the gl3engine::register_edge_transform() function.
   * edir can be either IN_EDGES, OUT_EDGES or ALL_EDGES.
   *
   * Can be called from within a vertex program or a parfor_all_local_vertices().
   */
  void edge_transform(size_t taskid,
                      edge_dir_type edir,
                      bool wait = true) {
    ASSERT_NE(lvid, (lvid_type)(-1));
    // can only be called from within an vertex program.
    // thus spawn_subtask, which requires the lvid to be locked, is safe
    edge_transform_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    engine->spawn_subtask(lvid, taskid, any(task_param), wait);
  }

  /**
   * Gets a requested set of keys from the DHT.
   * The resultant map will contain all the requested keys. Note that
   * the values may be empty. it is important to check with any::empty()
   * before accessing.
   *
   * Can be called anywhere.
   */
  boost::unordered_map<size_t, any>
      dht_gather(const std::vector<size_t>& entries) {
    boost::unordered_map<size_t, any> ret;
    if (entries.size() == 0) return ret;

    size_t nprocs = engine->rmi.numprocs();
    // construct all the gather entries
    std::vector<dht_gather_task_param> params(nprocs);
    for (size_t i = 0;i < entries.size(); ++i) {
      params[entries[i] % nprocs].gather_entries.push_back(entries[i]);
    }

    std::vector<procid_t> targetprocs;
    std::vector<any> targetparam;

    for (size_t i = 0;i < params.size(); ++i) {
      if (!params[i].gather_entries.empty()) {
        targetprocs.push_back(i);
        targetparam.push_back(any(params[i]));
      }
    }

    any res = engine->spawn_subtask(GL3_DHT_GATHER_TASK_ID,
                                 targetprocs,
                                 targetparam);

    if (res.empty()) return ret;

    std::vector<std::pair<size_t, any> >& retvec =
        res.as<std::vector<std::pair<size_t, any> > >();

    for (size_t i = 0;i < retvec.size(); ++i) {
      ret[retvec[i].first].swap(retvec[i].second);
    }
    return ret;
  }

  /**
   * Issues a DHT scatter operation. The operation must be registered with
   * gl3engine::register_dht_scatter() before usage.
   * See gl3engine::register_dht_scatter() for details.
   *
   * Can be called anywhere.
   */
  void dht_scatter(size_t taskid,
                   const boost::unordered_map<size_t, any>& entries) {
    if (entries.size() == 0) return;

    size_t nprocs = engine->rmi.numprocs();
    // construct all the scatter entries
    std::vector<dht_scatter_task_param> params(nprocs);
    boost::unordered_map<size_t, any>::const_iterator iter = entries.begin();
    while(iter != entries.end()) {
      params[iter->first % nprocs].scatter_entries.push_back(
          std::pair<size_t, any>(iter->first, iter->second));
      ++iter;
    }

    std::vector<procid_t> targetprocs;
    std::vector<any> targetparam;

    for (size_t i = 0;i < params.size(); ++i) {
      if (!params[i].scatter_entries.empty()) {
        targetprocs.push_back(i);
        targetparam.push_back(any(params[i]));
      }
    }

    any res = engine->spawn_subtask(taskid,
                                 targetprocs,
                                 targetparam,
                                 true);
  }

  /**
   * Issues an vertex delta operation on a particular vertex.
   * The taskid must have been registered with the engine before utilization.
   * See gl3engine::register_vertex_delta() for details.
   *
   * Deltas are immediately applied locally, as well as forwarded to the master.
   * However, changes to the master are not synchronized. You will need to
   * provide additional synchronization routines.
   *
   * Can be called from anywhere as long as a vertex_type is available.
   */
  template <typename T>
  void send_delta(size_t taskid,
                  const vertex_type& vtx,
                  T& delta_val) {
    vertex_delta_task_param param;
    param.delta_value = delta_val;
    engine->spawn_subtask_to_master(vtx.local_id(),
                                    taskid,
                                    any(param),
                                    true);
    // I want to apply the delta immediately locally as well.
    // this is ... a little intrusive.
    local_vertex_type lvtx(vtx);
    if (!lvtx.owned()) {
      lvid_type lvid = lvtx.id();
      engine->unlock_vlock_for_subtask_execution(lvid);
      any ret = engine->task_types[taskid]->exec(engine->graph,
                                                 lvtx.global_id(),
                                                 any(param),
                                                 engine);
      engine->relock_vlock_for_subtask_execution(lvid);
    }
  }

  /**
   * If the current vertex is a master vertex, synchronizes the data on the
   * master vertex with all mirrors. Has no effect if it is a slave vertex.
   *
   * Can be called from anywhere as long as a vertex type is avaiable.
   */
  void synchronize(const vertex_type& vtx) {
    local_vertex_type lvtx(vtx);
    if (lvtx.owned()) {
      engine->spawn_subtask(vtx.local_id(),
                            GL3_SYNCHRONIZE_MIRRORS_TASK_ID,
                            any(),
                            true);
    }
  }

};


}; // namespace graphlab

#undef FRESULT
#undef REMOVE_CONST_REF

#endif

