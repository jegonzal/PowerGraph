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
  size_t thread_id;

  void poll() {
    engine->poll(thread_id);
  }

  procid_t procid() const {
    return engine->procid();
  }

  EngineType& get_engine() const {
    return (*engine);
  }

  any map_reduce(size_t taskid,
                 edge_dir_type edir) {
    ASSERT_NE(lvid, (lvid_type)(-1));
    map_reduce_neighbors_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    return engine->spawn_task(lvid, taskid, any(task_param));
  }


  void broadcast_signal(edge_dir_type edir,
                        const message_type& msg = message_type()) {
    ASSERT_NE(lvid, (lvid_type)(-1));
    broadcast_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    task_param.message = msg;
    engine->spawn_task(lvid, GL3_BROADCAST_TASK_ID, any(task_param), true);
  }

  void edge_transform(size_t taskid,
                      edge_dir_type edir,
                      bool wait = true) {
    ASSERT_NE(lvid, (lvid_type)(-1));
    edge_transform_task_param task_param;
    task_param.in = (edir == IN_EDGES) || (edir == ALL_EDGES);
    task_param.out = (edir == OUT_EDGES) || (edir == ALL_EDGES);
    engine->spawn_task(lvid, taskid, any(task_param), wait);
  }


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

    any res = engine->spawn_task(GL3_DHT_GATHER_TASK_ID,
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

    any res = engine->spawn_task(taskid,
                                 targetprocs,
                                 targetparam,
                                 true);
  }



};


}; // namespace graphlab

#undef FRESULT
#undef REMOVE_CONST_REF

#endif

