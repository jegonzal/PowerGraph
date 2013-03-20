#ifndef GRAPHLAB_GL3_TASK_HPP
#define GRAPHLAB_GL3_TASK_HPP
#include <vector>
#include <boost/function.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {
#define GL3_BROADCAST_TASK_ID 255

template <typename GraphType, typename EngineType>
struct gl3task_descriptor {
  virtual any exec(GraphType& graph,
                   vertex_id_type,
                   const any& params,
                   EngineType* engine,
                   std::vector<simple_spinlock>& vlocks,
                   std::vector<simple_spinlock>& elocks) = 0;
  virtual void combine(any& a, const any& b, const any& params) = 0;
};


/**************************************************************************/
/*                                                                        */
/*                    Defines the MapReduce Task Type                     */
/*                                                                        */
/**************************************************************************/



struct map_reduce_neighbors_task_param{
  bool in;
  bool out;
  void save(oarchive& oarc) const {
    oarc << in << out;
  }
  void load(iarchive& iarc) {
    iarc >> in >> out;
  }
};

template <typename GraphType, typename EngineType, typename T>
struct map_reduce_neighbors_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;
  typedef typename GraphType::vertex_type vertex_type;
  typedef typename GraphType::edge_type edge_type;

  typedef boost::function<T (const vertex_type&,
                             edge_type&,
                             const vertex_type&)> map_fn_type;
  typedef boost::function<void (T&, const T&)> combine_fn_type;



  map_fn_type map_fn;
  combine_fn_type combine_fn;
  T zero;

  map_reduce_neighbors_task_descriptor(map_fn_type mapper,
                                       combine_fn_type combiner,
                                       const T& zero)
      :map_fn(mapper), combine_fn(combiner), zero(zero) { }

  virtual void combine(any& a, const any& b, const any& params) {
    combine_fn(a.as<T>(), b.as<const T>());
  }

  any exec(GraphType& graph, vertex_id_type vid, const any& params,
           EngineType* engine,
           std::vector<simple_spinlock>& vlocks,
           std::vector<simple_spinlock>& elocks) {
    const map_reduce_neighbors_task_param& task_param = params.as<map_reduce_neighbors_task_param>();

    bool in = task_param.in;
    bool out = task_param.out;

    typedef typename GraphType::local_edge_type local_edge_type;
    typedef typename GraphType::local_vertex_type local_vertex_type;
    typedef typename GraphType::edge_type edge_type;
    typedef typename GraphType::vertex_type vertex_type;

    lvid_type lvid = graph.local_vid(vid);
    local_vertex_type lvertex = graph.l_vertex(lvid);
    vertex_type vertex = vertex_type(lvertex);
    T agg = zero;
    if (in) {
      foreach(local_edge_type ledge, lvertex.in_edges()) {
        edge_type edge(ledge);
        vertex_type other(ledge.source());
        elocks[ledge.id()].lock();
        combine_fn(agg, map_fn(vertex, edge, other));
        elocks[ledge.id()].unlock();
      }
    }

    if (out) {
      foreach(local_edge_type ledge, lvertex.out_edges()) {
        edge_type edge(ledge);
        vertex_type other(ledge.target());
        elocks[ledge.id()].lock();
        combine_fn(agg, map_fn(vertex, edge, other));
        elocks[ledge.id()].unlock();
      }
    }
    any ret(agg);
    return ret;
  }
};







/**************************************************************************/
/*                                                                        */
/*                    Defines the Broadcast Task Type                     */
/*                                                                        */
/**************************************************************************/



struct broadcast_task_param {
  bool in;
  bool out;
  any message;
  void save(oarchive& oarc) const {
    oarc << in << out << message;
  }
  void load(iarchive& iarc) {
    iarc >> in >> out >> message;
  }
};

template <typename GraphType, typename EngineType>
struct broadcast_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;

  broadcast_task_descriptor() { }

  virtual void combine(any& a, const any& b, const any& params) {
  }

  any exec(GraphType& graph, vertex_id_type vid, const any& params,
           EngineType* engine,
           std::vector<simple_spinlock>& vlocks,
           std::vector<simple_spinlock>& elocks) {
    const broadcast_task_param& task_param = params.as<broadcast_task_param>();
    typedef typename EngineType::message_type message_type;
    typedef typename GraphType::local_edge_type local_edge_type;
    typedef typename GraphType::local_vertex_type local_vertex_type;
    typedef typename GraphType::vertex_type vertex_type;

    message_type msg;
    if (!task_param.message.empty()) {
      msg = task_param.message.as<message_type>();
    }
    lvid_type lvid = graph.local_vid(vid);
    local_vertex_type lvertex = graph.l_vertex(lvid);
    if (task_param.in) {
      foreach(local_edge_type edge, lvertex.in_edges()) {
        engine->internal_signal(vertex_type(edge.source()));
      }
    }

    if (task_param.out) {
      foreach(local_edge_type edge, lvertex.out_edges()) {
        engine->internal_signal(vertex_type(edge.target()));
      }
    }
    return any();
  }
};





} // graphlab

#include <graphlab/macros_undef.hpp>
#endif
