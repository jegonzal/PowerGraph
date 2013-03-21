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
#define GL3_DHT_GATHER_TASK_ID 254
#define GL3_DHT_SCATTER_TASK_ID 253

template <typename GraphType, typename EngineType>
struct gl3task_descriptor {
  virtual any exec(GraphType& graph,
                   vertex_id_type,
                   const any& params,
                   EngineType* engine) = 0;
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
           EngineType* engine) {
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
        engine->elocks[ledge.id()].lock();
        combine_fn(agg, map_fn(vertex, edge, other));
        engine->elocks[ledge.id()].unlock();
      }
    }

    if (out) {
      foreach(local_edge_type ledge, lvertex.out_edges()) {
        edge_type edge(ledge);
        vertex_type other(ledge.target());
        engine->elocks[ledge.id()].lock();
        combine_fn(agg, map_fn(vertex, edge, other));
        engine->elocks[ledge.id()].unlock();
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
           EngineType* engine) {
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


/**************************************************************************/
/*                                                                        */
/*                    Defines the DHT Gather Task Type                    */
/*                                                                        */
/**************************************************************************/

struct dht_gather_task_param {
  std::vector<size_t> gather_entries;
  void save(oarchive& oarc) const {
    oarc << gather_entries;
  }
  void load(iarchive& iarc) {
    iarc >> gather_entries;
  }
};



template <typename GraphType, typename EngineType>
struct dht_gather_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef std::vector<std::pair<size_t, any> > gather_type;
  any exec(GraphType& graph,
           vertex_id_type,
           const any& params,
           EngineType* engine) {
    const dht_gather_task_param& task_param = params.as<dht_gather_task_param>();
    gather_type res;
    for (size_t i = 0;i < task_param.gather_entries.size(); ++i) {
      size_t entryid = task_param.gather_entries[i];
      // now, the next 8 bits
      size_t table = (entryid >> 8) % engine->NUM_DHTS;
      engine->dht_lock[table].lock();
      res.push_back(std::pair<size_t, any>(entryid, engine->dht[table][entryid]));
      engine->dht_lock[table].unlock();
    }
    any anyres(res);
    return anyres;
  }

  virtual void combine(any& a, const any& b, const any& params) {
    gather_type& aval = a.as<gather_type>();
    const gather_type& bval = b.as<const gather_type>();
    for (size_t i = 0;i < bval.size(); ++i) {
      aval.push_back(bval[i]);
    }
  }
};



/**************************************************************************/
/*                                                                        */
/*                    Defines the DHT ScatterTask Type                    */
/*                                                                        */
/**************************************************************************/

struct dht_scatter_task_param {
  std::vector<std::pair<size_t, any> > scatter_entries;
 void save(oarchive& oarc) const {
    oarc << scatter_entries;
  }
  void load(iarchive& iarc) {
    iarc >> scatter_entries;
  }
};



template <typename GraphType, typename EngineType>
struct dht_scatter_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef boost::function<void(any&, const any&)> scatter_fn_type;
  scatter_fn_type scatter_fn;

  any exec(GraphType& graph,
           vertex_id_type,
           const any& params,
           EngineType* engine) {
    const dht_scatter_task_param& task_param = params.as<dht_scatter_task_param>();
    for (size_t i = 0;i < task_param.scatter_entries.size(); ++i) {
      size_t entryid = task_param.scatter_entries[i].first;
      // now, the next 8 bits
      size_t table = (entryid >> 8) % engine->NUM_DHTS;
      engine->dht_lock[table].lock();
      scatter_fn(engine->dht[table][entryid], task_param.scatter_entries[i].second);
      engine->dht_lock[table].unlock();
    }
    return any();
  }



  virtual void combine(any& a, const any& b, const any& params) {
  }
};







} // graphlab

#include <graphlab/macros_undef.hpp>
#endif
