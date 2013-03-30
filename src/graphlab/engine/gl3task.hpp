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
#define GL3_SYNCHRONIZE_MIRRORS_TASK_ID 253

/*
 * The gl3task_descriptor is the general interface for a GraphLab 3 subtask.
 * Launching subtasks must be done through the gl3engine::spawn_subtask()
 * or gl3engine::spawn_subtask_to_master() functions. These functions
 * issue calls to remote (or local) machines to run a particular subtask
 * object. When a subtask is run, the exec() function is called.
 *
 * Return values are then returned to the sender, where they are merged using
 * the combine() function.
 *
 * Basically every task type must be a friend of gl3engine.
 * Most task types require an implementation in three places.
 * First is the task logic itself. Second is the launching logic in gl3context,
 * third is the task_id registration in GL3 engine. This ought to be refactored
 * so that the logic is contained in a smaller number of places.
 */
template <typename GraphType, typename EngineType>
struct gl3task_descriptor {
  virtual any exec(GraphType& graph,
                   vertex_id_type,
                   const any& params,
                   EngineType* engine) = 0;
  virtual void combine(any& a, const any& b, const any& params) = 0;
  // if returns true, this task will complete inside the RPC call itself
  // rather than placing it in a scheduler.
  virtual bool fast_path(const any& params) { return false; }


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

/*
 * Defines the map reduce over neighborhood task.
 */
template <typename GraphType, typename EngineType, typename T>
struct map_reduce_neighbors_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;
  typedef typename GraphType::vertex_type vertex_type;
  typedef typename GraphType::edge_type edge_type;
  typedef gl3task_descriptor<GraphType, EngineType> base_type;

  typedef boost::function<T (const vertex_type&,
                             edge_type&,
                             const vertex_type&)> map_fn_type;
  typedef boost::function<void (T&, const T&)> combine_fn_type;


  /*
   * This mess here is to permit two different map function types.
   * First is a map function which just takes in the neighbor vertex,
   * and the other which takes in the edge and both vertex endpoints.
   */
  struct create_map_reduce_task_impl{
    // The more general map function type
    typedef boost::function<T (const vertex_type&,
                               edge_type&,
                               const vertex_type&)>  full_map_function_type;

    // the simpler map function type
    typedef boost::function<T (const vertex_type&)> basic_map_function_type;


    /*
     * This function is enabled if the map function type matches the
     * general map function type. In which case, a task descriptor is created
     * using the user provided map function type.
     */
    template <typename MapFn>
    static typename boost::enable_if_c<boost::is_same<MapFn, full_map_function_type>::value,
                    base_type*>::type
      create(full_map_function_type mapfn,
             boost::function<void (T&, const T&)> combinefn) {
      return new map_reduce_neighbors_task_descriptor<GraphType, EngineType, T>(mapfn, combinefn);
    }

    // Converts a complex map function type to a simpler map function type
    static T simple_map_function_dispatch(boost::function<T (const vertex_type&)> mapfn,
                                   const vertex_type&,
                                   edge_type&,
                                   const vertex_type& other) {
      return mapfn(other);
    }

    /*
     * This function is enabled if the map function type matches the
     * simpler map function type. In which case, a task descriptor is created
     * by rebinding the simpler map function to simple_map_function_dispatch.
     */
    template <typename MapFn>
    static typename boost::enable_if_c<boost::is_same<MapFn, basic_map_function_type>::value,
                    base_type*>::type
      create(basic_map_function_type mapfn,
             boost::function<void (T&, const T&)> combinefn) {
      return new map_reduce_neighbors_task_descriptor<GraphType, EngineType, T>(
          boost::bind(simple_map_function_dispatch, mapfn, _1, _2, _3), combinefn);
    }
  };



  map_fn_type map_fn;
  combine_fn_type combine_fn;

  map_reduce_neighbors_task_descriptor(map_fn_type mapper,
                                       combine_fn_type combiner)
      :map_fn(mapper), combine_fn(combiner) { }

  /*
   * Combine results
   */
  virtual void combine(any& a, const any& b, const any& params) {
    std::pair<T, bool>& ap = a.as<std::pair<T, bool> >();
    const std::pair<T, bool>& bp = b.as<std::pair<T, bool> >();
    // if both are filled, we can just combine
    // if first is not filled, the result is the second
    if (ap.second && bp.second) combine_fn(ap.first, bp.first);
    else if (!ap.second) ap = bp;
  }

  /* Runs the map function
   * on the selected subset of edges.
   */
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
    bool filled = false;
    T agg = T();
    if (in) {
      foreach(local_edge_type ledge, lvertex.in_edges()) {
        edge_type edge(ledge);
        vertex_type other(ledge.source());
        engine->elocks[ledge.id()].lock();
        if (filled) {
          combine_fn(agg, map_fn(vertex, edge, other));
        } else {
          agg = map_fn(vertex, edge, other);
          filled = true;
        }
        engine->elocks[ledge.id()].unlock();
      }
    }

    if (out) {
      foreach(local_edge_type ledge, lvertex.out_edges()) {
        edge_type edge(ledge);
        vertex_type other(ledge.target());
        engine->elocks[ledge.id()].lock();
        if (filled) {
          combine_fn(agg, map_fn(vertex, edge, other));
        } else {
          agg = map_fn(vertex, edge, other);
          filled = true;
        }
        engine->elocks[ledge.id()].unlock();
      }
    }
    any ret(std::pair<T, bool>(agg, filled));
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

/*
 * Schedules a selected of neighbors with a particular message.
 */
template <typename GraphType, typename EngineType>
struct broadcast_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;

  broadcast_task_descriptor() { }

  // No combine needed. This function does not return
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
/*             Defines the Edge Transform Task Type                       */
/*                                                                        */
/**************************************************************************/


struct edge_transform_task_param {
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

/*
 * Transforms neigbors.
 */
template <typename GraphType, typename EngineType>
struct edge_transform_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;
  typedef typename GraphType::vertex_type vertex_type;
  typedef typename GraphType::edge_type edge_type;

  typedef boost::function<void (const vertex_type&,
                                edge_type&,
                                const vertex_type&)> edge_transform_fn_type;

  edge_transform_fn_type edge_transform_fn;

  edge_transform_task_descriptor(edge_transform_fn_type fn):edge_transform_fn(fn) { }

  virtual void combine(any& a, const any& b, const any& params) { }

  any exec(GraphType& graph, vertex_id_type vid, const any& params,
           EngineType* engine) {
    const edge_transform_task_param& task_param = params.as<edge_transform_task_param>();

    bool in = task_param.in;
    bool out = task_param.out;

    typedef typename GraphType::local_edge_type local_edge_type;
    typedef typename GraphType::local_vertex_type local_vertex_type;

    lvid_type lvid = graph.local_vid(vid);
    local_vertex_type lvertex = graph.l_vertex(lvid);
    vertex_type vertex = vertex_type(lvertex);
    if (in) {
      foreach(local_edge_type ledge, lvertex.in_edges()) {
        edge_type edge(ledge);
        vertex_type other(ledge.source());
        engine->elocks[ledge.id()].lock();
        edge_transform_fn(vertex, edge, other);
        engine->elocks[ledge.id()].unlock();
      }
    }

    if (out) {
      foreach(local_edge_type ledge, lvertex.out_edges()) {
        edge_type edge(ledge);
        vertex_type other(ledge.target());
        engine->elocks[ledge.id()].lock();
        edge_transform_fn(vertex, edge, other);
        engine->elocks[ledge.id()].unlock();
      }
    }
    return any();
  }
};





/**************************************************************************/
/*                                                                        */
/*        Defines the Synchronize Mirrors Task Type                       */
/*        This is basically an empty task.                                */
/*                                                                        */
/**************************************************************************/


typedef graphlab::empty synchronize_mirrors_task_param;

template <typename GraphType, typename EngineType>
struct synchronize_mirrors_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;
  typedef typename GraphType::vertex_type vertex_type;
  typedef typename GraphType::edge_type edge_type;

  virtual void combine(any& a, const any& b, const any& params) { }

  any exec(GraphType& graph, vertex_id_type vid, const any& params,
           EngineType* engine) {
    return any();
  }

  // fast path this empty operation
  bool fast_path(const any& params) {
    return true;
  }
};




/**************************************************************************/
/*                                                                        */
/*               Defines the Vertex Delta Task Type                       */
/*                                                                        */
/**************************************************************************/


struct vertex_delta_task_param {
  any delta_value;

  void save(oarchive& oarc) const {
    oarc << delta_value;
  }
  void load(iarchive& iarc) {
    iarc >> delta_value;
  }
};

template <typename GraphType, typename EngineType, typename T>
struct vertex_delta_task_descriptor: public gl3task_descriptor<GraphType, EngineType> {
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;
  typedef typename GraphType::vertex_type vertex_type;
  typedef typename GraphType::edge_type edge_type;

  typedef boost::function<void (vertex_type&, const T&)> delta_fn_type;
  delta_fn_type delta_fn;

  vertex_delta_task_descriptor(delta_fn_type delta_fn):delta_fn(delta_fn) { }

  virtual void combine(any& a, const any& b, const any& params) { }

  any exec(GraphType& graph, vertex_id_type vid, const any& params,
           EngineType* engine) {
    const vertex_delta_task_param& delta_param = params.as<vertex_delta_task_param>();

    lvid_type lvid = graph.local_vid(vid);
    vertex_type vertex(graph.l_vertex(lvid));

    engine->vlocks[lvid].lock();
    delta_fn(vertex, delta_param.delta_value.as<T>());
    engine->vlocks[lvid].unlock();
    return any();
  }

  // fast path unlock operations
  bool fast_path(const any& params) {
    return true;
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
