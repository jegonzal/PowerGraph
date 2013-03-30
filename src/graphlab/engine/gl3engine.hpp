#ifndef GRAPHLAB_GL3ENGINE_HPP
#define GRAPHLAB_GL3ENGINE_HPP
#include <vector>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/functional/hash.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <qthread/io.h>

#include <graphlab/engine/execution_status.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/parallel/qthread_tools.hpp>
#include <graphlab/util/tracepoint.hpp>
#include <graphlab/util/memory_info.hpp>
#include <graphlab/util/hashstream.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/engine/gl3task.hpp>
#include <graphlab/engine/gl3context.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>
#include <graphlab/rpc/async_consensus.hpp>
#include <graphlab/util/inplace_lf_queue2.hpp>
#include <graphlab/parallel/qthread_tools.hpp>
#include <graphlab/util/empty.hpp>



#include <graphlab/macros_def.hpp>
#define REMOVE_CONST_REF(REF) typename boost::remove_const<typename boost::remove_reference<REF>::type>::type
#define FRESULT(F) REMOVE_CONST_REF(typename boost::function<typename boost::remove_pointer<F>::type>::result_type)

#define NORMALIZE_FUNCTION(F) typename boost::function<typename boost::remove_pointer<F>::type>

namespace graphlab {

  /*
   * The task local contents of each virtual thread.
   * Whether this is a "vthread" task, and whether this is
   * currently running a vertex program. See the Design section in gl3engine
   * for details.
   */
  namespace GL3TLS {
    struct worker_task_local {
      bool is_vthread_task;
      bool is_vertex_program;
      lvid_type locked_lvid;
    };

    inline worker_task_local* GET_TASK_LOCAL() {
      return (worker_task_local*) qthread_get_tasklocal(sizeof(worker_task_local));
    }

    inline bool IN_VTHREAD_TASK() {
      return GET_TASK_LOCAL()->is_vthread_task;
    }

    inline bool IN_VERTEX_PROGRAM() {
      return GET_TASK_LOCAL()->is_vertex_program;
    }
    inline lvid_type LOCKED_LVID() {
      return GET_TASK_LOCAL()->locked_lvid;
    }

    inline void SET_IN_VTHREAD_TASK(bool val) {
      GET_TASK_LOCAL()->is_vthread_task = val;
    }

    inline void SET_IN_VERTEX_PROGRAM(bool val) {
      GET_TASK_LOCAL()->is_vertex_program = val;

    }
    inline void SET_LOCKED_LVID(lvid_type val) {
      GET_TASK_LOCAL()->locked_lvid = val;
    }
  }



/**
 * The GraphLab 3 engine describes a new purely asynchronous
 * graph abstraction built around the concept of non-blocking tasks and
 * function calls. This is made possible through the use of user-mode
 * threads provided by the qthread library.
 * It has the unfortunate side effect of making debugging quite difficult
 * since GDB does not know about user mode threads.
 *
 * ## Using the GL3 Engine
 * The GL3 engine works in a very different way compared to all the other
 * engines. Tasks, vertex programs, and pretty much nearly everything,
 * are non-blocking, thus permitting an arbitrary amount of interleaving between
 * different types of computation. In a sense, it is subtantially lower level
 * the the GraphLab 2 engines, but its flexibility is also substantially greater.
 * There are <b>many</b> ways to do anything.
 *
 * ### Key Typing Differences
 * Firstly, due to its differences, the gl3engine does not inherit from
 * iengine, unlike all the other engines.
 *
 * Secondly, to track changes in the graph data, it is necessary to define a
 * hash function for the vertex data. This is accomplished simply by defining a
 * hash_value function in the same namespace as the vertex data type
 * \code
 * size_t hash_value(vertex_data const& v) {
 *    return ... hash of vertex ...
 * }
 * \endcode
 * Note that it is not necessary for this hash function to be a good hash
 * function. All that is necessary is for this hash value to be different
 * whenever the data changes. This can be accomplished as simply as keeping a
 * change counter on the vertex data, and making sure to increment the counter
 * whenever data changes.
 *
 * ### Usage
 * Usage is easiest explained with an example.
 * \code
 *  typedef distributed_graph<float, graphlab::empty> graph_type;
 *  typedef gl3engine<graph_type> engine_type;
 *
 *  // define a neighborhood map reduce operation
 *  #define PAGERANK_MAP_REDUCE 0
 *  float pagerank_map(const graph_type::vertex_type& v) {
 *    return v.data() / v.num_out_edges();
 *  }
 *  void pagerank_combine(float& v1, const float& v2) {
 *    v1 += v2;
 *  }
 *
 *  void pagerank_program(engine_type::context_type& context,
 *                       graph_type::vertex_type& vertex,
 *                       const engine_type::message_type& unused) {
 *    float prev = vertex.data();
 *    // map reduce over neighbors
 *    vertex.data() = 0.15 + 0.85 *
 *        context.map_reduce<float>(PAGERANK_MAP_REDUCE, IN_EDGES);
 *
 *    float last_change = std::fabs((vertex.data()- prev) / vertex.num_out_edges());
 *    if (last_change > TOLERANCE) {
 *      // signals out neighbors if I change substantially
 *      context.broadcast_signal(OUT_EDGES);
 *    }
 *  }
 *  ... in main() ...
 *  // register the map reduce operation before usage
 *  // Each task registration must have a distinct ID ranging fro 0 to 223
 *  engine.register_map_reduce(PAGERANK_MAP_REDUCE,
 *                            pagerank_map,
 *                            pagerank_combine);
 *
 *  engine.set_vertex_program(pagerank_program);
 *  engine.signal_all();
 *  engine.wait();
 * \endcode
 *
 * Now, it is important to note that computation happens immediately on
 * the signal_all() call. It is non-blocking and asynchronous. The wait()
 * call simply waits for all activity to stop. This also means that more
 * interesting things can be done. For instance, a "sync" operation
 * be implemented trivially.
 * \code
 *
 * float pagerank_sum(graph_type::vertex_type v) { return v.data(); }
 * // ... instead of engine.wait... we substitute with the following code ...
 * while(!engine.try_wait()) {
 *   timer::sleep(1);
 *   float totalpr = graph.map_reduce_vertices<float>(pagerank_sum);
 *    dc.cout() << "Total PageRank = " << totalpr << "\n";
 * }
 * engine.wait();
 * \endcode
 *
 * As an example of There Is More Than One Way to Do It, we can also use
 * parfor_all_local_vertices() to implement non-dynamic asynchronous pagerank.
 *
 * \code
 * void pagerank_function(engine_type::context_type& context,
 *                       graph_type::vertex_type& vertex) {
 *  vertex.data() = 0.15 + 0.85 *
 *      context.map_reduce<float>(PAGERANK_MAP_REDUCE, IN_EDGES);
 * }
 * // and in main ...
 * for (size_t i = 0;i < iterations; ++i) {
 *    engine.parfor_all_local_vertices(pagerank_function);
 *    engine.wait();
 * }
 * \endcode
 *
 * Another really cool example. SGD Matrix Factorization on a bipartite graph.
 * The issue is when and how often to synchronize the parameters on the
 * vertices. Since everything interleave... We could simply run a synchronization
 * function while performing the SGD at the same time.
 * \code
 * void sgd_function(engine_type::context_type& context,
 *                    graph_type::edge_type& edge) {
 *   double pred = edge.source().data().pvec.dot(edge.target().data().pvec);
 *
 *   const float err = edge.data().obs - pred;
 *   vec_type delta;
 *   delta = [gradient of target() with respect to edge value]
 *   // send the changes as a delta operation
 *   context.send_delta(VERTEX_DELTA_TASK_ID, edge.target(), delta);
 *
 *   delta = [gradient of source() with respect to edge value]
 *   // send the changes as a delta operation
 *   context.send_delta(VERTEX_DELTA_TASK_ID, edge.source(), delta);
 * }
 * void sync_function(engine_type::context_type& context,
 *                   graph_type::vertex_type& vertex) {
 *   context.synchronize(vertex);
 * }
 *
 * ... and in main ...
 * for (size_t i = 0;i < ITERATIONS; ++i) {
 *   // interleave synchronization and SGD computation
 *   // both parfor functions take additional arguments to control
 *   // the number of user threads and the time slice size of each thread.
 *   engine.parfor_all_local_edges(sgd_function);
 *   engine.parfor_all_local_vertices(sync_function);
 *   engine.wait();
 *   double rmse = graph.map_reduce_edges<double>(extract_l2_error);
 *   dc.cout() << "RMSE = " << sqrt(rmse / graph.num_edges()) << std::endl;
 *   GAMMA = GAMMA * 0.9;
 * }
 * \endcode
 * ### Race Conditions
 * The arbitrary interleaving of tasks provide enormous flexibility, though
 * it comes at a price. It is up you to make sure race conditions on data do
 * not emerge.
 * Within a vertex program, we guarantee that the vertex data may not change
 * while your code is running, regardless of any other operations that are
 * going on. However, vertex data could change after "blocking" operations
 * through the context (such as map_reduce). Basically, blocking operations
 * unlock the vertex data, allowing other changes to happen.
 *
 * \code
 *  void pagerank_program(engine_type::context_type& context,
 *                       graph_type::vertex_type& vertex,
 *                       const engine_type::message_type& unused) {
 *    // the value of vertex.data() here...
 *    float ret = context.map_reduce<float>(PAGERANK_MAP_REDUCE, IN_EDGES);
 *    // could be different from the value of vertex.data() here.
 * \endcode
 *
 * ## Design
 * ### Work Sources and Threads
 * The engine manages work from 3 different sources.
 *
 * First are scheduler tasks which run user defined <b>vertex programs</b>.
 * Only one vertex program type can be active at any one time. This is
 * defined by the set_vertex_program() function. These are managed by
 * what the engine call "vthreads" which manage vertex program and
 * scheduler activity.
 *
 * The second are <b>user task threads</b> These are issued via the
 * launch_other_task() function, and can be any other task. Tasks launched via
 * the launch_other_task() function behave like a thread, with the exception
 * that the task should call gl3engine::poll() or gl3context::poll() every
 * now and then.
 *
 * The third type of tasks are subtasks. These are tasks created by the
 * spawn_subtask* family of functions. These describe remote
 * operations. For instance map_reduce_neighbors, transform neighbors,
 * DHT gather, etc. Each subtask type is numbered with a ID ranging from
 * 0 to 223. (IDs above 223 are reserved). Subtask types are implemented
 * by inheriting from gl3task_descriptor.
 *
 */
template <typename GraphType, typename MessageType = empty>
class gl3engine {
 public:
  typedef GraphType graph_type;
  typedef typename GraphType::vertex_type vertex_type;
  typedef typename GraphType::edge_type edge_type;
  typedef typename GraphType::vertex_data_type vertex_data_type;
  typedef typename GraphType::edge_data_type edge_data_type;
  typedef typename GraphType::local_vertex_type    local_vertex_type;
  typedef typename GraphType::local_edge_type      local_edge_type;
  typedef typename GraphType::lvid_type            lvid_type;

  typedef gl3engine<GraphType> engine_type;
  typedef gl3context<engine_type> context_type;
  typedef gl3task_descriptor<GraphType, engine_type> task_descriptor_type;
  typedef MessageType message_type;

  typedef boost::function<void (context_type&,
                                vertex_type&,
                                const message_type&)> vertex_program_type;



  // \internal \brief The base type of all schedulers
  typedef ischeduler<message_type> ischeduler_type;

  // sub task
  struct task {
    task* next;
    vertex_id_type vid;
    any param;
    unsigned char task_id;
    procid_t origin;
    size_t handle;
  };
 private:
  dc_dist_object<gl3engine> rmi;

  // maximum number of vthreads to support
  size_t num_vthreads;
  // number of worker threads
  size_t ncpus;
  // A reference to the graph
  graph_type& graph;
  // the subtask types
  task_descriptor_type* task_types[256];

  // The vertex program currently being scheduled
  vertex_program_type active_function;
  // A schedule of the tasks to be executed
  std::vector<inplace_lf_queue2<task>* > local_tasks;

  struct vertex_local_data_type {
    size_t hash;
    message_type buffered_messages;
    bool is_running; bool has_message;
    simple_spinlock lock;
  };
  std::vector<vertex_local_data_type> vertex_local_data;

  atomic<size_t> programs_completed;
  atomic<size_t> tasks_completed;

  size_t total_programs_completed;
  size_t total_tasks_completed;

  atomic<size_t> active_vthread_count; // number of vthreads launched
  atomic<size_t> thread_counter; // just used to increment and load balance through queues.
  atomic<size_t> num_working_threads; // number of vthreads not in a yield loop

  std::vector<mutex> worker_mutex;
  //! The scheduler
  ischeduler_type* scheduler_ptr;

  async_consensus* consensus;


  // The DHT.
  static const size_t NUM_DHTS = 256;
  mutex dht_lock[NUM_DHTS];
  boost::unordered_map<size_t, any> dht[NUM_DHTS]; // a large dht

  /*
   * \brief The vertex locks protect access to vertex specific
   * data-structures including
   * \ref graphlab::synchronous_engine::gather_accum
   * and \ref graphlab::synchronous_engine::messages.
   */
  std::vector<simple_spinlock> vlocks;

  std::vector<mutex> subtask_thread_alive_lock;
  std::vector<char> subtask_thread_alive; // if 0, it is dead, if 1 it is quiting, if 2 it is alive

  /*
   * \brief The elocks protect individual edges during gather and
   * scatter.  Technically there is a potential race since gather
   * and scatter can modify edge values and can overlap.  The edge
   * lock ensures that only one gather or scatter occurs on an edge
   * at a time.
   */
  std::vector<simple_spinlock> elocks;




  graphlab::qthread_group execution_group;

  // we have to friend basically every subtask
  friend struct gl3task_descriptor<GraphType, engine_type>;
  template <typename, typename, typename>
  friend struct map_reduce_neighbors_task_descriptor;
  friend struct broadcast_task_descriptor<GraphType, engine_type>;
  friend struct edge_transform_task_descriptor<GraphType, engine_type>;
  friend struct dht_gather_task_descriptor<GraphType, engine_type>;
  friend struct dht_scatter_task_descriptor<GraphType, engine_type>;
  friend struct synchronize_mirrors_task_descriptor<GraphType, engine_type>;
  template <typename, typename, typename>
  friend struct vertex_delta_task_descriptor;

  friend struct gl3context<engine_type>;
  context_type base_context;


/**************************************************************************/
/*                                                                        */
/*                        Basic Interface Function                        */
/*                                                                        */
/**************************************************************************/

 public:
  gl3engine(distributed_control& dc, graph_type& graph,
            graphlab_options opts = graphlab_options()):
      rmi(dc, this), graph(graph){
    rmi.barrier();
    num_vthreads = 10000;
    ncpus = opts.get_ncpus();
    worker_mutex.resize(ncpus);


    // read the options
    std::vector<std::string> keys = opts.get_engine_args().get_option_keys();
    foreach(std::string opt, keys) {
      if (opt == "num_vthreads") {
        opts.get_engine_args().get_option("num_vthreads", num_vthreads);
        if (rmi.procid() == 0)
          logstream(LOG_EMPH) << "Engine Option: num_vthreads = "
                              << num_vthreads << std::endl;
      } else {
        logstream(LOG_FATAL) << "Unexpected Engine Option: " << opt << std::endl;
      }
    }

    if (opts.get_scheduler_type() == "") {
      opts.set_scheduler_type("fifo");
    }
    // create the scheduler
    scheduler_ptr = scheduler_factory<message_type>::
                    new_scheduler(graph.num_local_vertices(),
                                  opts);
    // construct the termination consensus object
    // 1 thread
    consensus = new async_consensus(rmi.dc(), 1);
    // construct the locks
    vlocks.resize(graph.num_local_vertices());
    elocks.resize(graph.num_local_edges());
    vertex_local_data.resize(graph.num_local_vertices());
    // construct the queues
    local_tasks.resize(ncpus);
    subtask_thread_alive.resize(ncpus, 0);
    subtask_thread_alive_lock.resize(ncpus);
    for (size_t i = 0;i < local_tasks.size(); ++i) {
      local_tasks[i] = new inplace_lf_queue2<task>();
    }

    // make default registrations
    task_types[GL3_BROADCAST_TASK_ID] = new broadcast_task_descriptor<GraphType, engine_type>();
    task_types[GL3_DHT_GATHER_TASK_ID] = new dht_gather_task_descriptor<GraphType, engine_type>();
    task_types[GL3_SYNCHRONIZE_MIRRORS_TASK_ID] = new synchronize_mirrors_task_descriptor<GraphType, engine_type>();
    //
    // pre init
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)graph.num_local_vertices(); ++i) {
      vertex_local_data[i].is_running = false;
      vertex_local_data[i].has_message = false;
      vertex_local_data[i].hash = get_vertex_data_hash_lvid((lvid_type)i);
    }
    scheduler_ptr->start();

    //reset counters
    programs_completed = 0;
    tasks_completed = 0;
    total_programs_completed = 0;
    total_tasks_completed = 0;
    active_vthread_count = 0;
    num_working_threads = 0;


    graphlab::qthread_tools::init(ncpus, 128 * 1024);

    base_context.engine = this;
    base_context.lvid = (lvid_type)(-1);
    rmi.full_barrier();
    qt_begin_blocking_action();
  }

  ~gl3engine() {
    qt_end_blocking_action();
  }

  procid_t procid() const {
    return rmi.procid();
  }


  context_type& get_context() {
    return base_context;
  }


/**************************************************************************/
/*                                                                        */
/*                          Scheduler Functions                           */
/*                                                                        */
/**************************************************************************/
 public:
  void set_vertex_program(vertex_program_type uf) {
    active_function = uf;
    rmi.barrier();
  }

 private:
  void internal_signal(const vertex_type& vtx,
                       const message_type& message = message_type()) {
    internal_schedule(vtx.local_id(), message);
  } // end of schedule


  void internal_schedule(const lvid_type lvid,
                         const message_type& message) {
    scheduler_ptr->schedule(lvid, message);
    size_t target_queue = thread_counter++;

    if (active_vthread_count < num_vthreads
        && num_working_threads < ncpus) {
      active_vthread_count.inc();
      num_working_threads.inc();
      execution_group.launch(boost::bind(&gl3engine::vthread_start,
                                         this,
                                         target_queue));
      consensus->cancel();
    }
  }

 public:

  void signal(vertex_id_type gvid,
              const message_type& message = message_type()) {
    internal_signal_gvid(gvid, message);
  }

  void internal_signal_gvid(vertex_id_type gvid,
                            const message_type& message = message_type()) {
    if (graph.is_master(gvid)) {
      internal_signal(graph.vertex(gvid), message);
    }
  }



  void rpc_signal(vertex_id_type gvid,
                  const message_type& message) {
    internal_signal(graph.vertex(gvid), message);
  }


  void signal_all(const message_type& message = message_type(),
                  const std::string& order = "shuffle") {
    logstream(LOG_DEBUG) << rmi.procid() << ": Schedule All" << std::endl;
    // allocate a vector with all the local owned vertices
    // and schedule all of them.
    std::vector<vertex_id_type> vtxs;
    vtxs.reserve(graph.num_local_own_vertices());
    for(lvid_type lvid = 0;
        lvid < graph.get_local_graph().num_vertices();
        ++lvid) {
      if (graph.l_vertex(lvid).owner() == rmi.procid()) {
        vtxs.push_back(lvid);
      }
    }

    if(order == "shuffle") {
      graphlab::random::shuffle(vtxs.begin(), vtxs.end());
    }
    foreach(lvid_type lvid, vtxs) {
      internal_schedule(lvid, message);
    }
    rmi.barrier();
  } // end of schedule all


  void signal_vset(const vertex_set& vset,
                   const message_type& message = message_type(),
                   const std::string& order = "shuffle") {
    logstream(LOG_DEBUG) << rmi.procid() << ": Schedule All" << std::endl;
    // allocate a vector with all the local owned vertices
    // and schedule all of them.
    std::vector<vertex_id_type> vtxs;
    vtxs.reserve(graph.num_local_own_vertices());
    for(lvid_type lvid = 0;
        lvid < graph.get_local_graph().num_vertices();
        ++lvid) {
      if (graph.l_vertex(lvid).owner() == rmi.procid() &&
          vset.l_contains(lvid)) {
        vtxs.push_back(lvid);
      }
    }

    if(order == "shuffle") {
      graphlab::random::shuffle(vtxs.begin(), vtxs.end());
    }
    foreach(lvid_type lvid, vtxs) {
      internal_schedule(lvid, message);
    }
    rmi.barrier();
  }


/**************************************************************************/
/*                                                                        */
/*                      Task Registration Functions                       */
/*                                                                        */
/**************************************************************************/

 public:

  /**
   * Registers a map_reduce over neighborhood task. This task can then
   * be called from within an vertex program. After registration, the
   * map_reduce can be called from within a  vertex program using the
   * gl3context::map_reduce() function.
   *
   * For instance, a PageRank operation can be defined as follows.
   * \code
   * typedef graphlab::distributed_graph<float, graphlab::empty> graph_type;
   *
   * const unsigned char PAGERANK_MAP_REDUCE = 0;
   *
   * float pagerank_map(const graph_type::vertex_type& v) {
   *   return 0.85 * v.data() / v.num_out_edges();
   * }
   * void pagerank_combine(float& v1, const float& v2) {
   *   v1 += v2;
   * }
   * ... in main() ...
   * engine.register_map_reduce(PAGERANK_MAP_REDUCE,
   *                            pagerank_map,
   *                            pagerank_combine);
   * ... in the vertex program ...
   * vertex.data() = 0.15 + 0.85 *
   *   context.map_reduce<float>(PAGERANK_MAP_REDUCE, IN_EDGES);
   *
   * \endcode
   *
   * You can also use a more general map function type that takes in the
   * edge and both vertex endpoints. Note that the edge is not const and
   * modifiable.
   * \code
   * float pagerank_map(const graph_type::vertex_type& center,
   *                    graph_type::edge_type& edge,
   *                    const graph_type::vertex_type& other) {
   *   return 0.85 * other.data() / other.num_out_edges();
   * }
   * void pagerank_combine(float& v1, const float& v2) {
   *   v1 += v2;
   * }
   * \endcode
   *
   * Every registered task must have a distinct ID ranging from 0 to 223
   */
  template <typename MapFn, typename CombineFn>
  void register_map_reduce(size_t id,
                           MapFn mapfn,
                           CombineFn combinefn) {
    ASSERT_LT(id, 224);
    task_types[id] = map_reduce_neighbors_task_descriptor<GraphType, engine_type, FRESULT(MapFn)>::
        create_map_reduce_task_impl::template create<NORMALIZE_FUNCTION(MapFn)>(mapfn, combinefn);
    rmi.barrier();
  }

  /**
   * Registers an edge transform that can
   * be called from within an vertex program. After registration, the task
   * can be called from within a vertex program using the
   * gl3context::edge_transform() function in a blocking, or a
   * non-blocking fashion.
   *
   * For instance, an operation that sets all the edge values to 0
   * can be defined as follows.
   *
   * \code
   * typedef graphlab::distributed_graph<float, int> graph_type;
   *
   * const unsigned char SET_TO_ZERO = 0;
   *
   * float set_edge_to_zero(const graph_type::vertex_type& center,
   *                        graph_type::edge_type& edge,
   *                        const graph_type::vertex_type& other) {
   *   edge.data() = 0;
   * }
   *
   * ... in main() ...
   * engine.register_edge_transform(SET_TO_ZERO,
   *                                set_edge_to_zero);
   * ... in the vertex program ...
   *   // This will block until the task is complete
   *   context.edge_transform(SET_TO_ZERO, IN_EDGES);
   *   // Alternative, This will not block
   *   context.edge_transform(SET_TO_ZERO, IN_EDGES, false);
   * \endcode
   *
   * Every registered task must have a distinct ID ranging from 0 to 223
   */
  void register_edge_transform(size_t id,
                               boost::function<void (const vertex_type&,
                                                     edge_type&,
                                                     const vertex_type&)> etfn) {
    ASSERT_LT(id, 224);
    task_types[id] = new edge_transform_task_descriptor<GraphType, engine_type>(etfn);
    rmi.barrier();
  }

  /**
   * Registers a Scatter Delta operation on the builtin DHT.
   * This allows you to send an object to a DHT key, and have that key
   * receive it and run a user defined function to combine it. The scatter
   * delta operation can be used anywhere with the gl3context::dht_scatter()
   * operation.
   * Example:
   *
   * \code
   * // DHT should contain integers, and I send increments to the integer
   * void scatter_increment(any& val, const any& other) {
   *   if (val.empty()) val = other;
   *   else val.as<int>() += other.as<int>();
   * }
   *
   * .. in main ..
   * // build up a change set
   * boost::unordered_map<size_t, any> changes;
   * // increment dht[0] by 1, and increment dht[10] by 5
   * changes[0] = int(1);
   * changes[10] = int(5);
   * // issue the changes
   * engine.get_context().dht_scatter(changes);
   * engine.wait();
   *
   * // now read the values of dht[0] and dht[10]
   * std::vector<size_t> entries; entries.push_back(0); entries.push_back(10);
   * boost::unordered_map<size_t, any> dhtvalues = engine.get_context.dht_gather(entries);
   * // dht_values[0] and dht_values[10] will now contain the dht values after
   * // the change
   * \endcode
   */
  void register_dht_scatter(size_t id,
                           boost::function<void (any&, const any&)> scatter_fn) {
    dht_scatter_task_descriptor<GraphType, engine_type>* desc =
        new dht_scatter_task_descriptor<GraphType, engine_type>;
    desc->scatter_fn = scatter_fn;
    task_types[id] = desc;
    rmi.barrier();
  }

  /**
   * Registers a task which performs a delta operation on a vertex.
   * A delta operation allows you to send an object to a (master) vertex (from the
   * gl3context::send_delta() function ) and have that (master) vertex receive it,
   * running a user defined function to combine it. The delta operation
   * can be used anytime you have a reference to a vertex_type and a context
   * object. Thus this can be used from within parfor_all_local_edges()
   * or parfor_all_local_vertices() as well as vertex programs.
   *
   * This example uses deltas to count the number of edges on each
   * vertex, writing it as the vertex value.
   *
   * Deltas are not required to be additive. Deltas are not combined,
   * and are each sent individually as a seperate message.
   *
   * Also, if the delta is sent to a mirror vertex, it is combined immediately
   * with the local vertex data, and it also sent to the master vertex.
   * Synchronization of the master vertex with its mirrors does not occur
   * automatically, but must be managed externally via additional synchronize
   * operations. (Either in a non-blocking fashion using
   * parfor_all_local_vertices() and gl3context::synchronize(), or synchronously
   * using distributed_graph::synchronize() )
   * \code
   * typedef distributed_graph<int, empty> graph_type;
   *
   * #define VERTEX_DELTA_TASK_ID 0
   *
   * void delta_function(graph_type::vertex_type& vtx, int delta) {
   *   vtx.data() += delta;
   * }
   *
   * void count_function(engine_type::context_type& context,
   *                     graph_type::edge_type& edge) {
   *   context.send_delta(VERTEX_DELTA_TASK_ID, edge.source(), (int)1);
   *   context.send_delta(VERTEX_DELTA_TASK_ID, edge.target(), (int)1);
   * }
   *
   * ... in main ...
   * engine.register_vertex_delta<int>(VERTEX_DELTA_TASK_ID, delta_function);
   * engine.parfor_all_local_edges(count_function);
   * engine.sync();
   * \endcode
   *
   *
   * Every registered task must have a distinct ID ranging from 0 to 223
   */
  template <typename T>
  void register_vertex_delta(size_t id,
                             boost::function<void (vertex_type&, const T&)> deltafn) {
   ASSERT_LT(id, 224);
   vertex_delta_task_descriptor<GraphType, engine_type, T>* desc =
        new vertex_delta_task_descriptor<GraphType, engine_type, T>(deltafn);
    task_types[id] = desc;
    rmi.barrier();
  }


/**************************************************************************/
/*                                                                        */
/*               NonBlocking Parallel for over local edges                */
/*                                                                        */
/**************************************************************************/

 private:
  void parfor_local_edges_loop(context_type& context,
                               size_t start,
                               size_t skip,
                               size_t timeslice_ms,
                               boost::function<void (context_type&, edge_type&)> edgefn) {
    timer ti;
    ti.start();
    for (lvid_type i = start; i < graph.num_local_vertices(); i += skip) {
      local_vertex_type lvertex = graph.l_vertex(i);
      foreach(local_edge_type ledge, lvertex.in_edges()) {
        edge_type edge(ledge);
        elocks[ledge.id()].lock();
        edgefn(context, edge);
        elocks[ledge.id()].unlock();
      }
      if (ti.current_time_millis() > timeslice_ms) {
        poll();
        ti.start();
      }
    }
  }
 public:
  /**
   * Loops a function over all local edges. The function could transform
   * the data on the edge. Non-blocking.
   *
   * For instance, we could write a matrix factorization SGD:
   * \code
   * void sgd_function(engine_type::context_type& context,
   *                    graph_type::edge_type& edge) {
   *   double pred = edge.source().data().pvec.dot(edge.target().data().pvec);
   *
   *   const float err = edge.data().obs - pred;
   *   vec_type delta;
   *   delta = [gradient of target() with respect to edge value]
   *   // send the changes as a delta operation
   *   context.send_delta(VERTEX_DELTA_TASK_ID, edge.target(), delta);
   *
   *   delta = [gradient of source() with respect to edge value]
   *   // send the changes as a delta operation
   *   context.send_delta(VERTEX_DELTA_TASK_ID, edge.source(), delta);
   * }
   * void sync_function(engine_type::context_type& context,
   *                   graph_type::vertex_type& vertex) {
   *   context.synchronize(vertex);
   * }
   * \endcode
   *
   * The additional arguments control the number of user threads to spawn
   * and the time slice given to each thread. The default number of user threads
   * is equal tot num_vthreads and the default timeslice is 10ms
   */
  void parfor_all_local_edges(boost::function<void (context_type&, edge_type&)> edgefn,
                              size_t num_threads = (size_t)-1,
                              size_t timeslice_ms = 10) {
    if (num_threads == (size_t)(-1)) num_threads = num_vthreads;
    for (size_t i = 0;i < num_threads ; ++i) {
      launch_other_task(boost::bind(&gl3engine::parfor_local_edges_loop,
                                    this,
                                    _1,
                                    i,
                                    num_threads,
                                    timeslice_ms,
                                    edgefn));
    }
  }

/**************************************************************************/
/*                                                                        */
/*               NonBlocking Parallel for over vertices                   */
/*                                                                        */
/**************************************************************************/

// This is almost like a vertex program. But without the
// scheduling capabilities. We will emulate it as such

 private:
  void parfor_local_vertices_loop(context_type& context,
                                  size_t start,
                                  size_t skip,
                                  size_t timeslice_ms,
                                  boost::function<void (context_type&, vertex_type&)> vertex_fn) {
    GL3TLS::SET_IN_VERTEX_PROGRAM(true);
    timer ti;
    ti.start();
    for (lvid_type i = start; i < graph.num_local_vertices(); i += skip) {
      local_vertex_type lvertex = graph.l_vertex(i);
      if (lvertex.owned()) {
        vertex_type vertex(lvertex);
        // set up the context
        context.lvid = i;
        // lock and call the function
        while (!vlocks[i].try_lock()) qthread_yield();
        GL3TLS::SET_LOCKED_LVID(i);
        vertex_fn(context, vertex);

        size_t newhash = get_vertex_data_hash_lvid(i);
        // if the hash changed, broadcast
        if (newhash != vertex_local_data[i].hash) {
          vertex_local_data[i].hash = newhash;
          local_vertex_type lvertex(graph.l_vertex(i));
          rmi.remote_call(lvertex.mirrors().begin(), lvertex.mirrors().end(),
                          &engine_type::sync_vdata,
                          lvertex.global_id(),
                          lvertex.data());
        }
        vlocks[i].unlock();
        GL3TLS::SET_LOCKED_LVID((lvid_type)(-1));
      }
      if (ti.current_time_millis() > timeslice_ms) {
        poll();
        ti.start();
      }
    }
    GL3TLS::SET_IN_VERTEX_PROGRAM(false);
  }
 public:
  /**
   * Runs a function over all vertices. Non-blocking.
   * Behaves very much like like a vertex program, but without the scheduling
   * capabilities.
   *
   * This provides <b> yet another </b> way to implement pagerank.
   * \code
   * void update_function(engine_type::context_type& context,
   *                      graph_type::vertex_type& vertex) {
   *   vertex.data() = 0.15 + 0.85 *
   *            context.map_reduce<float>(PAGERANK_MAP_REDUCE, IN_EDGES);
   * }
   * ... in main ...
   * for (size_t i = 0;i < iterations; ++i) {
   *    engine.parfor_all_local_vertices(update_function);
   *    engine.wait();
   * }
   * \endcode
   *
   * The additional arguments control the number of user threads to spawn
   * and the time slice given to each thread. The default number of user threads
   * is equal tot num_vthreads and the default timeslice is 10ms
   */
  void parfor_all_local_vertices(boost::function<void (context_type&, vertex_type&)> vertex_fn,
                                 size_t num_threads = (size_t)-1,
                                 size_t timeslice_ms = 10) {
    if (num_threads == (size_t)(-1)) num_threads = num_vthreads;
    for (size_t i = 0;i < num_threads ; ++i) {
      launch_other_task(boost::bind(&gl3engine::parfor_local_vertices_loop,
                                    this,
                                    _1,
                                    i,
                                    num_threads,
                                    timeslice_ms,
                                    vertex_fn));
    }
  }

/**************************************************************************/
/*                                                                        */
/*                      Custom User Task Management                       */
/*                                                                        */
/**************************************************************************/

 private:
  void task_start(boost::function<void (context_type&)> fn, size_t id) {
    GL3TLS::SET_IN_VTHREAD_TASK(false);
    GL3TLS::SET_IN_VERTEX_PROGRAM(false);
    context_type context;
    context.engine = this;
    context.lvid = lvid_type(-1);
    context.thread_id = id;
    fn(context);
  }
 public:
  void launch_other_task(boost::function<void (context_type&)> fn) {
    size_t target_queue = thread_counter++;
    execution_group.launch(boost::bind(&gl3engine::task_start,
                                       this,
                                       fn,
                                       target_queue));
  }


/**************************************************************************/
/*                                                                        */
/*                           Subtask Management                           */
/*                                                                        */
/**************************************************************************/

 private:
  struct __attribute((__may_alias__)) future_combiner {
    any param;
    any result;
    atomic<procid_t> count_down;
    unsigned char task_id;
    simple_spinlock lock;
  };

  /**
   * Spawns a subtask unassociated with the graph on a collection of target
   * machines, sending each a task parameter object. target_machines must be
   * the same length as task_param.
   */
  any spawn_subtask(unsigned char task_id,
                 const std::vector<procid_t>& target_machines,
                 const std::vector<any>& task_param,
                 bool no_reply = false) {
    if (target_machines.size() == 1 && target_machines[0] == rmi.procid()) {
      // fast path for purely local task
      return task_types[task_id]->exec(graph, vertex_id_type(-1), task_param[0], this);
    }
    future_combiner combiner;
    combiner.count_down = target_machines.size();
    combiner.task_id = task_id;
    combiner.param = task_param;

    conditional_serialize<vertex_data_type> cs;
    size_t cb = reinterpret_cast<size_t>(&combiner);
    if (no_reply) cb = 0;

    bool has_self_task = false;
    size_t self_task_id = 0;
    for (size_t i = 0;i < target_machines.size(); ++i) {
      if (target_machines[i] == rmi.procid()) {
        has_self_task = true;
        self_task_id = i;
      } else {
        ASSERT_LT(target_machines[i], rmi.numprocs());
        rmi.remote_call(target_machines[i],
                        &engine_type::rpc_receive_task,
                        task_id,
                        vertex_id_type(-1),
                        cs,
                        task_param[i],
                        rmi.procid(),
                        cb);
      }
    }
    if (has_self_task) {
      // we execute my own subtasks inplace
      any ret = task_types[task_id]->exec(graph, vertex_id_type(-1), task_param[self_task_id], this);
      if (!no_reply) task_reply(&combiner, ret);
    }
    if (!no_reply) {
      if (GL3TLS::IN_VTHREAD_TASK()) num_working_threads.dec();
      while(combiner.count_down != 0) {
        qthread_yield();
      }
      if (GL3TLS::IN_VTHREAD_TASK()) num_working_threads.inc();
    }
    return combiner.result;
  }

  void unlock_vlock_for_subtask_execution(lvid_type lvid) {
    bool in_vertex_program = GL3TLS::IN_VERTEX_PROGRAM();
    lvid_type locked_lvid = GL3TLS::LOCKED_LVID();
    if (in_vertex_program && lvid == locked_lvid) vlocks[lvid].unlock();
  }

  void relock_vlock_for_subtask_execution(lvid_type lvid) {
    bool in_vertex_program = GL3TLS::IN_VERTEX_PROGRAM();
    lvid_type locked_lvid = GL3TLS::LOCKED_LVID();
    if (in_vertex_program && lvid == locked_lvid) {
      while (!vlocks[lvid].try_lock()) qthread_yield();
    }
  }
  /**
   * Spawns a single subtask on the corresponding master vertex of a given lvid.
   * If lvid already the master, the subtask is executed immediately.
   */
  any spawn_subtask_to_master(lvid_type lvid,
                           unsigned char task_id,
                           const any& task_param,
                           bool no_reply = false) {
    local_vertex_type lvertex(graph.l_vertex(lvid));
    if (lvertex.owned()) {
      // fast path if I am the master
      unlock_vlock_for_subtask_execution(lvid);
      any ret = task_types[task_id]->exec(graph, lvertex.global_id(), task_param, this);
      relock_vlock_for_subtask_execution(lvid);
      return ret;
    }
    future_combiner combiner;
    combiner.count_down = 1;
    combiner.task_id = task_id;
    combiner.param = task_param;
    /*
    logstream(LOG_EMPH) << "Creating Subtask type "<< (int)task_id
                        << " on vertex " << graph.l_vertex(lvid).global_id() << " Handle " << combiner
                        << "\n";
    */
    // do not send vertex data
    conditional_serialize<vertex_data_type> cs;

    size_t cb = reinterpret_cast<size_t>(&combiner);
    if (no_reply) cb = 0;
    rmi.remote_call(lvertex.owner(),
                    &engine_type::rpc_receive_task,
                    task_id,
                    lvertex.global_id(),
                    cs,
                    task_param,
                    rmi.procid(),
                    cb);

    // we execute my own subtasks inplace
    /*
    logstream(LOG_EMPH) << "Execing subtask type " << (int)(task_id)
                        << " on vertex " << lvertex.global_id() << "\n";
    */
    // if we are in an update function
    // unlock the task locks so we don't dead-lock with the subtask
    if (!no_reply) {
      unlock_vlock_for_subtask_execution(lvid);
      if (GL3TLS::IN_VTHREAD_TASK()) num_working_threads.dec();
      while(combiner.count_down != 0) {
        qthread_yield();
      }
      if (GL3TLS::IN_VTHREAD_TASK()) num_working_threads.inc();
      relock_vlock_for_subtask_execution(lvid);
    }
    return combiner.result;
  }


  // sends a task to all mirrors
  // vlocks[lvid] must be locked
  /**
   * Spawns a single subtask on all mirror vertices of a given master vertex.
   * lvid must be a master.
   */
  any spawn_subtask(lvid_type lvid,
                 unsigned char task_id,
                 const any& task_param,
                 bool no_reply = false) {
    ASSERT_TRUE(graph.l_is_master(lvid));
    local_vertex_type lvertex(graph.l_vertex(lvid));
    if (graph.l_vertex(lvid).num_mirrors() == 0) {
      // fast path for purely local task
      unlock_vlock_for_subtask_execution(lvid);
      any ret = task_types[task_id]->exec(graph, lvertex.global_id(), task_param, this);
      relock_vlock_for_subtask_execution(lvid);
      return ret;
    }
    future_combiner combiner;
    combiner.count_down = graph.l_vertex(lvid).num_mirrors() + 1;
    combiner.task_id = task_id;
    combiner.param = task_param;
    /*
    logstream(LOG_EMPH) << "Creating Subtask type "<< (int)task_id
                        << " on vertex " << graph.l_vertex(lvid).global_id() << " Handle " << combiner
                        << "\n";
    */
    size_t newhash = get_vertex_data_hash_lvid(lvid);
    conditional_serialize<vertex_data_type> cs;
    if (newhash != vertex_local_data[lvid].hash) {
      cs.hasval = true;
      cs.val = lvertex.data();
      vertex_local_data[lvid].hash = newhash;
    }
    size_t cb = reinterpret_cast<size_t>(&combiner);
    if (no_reply) cb = 0;
    rmi.remote_call(lvertex.mirrors().begin(), lvertex.mirrors().end(),
                    &engine_type::rpc_receive_task,
                    task_id,
                    lvertex.global_id(),
                    cs,
                    task_param,
                    rmi.procid(),
                    cb);

    // we execute my own subtasks inplace
    /*
    logstream(LOG_EMPH) << "Execing subtask type " << (int)(task_id)
                        << " on vertex " << lvertex.global_id() << "\n";
    */
    // if we are in an update function,
    // unlock the task locks so we don't dead-lock with the subtask
    unlock_vlock_for_subtask_execution(lvid);
    any ret = task_types[task_id]->exec(graph, lvertex.global_id(), task_param, this);
    if (!no_reply) {
      if (GL3TLS::IN_VTHREAD_TASK()) num_working_threads.dec();
      task_reply(&combiner, ret);
      while(combiner.count_down != 0) {
        qthread_yield();
      }
      if (GL3TLS::IN_VTHREAD_TASK()) num_working_threads.inc();
    }
    relock_vlock_for_subtask_execution(lvid);
    return combiner.result;
  }




  /// Task replies come through here
  void task_reply_rpc(size_t handle, any& val) {
    future_combiner* combiner = reinterpret_cast<future_combiner*>(handle);
    task_reply(combiner, val);
  }


  /// Task replies come through here
  void task_reply(future_combiner* combiner, any& val) {
    //logstream(LOG_EMPH) << "some subtask completion on handle " << combiner << "\n";
    combiner->lock.lock();
    ASSERT_GT(combiner->count_down, 0);
    if (!combiner->result.empty()) {
      task_types[combiner->task_id]->combine((combiner->result),
                                             val, combiner->param);
    } else {
      combiner->result = val;
    }
    combiner->lock.unlock();
    combiner->count_down.dec();
  }

  /// Receives a subtask from a remote machine.
  void rpc_receive_task(unsigned char task_id,
                        vertex_id_type vid,
                        conditional_serialize<vertex_data_type>& vdata,
                        const any& param,
                        procid_t caller,
                        size_t handle) {
    //logstream(LOG_EMPH) << "Receiving subtask on handle " << (void*)(handle) << "\n";
    if (vid != vertex_id_type(-1)) {
      lvid_type lvid = graph.local_vid(vid);
      if (vdata.hasval) {
        vlocks[lvid].lock();
        graph.l_vertex(lvid).data() = vdata.val;
        vlocks[lvid].unlock();
      }
    }
    // do we fast path it?
    if (task_types[task_id]->fast_path(param)) {
      // fast path
      any ret = task_types[task_id]->exec(graph, vid, param, this);
      // return to origin
      if (handle != 0) {
        rmi.remote_call(caller,
                        &gl3engine::task_reply_rpc,
                        handle,
                        ret);
      }
      tasks_completed.inc();
    } else {
      // regular path
      task* t = new task;
      t->origin = caller;
      t->param = param;
      t->task_id = task_id;
      t->handle = handle;
      t->vid = vid;
      size_t target_queue = (thread_counter++) % ncpus;
      local_tasks[target_queue]->enqueue(t);

      // if the subtask thread handing this task is not alive any more,
      // re-create it.
      //std::cout << "subtask thread " << target_queue << " stat: " << (int)subtask_thread_alive[target_queue] << "\n";
      if (subtask_thread_alive[target_queue] != 2) {
        subtask_thread_alive_lock[target_queue].lock();
        if (subtask_thread_alive[target_queue] == 0) {
          subtask_thread_alive[target_queue] = 2;
          subtask_thread_alive_lock[target_queue].unlock();
          execution_group.launch(boost::bind(&gl3engine::subtask_thread_start,
                                             this,
                                             target_queue));
        } else {
          subtask_thread_alive_lock[target_queue].unlock();
        }
        consensus->cancel();
      }
    }
  }

  /**
   * Reads a subtask from a queue and runs it.
   */
  bool exec_subtasks(size_t worker) {
    GL3TLS::SET_IN_VERTEX_PROGRAM(false);
    bool ret = false;
    if (worker_mutex[worker].try_lock()) {
      //rmi.dc().handle_incoming_calls(worker, ncpus);
      //
      task* tasks = local_tasks[worker]->dequeue_all();
      if (tasks != NULL) {
        ret = true;
        // execute tasks
        while (!local_tasks[worker]->end_of_dequeue_list(tasks)) {
          task* cur = tasks;
          // execute cur
          /*
             logstream(LOG_EMPH) << "Execing subtask type " << (int)(cur->task_id)
                                 << " on vertex " << cur->vid << "\n";
                                 */
          any ret = task_types[cur->task_id]->exec(graph, cur->vid, cur->param, this);
          // return to origin
          if (cur->handle != 0) {
            rmi.remote_call(cur->origin,
                            &gl3engine::task_reply_rpc,
                            cur->handle,
                            ret);
          }
          tasks_completed.inc();
          // get the next task in the queue
          while(inplace_lf_queue2<task>::get_next(tasks) == NULL) asm volatile ("" : : : "memory");
          tasks = inplace_lf_queue2<task>::get_next(tasks);
          delete cur;
        }
      }
      worker_mutex[worker].unlock();
    }
    return ret;
  }

/**************************************************************************/
/*                                                                        */
/*                       Scheduler Task Management                        */
/*                                                                        */
/**************************************************************************/
 private:
  /// Computes the hash of a given vertex
  size_t get_vertex_data_hash_lvid(lvid_type lvid) {
    boost::hash<vertex_data_type> hasher;
    return hasher(graph.l_vertex(lvid).data());
  }

  /// Writes the vertex data into.
  void sync_vdata(vertex_id_type vid,
                  vertex_data_type& vdata) {
    local_vertex_type lvertex(graph.vertex(vid));
    ASSERT_FALSE(lvertex.owned());
    vlocks[lvertex.id()].lock();
    lvertex.data() = vdata;
    vlocks[lvertex.id()].unlock();
  }


  /**
   * Executed before the scheduled function is actually run.
   * Verifies that no other thread is currently processing this task,
   * and buffers the message otherwise.
   * Returns true if this task can be executed now.
   * If false is returned, task should not be executed and
   * scheduler_task_epilogue should not be executed.
   *
   * This may be called outside of the task vlock.
   */
  bool scheduler_task_prologue(lvid_type lvid, const message_type& msg) {
    bool ok_to_run = false;
    vertex_local_data[lvid].lock.lock();
    if (vertex_local_data[lvid].is_running) {
      //this is already running!. buffer it
      if (vertex_local_data[lvid].has_message) {
        vertex_local_data[lvid].buffered_messages += msg;
      } else {
        vertex_local_data[lvid].buffered_messages = msg;
        vertex_local_data[lvid].has_message = true;
      }
    } else {
      vertex_local_data[lvid].is_running = true;
      ok_to_run = true;
    }
    vertex_local_data[lvid].lock.unlock();
    return ok_to_run;
  }

  /**
   *
   * Executes after the scheduled function is run.
   * Releases any messages buffered inside the object.
   * Must be called before releasing the task vlock.
   */
  void scheduler_task_epilogue(lvid_type lvid) {
    vertex_local_data[lvid].lock.lock();
    vertex_local_data[lvid].is_running = false;
    // there was a message. signal it
    if (vertex_local_data[lvid].has_message) {
      internal_signal(graph.l_vertex(lvid),
                      vertex_local_data[lvid].buffered_messages);
      vertex_local_data[lvid].buffered_messages = message_type();
      vertex_local_data[lvid].has_message = false;
    }
    vertex_local_data[lvid].lock.unlock();
  }

  /**
   * Reads a scheduler task from a given queue ID and runs it
   */
  bool exec_scheduler_task(size_t id) {
    GL3TLS::SET_IN_VERTEX_PROGRAM(true);
    context_type context;
    context.engine = this;
    context.thread_id = id;
    lvid_type lvid;
    message_type msg;
    sched_status::status_enum stat =
        scheduler_ptr->get_next(id, lvid, msg);
    // get a task from the scheduler
    // if no task... quit
    if (stat == sched_status::EMPTY) return false;
    // otherwise run the task
    //
    // lock the vertex
    // if this is not the master, we forward it
    if (!graph.l_is_master(lvid)) {
      const procid_t vowner = graph.l_get_vertex_record(lvid).owner;
      rmi.remote_call(vowner,
                      &gl3engine::rpc_signal,
                      graph.global_vid(lvid),
                      msg);
      return true;
    }

    bool ok_to_run = scheduler_task_prologue(lvid, msg);

    // not ok to run. quit
    if (!ok_to_run) return true;

    // now, can I run this task? is it already running?

    while (!vlocks[lvid].try_lock()) qthread_yield();
    GL3TLS::SET_LOCKED_LVID(lvid);
    vertex_type vertex(graph.l_vertex(lvid));
    context.lvid = lvid;


    // logger(LOG_EMPH, "Running vertex %ld", vertex.id());
    active_function(context, vertex, msg);

    scheduler_task_epilogue(lvid);

    programs_completed.inc();
    size_t newhash = get_vertex_data_hash_lvid(lvid);
    // if the hash changed, broadcast
    if (newhash != vertex_local_data[lvid].hash) {
      vertex_local_data[lvid].hash = newhash;
      local_vertex_type lvertex(graph.l_vertex(lvid));
      rmi.remote_call(lvertex.mirrors().begin(), lvertex.mirrors().end(),
                      &engine_type::sync_vdata,
                      lvertex.global_id(),
                      lvertex.data());
    }
    vlocks[lvid].unlock();
    GL3TLS::SET_LOCKED_LVID((lvid_type)(-1));
    return true;
  }


  /**
   * Creates a vthread with a given worker ID
   */
  void subtask_thread_start(size_t id) {
    //std::cout << "Subtask thread " << id << " started\n";
    GL3TLS::SET_IN_VTHREAD_TASK(false);
    size_t ctr = timer::approx_time_millis();
    while(1) {
      rmi.dc().handle_incoming_calls(id % ncpus, ncpus);
      bool haswork = exec_subtasks(id % ncpus);
      // we should yield every so often
      if (timer::approx_time_millis() >= ctr + 100) {
        qthread_yield();
        ctr = timer::approx_time_millis();
      }
      if (!haswork) {
        subtask_thread_alive[id] = 1;
        asm volatile ("" : : : "memory");
        subtask_thread_alive_lock[id].lock();
        // double check
        haswork = exec_subtasks(id % ncpus);
        if (haswork == false) {
          subtask_thread_alive[id] = 0;
          subtask_thread_alive_lock[id].unlock();
          break;
        } else {
          subtask_thread_alive[id] = 2;
          subtask_thread_alive_lock[id].unlock();
        }
      }
    }
    //std::cout << "Subtask thread " << id << " ended\n";
  }


  /**
   * Creates a vthread with a given worker ID
   */
  void vthread_start(size_t id) {
    GL3TLS::SET_IN_VTHREAD_TASK(true);
    size_t ctr = timer::approx_time_millis();
    while(1) {
      bool haswork = exec_scheduler_task(id % ncpus);
      // we should yield every so often
      if (timer::approx_time_millis() >= ctr + 100) {
        qthread_yield();
        ctr = timer::approx_time_millis();
      }
      if (!haswork) {
        active_vthread_count.dec();
        // double check
        haswork = exec_scheduler_task(id % ncpus);
        if (haswork == false) break;
        else active_vthread_count.inc();
      }
      if (programs_completed.value > 0) {
        logger_ontick(1, LOG_EMPH, "programs completed: %ld", programs_completed.value);
      }
    }
  }


/**************************************************************************/
/*                                                                        */
/*                Other Public Thread Management Routines                 */
/*                                                                        */
/**************************************************************************/

 public:
  /**
   * Called by user tasks to yield to another thread.
   */
  void poll(size_t id = (size_t)(-1)) {
    qthread_yield();
  }


  /**
   * Like wait, but non-blocking.
   * Returns true if we are probably near completion,
   * and false otherwise. All machines must call this together,
   * and the function will return the same value to all machines.
   *
   * After calling try_wait(), it is still important to call wait()
   */
  bool try_wait() {
    size_t cvthreads = active_vthread_count.value;
    rmi.all_reduce(cvthreads);
    return cvthreads == 0;
  }



  /**
   * Waits for all user tasks and scheduler tasks to complete.
   */
  execution_status::status_enum wait() {
    qt_end_blocking_action();
    timer ti;
    ti.start();
    while(1) {
      execution_group.join();
      consensus->begin_done_critical_section(0);
      if (!execution_group.empty()) {
        consensus->cancel_critical_section(0);
      } else {
        if (consensus->end_done_critical_section(0)) break;
      }
    }

    size_t ctasks = programs_completed.value;
    rmi.all_reduce(ctasks);
    total_programs_completed = ctasks;


    ctasks = tasks_completed.value;
    rmi.all_reduce(ctasks);
    total_tasks_completed = ctasks;

    consensus->reset();
    rmi.barrier();

    qt_begin_blocking_action();
    return execution_status::TASK_DEPLETION;
  }

  size_t num_updates() const {
    return total_programs_completed;
  }
};

} // graphlab

#undef EXTERNAL_TASK
#undef VTHREAD_TASK
#undef FRESULT
#undef REMOVE_CONST_REF

#include <graphlab/macros_undef.hpp>
#endif
