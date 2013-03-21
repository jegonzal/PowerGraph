#ifndef GRAPHLAB_GL3ENGINE_HPP
#define GRAPHLAB_GL3ENGINE_HPP
#include <vector>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/functional/hash.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

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



  /// \internal \brief The base type of all schedulers
  typedef ischeduler<message_type> ischeduler_type;
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
  size_t num_vthreads;
  size_t ncpus;
  graph_type& graph;
  task_descriptor_type* task_types[256];
  bool finished;
  float engine_runtime;

  vertex_program_type active_function;
  std::vector<inplace_lf_queue2<task>* > local_tasks;
  std::vector<size_t> vdata_hash;
  atomic<size_t> programs_completed;
  atomic<size_t> tasks_completed;

  size_t total_programs_completed;
  size_t total_tasks_completed;

  atomic<size_t> active_vthread_count;
  atomic<size_t> thread_counter;

  std::vector<mutex> worker_mutex;
  //! The scheduler
  ischeduler_type* scheduler_ptr;

  bool block_vthread_launches;
  async_consensus* consensus;

  static const size_t NUM_DHTS = 256;
  mutex dht_lock[NUM_DHTS];
  boost::unordered_map<size_t, any> dht[NUM_DHTS]; // a large dht

  /**
   * \brief The vertex locks protect access to vertex specific
   * data-structures including
   * \ref graphlab::synchronous_engine::gather_accum
   * and \ref graphlab::synchronous_engine::messages.
   */
  std::vector<simple_spinlock> vlocks;


  /**
   * \brief The elocks protect individual edges during gather and
   * scatter.  Technically there is a potential race since gather
   * and scatter can modify edge values and can overlap.  The edge
   * lock ensures that only one gather or scatter occurs on an edge
   * at a time.
   */
  std::vector<simple_spinlock> elocks;

  bool rpc_handlers_disabled;

  graphlab::qthread_group execution_group;


  friend struct gl3task_descriptor<GraphType, engine_type>;
  template <typename, typename, typename>
  friend struct map_reduce_neighbors_task_descriptor;
  friend struct broadcast_task_descriptor<GraphType, engine_type>;
  friend struct dht_gather_task_descriptor<GraphType, engine_type>;
  friend struct dht_scatter_task_descriptor<GraphType, engine_type>;

  friend struct gl3context<engine_type>;
  context_type base_context;

 public:
  gl3engine(distributed_control& dc, graph_type& graph,
            const graphlab_options& opts = graphlab_options()):
      rmi(dc, this), graph(graph){
    rmi.barrier();
    num_vthreads = 1000;
    ncpus = opts.get_ncpus();
    worker_mutex.resize(ncpus);
    rpc_handlers_disabled = false;
    block_vthread_launches = false;


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
    vdata_hash.resize(graph.num_local_vertices());

    // construct the queues
    local_tasks.resize(ncpus);
    for (size_t i = 0;i < local_tasks.size(); ++i) {
      local_tasks[i] = new inplace_lf_queue2<task>();
    }

    // make default registrations
    task_types[GL3_BROADCAST_TASK_ID] = new broadcast_task_descriptor<GraphType, engine_type>();
    task_types[GL3_DHT_GATHER_TASK_ID] = new dht_gather_task_descriptor<GraphType, engine_type>();
    //
    // pre init
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)graph.num_local_vertices(); ++i) {
      vdata_hash[i] = get_vertex_data_hash_lvid((lvid_type)i);
    }
    scheduler_ptr->start();

    //reset counters
    finished = false;
    programs_completed = 0;
    tasks_completed = 0;
    total_programs_completed = 0;
    total_tasks_completed = 0;
    active_vthread_count = 0;



    graphlab::qthread_tools::init(ncpus, 128 * 1024);

    base_context.engine = this;
    base_context.lvid = (lvid_type)(-1);
    rmi.full_barrier();
  }


  template <typename T>
  struct create_map_reduce_task_impl{
    typedef boost::function<T (const vertex_type&,
                               edge_type&,
                               const vertex_type&)>  full_map_function_type;

    typedef boost::function<T (const vertex_type&)> basic_map_function_type;


    template <typename MapFn>
    static typename boost::enable_if_c<boost::is_same<MapFn, full_map_function_type>::value,
                    task_descriptor_type*>::type
      create(full_map_function_type mapfn,
             boost::function<void (T&, const T&)> combinefn,
             T zero = T()) {
      return new map_reduce_neighbors_task_descriptor<GraphType, engine_type, T>(mapfn, combinefn, zero);
    }

    static T simple_map_function_dispatch(boost::function<T (const vertex_type&)> mapfn,
                                   const vertex_type&,
                                   edge_type&,
                                   const vertex_type& other) {
      return mapfn(other);
    }

    template <typename MapFn>
    static typename boost::enable_if_c<boost::is_same<MapFn, basic_map_function_type>::value,
                    task_descriptor_type*>::type
      create(basic_map_function_type mapfn,
             boost::function<void (T&, const T&)> combinefn,
             T zero = T()) {
      return new map_reduce_neighbors_task_descriptor<GraphType, engine_type, T>(
          boost::bind(simple_map_function_dispatch, mapfn, _1, _2, _3), combinefn, zero);
    }
  };

  template <typename MapFn, typename CombineFn>
  void register_map_reduce(size_t id,
                           MapFn mapfn,
                           CombineFn combinefn,
                           FRESULT(MapFn) zero = FRESULT(MapFn) () ) {
    task_types[id] = create_map_reduce_task_impl<FRESULT(MapFn)>
        ::template create<NORMALIZE_FUNCTION(MapFn)>(mapfn, combinefn, zero);
    rmi.barrier();
  }

  void register_dht_scatter(size_t id,
                           boost::function<void (any&, const any&)> scatter_fn) {
    dht_scatter_task_descriptor<GraphType, engine_type>* desc =
        new dht_scatter_task_descriptor<GraphType, engine_type>;
    desc->scatter_fn = scatter_fn;
    task_types[id] = desc;
    rmi.barrier();
  }




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


  context_type& get_context() {
    return base_context;
  }

  void set_vertex_program(vertex_program_type uf) {
    active_function = uf;
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

  void internal_signal(const vertex_type& vtx,
                       const message_type& message = message_type()) {
    internal_schedule(vtx.local_id(), message);
  } // end of schedule


  void internal_schedule(const lvid_type lvid,
                         const message_type& message) {
    scheduler_ptr->schedule(lvid, message);
    size_t target_queue = thread_counter++;

    if (!block_vthread_launches && active_vthread_count < num_vthreads) {
      active_vthread_count.inc();
      execution_group.launch(boost::bind(&gl3engine::vthread_start,
                                         this,
                                         target_queue));
      consensus->cancel();
    }
  }

  void launch_other_task(boost::function<void (context_type&)> fn) {
    size_t target_queue = thread_counter++;
    execution_group.launch(boost::bind(&gl3engine::task_start,
                                       this,
                                       fn, 
                                       target_queue));
  }

  void task_start(boost::function<void (context_type&)> fn, size_t id) {
    context_type context;
    context.engine = this;
    context.lvid = lvid_type(-1);
    context.thread_id = id;    
    fn(context);
  }


  struct __attribute((__may_alias__)) future_combiner {
    any param;
    any result;
    atomic<procid_t> count_down;
    unsigned char task_id;
    simple_spinlock lock;
  };

  // a task which is not associated with any vertices
  any spawn_task(unsigned char task_id,
                 std::vector<procid_t>& target_machines,
                 const std::vector<any>& task_param,
                 bool no_reply = false) {
    future_combiner combiner;
    combiner.count_down = target_machines.size();
    combiner.task_id = task_id;
    combiner.param = task_param;

    conditional_serialize<vertex_data_type> cs;
    size_t cb = reinterpret_cast<size_t>(&combiner);
    if (no_reply) cb = 0;

    bool has_self_task = false;
    size_t self_task_id = 0;
    // if I am in the target_machines list, move myself to the front
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
      while(combiner.count_down != 0) {
        qthread_yield();
      }
    }
    return combiner.result;
  }

  any spawn_task(lvid_type lvid,
                 unsigned char task_id,
                 const any& task_param,
                 bool no_reply = false) {
    ASSERT_TRUE(graph.l_is_master(lvid));
    future_combiner combiner;
    combiner.count_down = graph.l_vertex(lvid).num_mirrors() + 1;
    combiner.task_id = task_id;
    combiner.param = task_param;
    local_vertex_type lvertex(graph.l_vertex(lvid));
    /*
    logstream(LOG_EMPH) << "Creating Subtask type "<< (int)task_id
                        << " on vertex " << graph.l_vertex(lvid).global_id() << " Handle " << combiner
                        << "\n";
    */
    size_t newhash = get_vertex_data_hash_lvid(lvid);
    conditional_serialize<vertex_data_type> cs;
    if (newhash != vdata_hash[lvid]) {
      cs.hasval = true;
      cs.val = lvertex.data();
      vdata_hash[lvid] = newhash;
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
    // unlock the task locks so we don't dead-lock with the subtask
    vlocks[lvid].unlock();
    any ret = task_types[task_id]->exec(graph, lvertex.global_id(), task_param, this);
    if (!no_reply) {
      task_reply(&combiner, ret);
      while(combiner.count_down != 0) {
        qthread_yield();
      }
    }

    while (!vlocks[lvid].try_lock()) qthread_yield();
    return combiner.result;
  }
  void task_reply_rpc(size_t handle, any& val) {
    future_combiner* combiner = reinterpret_cast<future_combiner*>(handle);
    task_reply(combiner, val);
  }

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

  void rpc_receive_task(unsigned char task_id,
                        vertex_id_type vid,
                        conditional_serialize<vertex_data_type>& vdata,
                        const any& param,
                        procid_t caller,
                        size_t handle) {
    //logstream(LOG_EMPH) << "Receiving subtask on handle " << (void*)(handle) << "\n";
    if (vid != vertex_id_type(-1)) {
      lvid_type lvid = graph.local_vid(vid);
      ASSERT_FALSE(graph.l_is_master(lvid));
      if (vdata.hasval) {
        vlocks[lvid].lock();
        graph.l_vertex(lvid).data() = vdata.val;
        vlocks[lvid].unlock();
      }
    }
    task* t = new task;
    t->origin = caller;
    t->param = param;
    t->task_id = task_id;
    t->handle = handle;
    t->vid = vid;
    size_t target_queue = (thread_counter++) % ncpus;
    local_tasks[target_queue]->enqueue(t);

    if (!block_vthread_launches && active_vthread_count < num_vthreads) {
      active_vthread_count.inc();
      execution_group.launch(boost::bind(&gl3engine::vthread_start,
                                         this,
                                         target_queue));
      consensus->cancel();
    }

  }

  size_t get_vertex_data_hash_lvid(lvid_type lvid) {
    return boost::hash_value(graph.l_vertex(lvid).data());
  }

  void sync_vdata(vertex_id_type vid,
                  vertex_data_type& vdata) {
    vertex_type vertex(graph.vertex(vid));
    lvid_type lvid = vertex.local_id();
    vlocks[lvid].lock();
    vertex.data() = vdata;
    vlocks[lvid].unlock();
  }


  bool exec_scheduler_task(size_t id) {
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
    while (!vlocks[lvid].try_lock()) qthread_yield();
    vertex_type vertex(graph.l_vertex(lvid));
    context.lvid = lvid;

    // logger(LOG_EMPH, "Running vertex %ld", vertex.id());
    active_function(context, vertex, msg);
    programs_completed.inc();
    size_t newhash = get_vertex_data_hash_lvid(lvid);
    // if the hash changed, broadcast
    if (newhash != vdata_hash[lvid]) {
      vdata_hash[lvid] = newhash;
      local_vertex_type lvertex(graph.l_vertex(lvid));
      rmi.remote_call(lvertex.mirrors().begin(), lvertex.mirrors().end(),
                      &engine_type::sync_vdata,
                      lvertex.global_id(),
                      lvertex.data());
    }
    vlocks[lvid].unlock();
    return true;
  }

  void poll(size_t id) {
    rmi.dc().handle_incoming_calls(id % ncpus, ncpus);
    exec_subtasks(id % ncpus);
    exec_scheduler_task(id % ncpus);
    qthread_yield();
  }

  void vthread_start(size_t id) {
    while(!finished) {
      rmi.dc().handle_incoming_calls(id % ncpus, ncpus);
      bool haswork = exec_subtasks(id % ncpus);
      haswork |= exec_scheduler_task(id % ncpus);
      if (!haswork) {
        active_vthread_count.dec();
        // double check
        haswork = exec_subtasks(id % ncpus);
        haswork |= exec_scheduler_task(id % ncpus);
        if (haswork == false) break;
        else active_vthread_count.inc();
      }
      qthread_yield();
      if (programs_completed.value > 0) {
        logger_ontick(1, LOG_EMPH, "programs completed: %ld", programs_completed.value);
      }
    }
  }

  void ping() {
  }

  atomic<size_t> pingid;

  bool exec_subtasks(size_t worker) {
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

  execution_status::status_enum start() {
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
    engine_runtime = ti.current_time();

    size_t ctasks = programs_completed.value;
    rmi.all_reduce(ctasks);
    total_programs_completed = ctasks;


    ctasks = tasks_completed.value;
    rmi.all_reduce(ctasks);
    total_tasks_completed = ctasks;

    consensus->reset();
    rmi.barrier();
    return execution_status::TASK_DEPLETION;
  }

  size_t num_updates() const {
    return programs_completed;
  }

  float elapsed_seconds() const  {
    return engine_runtime;
  }


};

} // graphlab

#undef FRESULT
#undef REMOVE_CONST_REF

#include <graphlab/macros_undef.hpp>
#endif
