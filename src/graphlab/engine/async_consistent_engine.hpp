/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */




#ifndef GRAPHLAB_ASYNC_CONSISTENT_ENGINE
#define GRAPHLAB_ASYNC_CONSISTENT_ENGINE

#include <deque>
#include <boost/bind.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>
#include <graphlab/vertex_program/ivertex_program.hpp>
#include <graphlab/vertex_program/icontext.hpp>
#include <graphlab/vertex_program/context.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/distributed_chandy_misra.hpp>

#include <graphlab/util/tracepoint.hpp>
#include <graphlab/util/memory_info.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/rpc/async_consensus.hpp>
#include <graphlab/aggregation/distributed_aggregator.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {
  
  
  /**
   * \ingroup engines
   *
   * \brief The asynchronous consistent engine executed vertex programs
   * asynchronously but ensures mutual exclusion such that adjacent vertices
   * are never executed simultaneously.
   *
   * \tparam VertexProgram
   * The user defined vertex program type which should implement the
   * \ref graphlab::ivertex_program interface.
   *
   * ### Execution Semantics
   * 
   * On start() the \ref graphlab::ivertex_program::init function is invoked
   * on all vertex programs in parallel to initialize the vertex program,
   * vertex data, and possibly signal vertices.
   *
   * After which, the engine spawns a collection of threads where each thread
   * individually performs the following tasks:
   * \li Extract a message from the scheduler. 
   * \li Perform distributed lock acquisition on the vertex which is supposed
   * to receive the message. The lock system enforces that no neighboring
   * vertex is executing at the same time. The implementation is based
   * on the Chandy-Misra solution to the dining philosophers problem. 
   * (Chandy, K.M.; Misra, J. (1984). The Drinking Philosophers Problem.
   *  ACM Trans. Program. Lang. Syst)
   * \li Once lock acquisition is complete,
   *  \ref graphlab::ivertex_program::recv_message is called on the vertex
   * program. As an optimization, any messages sent to this vertex
   * before completion of lock acquisition is merged into original message
   * extracted from the scheduler.
   * \li Execute the gather on the vertex program by invoking
   * the user defined \ref graphlab::ivertex_program::gather function
   * on the edge direction returned by the
   * \ref graphlab::ivertex_program::gather_edges function.  The gather
   * functions can modify edge data but cannot modify the vertex
   * program or vertex data and can be executed on multiple
   * edges in parallel. 
   * * \li Execute the apply function on the vertex-program by
   * invoking the user defined \ref graphlab::ivertex_program::apply
   * function passing the sum of the gather functions.  If \ref
   * graphlab::ivertex_program::gather_edges returns no edges then
   * the default gather value is passed to apply.  The apply function
   * can modify the vertex program and vertex data.
   * \li Execute the scatter on the vertex program by invoking
   * the user defined \ref graphlab::ivertex_program::scatter function
   * on the edge direction returned by the
   * \ref graphlab::ivertex_program::scatter_edges function.  The scatter
   * functions can modify edge data but cannot modify the vertex
   * program or vertex data and can be executed on multiple
   * edges in parallel.
   * \li Release all locks acquired in the lock acquisition stage,
   * and repeat until the scheduler is empty.
   *
   * The engine threads multiplexes the above procedure through a secondary
   * internal queue, allowing an arbitrary large number of vertices to
   * begin processing at the same time.
   *
   * ### Construction
   *
   * The asynchronous consistent engine is constructed by passing in a
   * \ref graphlab::distributed_control object which manages coordination
   * between engine threads and a \ref graphlab::distributed_graph object
   * which is the graph on which the engine should be run.  The graph should
   * already be populated and cannot change after the engine is constructed.
   * In the distributed setting all program instances (running on each machine)
   * should construct an instance of the engine at the same time.
   *
   * Computation is initiated by signaling vertices using either
   * \ref graphlab::async_consistent_engine::signal or
   * \ref graphlab::async_consistent_engine::signal_all.  In either case all
   * machines should invoke signal or signal all at the same time.  Finally,
   * computation is initiated by calling the
   * \ref graphlab::async_consistent_engine::start function.
   *
   * ### Example Usage
   *
   * The following is a simple example demonstrating how to use the engine:
   * \code
   * #include <graphlab.hpp>
   *
   * struct vertex_data {
   *   // code
   * };
   * struct edge_data {
   *   // code
   * };
   * typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
   * typedef float gather_type;
   * struct pagerank_vprog :
   *   public graphlab::ivertex_program<graph_type, gather_type> {
   *   // code
   * };
   *
   * int main(int argc, char** argv) {
   *   // Initialize control plain using mpi
   *   graphlab::mpi_tools::init(argc, argv);
   *   graphlab::distributed_control dc;
   *   // Parse command line options
   *   graphlab::command_line_options clopts("PageRank algorithm.");
   *   std::string graph_dir;
   *   clopts.attach_option("graph", &graph_dir, graph_dir,
   *                        "The graph file.");
   *   if(!clopts.parse(argc, argv)) {
   *     std::cout << "Error in parsing arguments." << std::endl;
   *     return EXIT_FAILURE;
   *   }
   *   graph_type graph(dc, clopts);
   *   graph.load_format(graph_dir, "tsv");
   *   graph.finalize();
   *   std::cout << "#vertices: " << graph.num_vertices()
   *             << " #edges:" << graph.num_edges() << std::endl;
   *   graphlab::async_consistent_engine<pagerank_vprog> engine(dc, graph, clopts);
   *   engine.signal_all();
   *   engine.start();
   *   std::cout << "Runtime: " << engine.elapsed_time();
   *   graphlab::mpi_tools::finalize();
   * }
   * \endcode
   *
   * \see graphlab::synchronous_engine
   */
   */
  template<typename VertexProgram>
  class async_consistent_engine: public iengine<VertexProgram> {
      
  public:
    /**
     * \brief The user defined vertex program type. Equivalent to the
     * VertexProgram template argument.
     *
     * The user defined vertex program type which should implement the
     * \ref graphlab::ivertex_program interface.
     */
    typedef VertexProgram vertex_program_type;

    /**
     * \brief The user defined type returned by the gather function.
     *
     * The gather type is defined in the \ref graphlab::ivertex_program
     * interface and is the value returned by the
     * \ref graphlab::ivertex_program::gather function.  The
     * gather type must have an <code>operator+=(const gather_type&
     * other)</code> function and must be \ref serializable.
     */
    typedef typename VertexProgram::gather_type gather_type;

    /**
     * \brief The user defined message type used to signal neighboring
     * vertex programs.
     *
     * The message type is defined in the \ref graphlab::ivertex_program
     * interface and used in the call to \ref graphlab::icontext::signal.
     * The message type must have an
     * <code>operator+=(const gather_type& other)</code> function and
     * must be \ref serializable.
     */
    typedef typename VertexProgram::message_type message_type;

    /**
     * \brief The type of data associated with each vertex in the graph
     *
     * The vertex data type must be \ref serializable.
     */
    typedef typename VertexProgram::vertex_data_type vertex_data_type;

    /**
     * \brief The type of data associated with each edge in the graph
     *
     * The edge data type must be \ref serializable.
     */
    typedef typename VertexProgram::edge_data_type edge_data_type;

    /**
     * \brief The type of graph supported by this vertex program
     *
     * See graphlab::distributed_graph
     */
    typedef typename VertexProgram::graph_type graph_type;

     /**
     * \brief The type used to represent a vertex in the graph.
     * See \ref graphlab::distributed_graph::vertex_type for details
     *
     * The vertex type contains the function
     * \ref graphlab::distributed_graph::vertex_type::data which
     * returns a reference to the vertex data as well as other functions
     * like \ref graphlab::distributed_graph::vertex_type::num_in_edges
     * which returns the number of in edges.
     *
     */
    typedef typename graph_type::vertex_type          vertex_type;

    /**
     * \brief The type used to represent an edge in the graph.
     * See \ref graphlab::distributed_graph::edge_type for details.
     *
     * The edge type contains the function
     * \ref graphlab::distributed_graph::edge_type::data which returns a
     * reference to the edge data.  In addition the edge type contains
     * the function \ref graphlab::distributed_graph::edge_type::source and
     * \ref graphlab::distributed_graph::edge_type::target.
     *
     */
    typedef typename graph_type::edge_type            edge_type;

    /**
     * \brief The type of the callback interface passed by the engine to vertex
     * programs.  See \ref graphlab::icontext for details.
     *
     * The context callback is passed to the vertex program functions and is
     * used to signal other vertices, get the current iteration, and access
     * information about the engine.
     */
    typedef typename context<async_consistent_engine>::icontext_type icontext_type;
    
  private:
    /**
     * \internal
     * \brief A wrapper around the gather_type to automatically manage
     * an "empty" state.
     *
     * \see conditional_addition_wrapper
     */
    typedef conditional_addition_wrapper<gather_type> conditional_gather_type;

    /// \internal \brief The base type of all schedulers
    typedef ischeduler<message_type> ischeduler_type;

    /** \internal
     * \brief The true type of the callback context interface which
     * implements icontext. \see graphlab::icontext graphlab::context
     */
    typedef context<async_consistent_engine> context_type;

    // context needs access to internal functions
    friend class context<async_consistent_engine>;

    /// \internal \brief The type used to refer to vertices in the local graph
    typedef typename graph_type::local_vertex_type    local_vertex_type;
    /// \internal \brief The type used to refer to edges in the local graph
    typedef typename graph_type::local_edge_type      local_edge_type;
    /// \internal \brief The type used to refer to vertex IDs in the local graph
    typedef typename graph_type::lvid_type            lvid_type;

    /// \internal \brief The type of the current engine instantiation
    typedef async_consistent_engine<VertexProgram> engine_type;

    /**
     * \internal
     * \brief States in the vertex state machine
     *
     * On a master vertex, it transitions sequentially
     * between the following 5 states
     * \li \c NONE Nothing is going on and nothing is scheduled on this
     * vertex. The NONE state transits into the LOCKING state when
     * a message is popped from the scheduler destined for this vertex.
     * The popped message is stored in
     * <code>vertex_state[lvid].current_message</code> and a call to the lock
     * implementation \ref graphlab::distributed_chandy_misra is made to
     * acquire exclusive locks on this vertex ( master_broadcast_locking(),
     * rpc_begin_locking() ) .
     * \li \c LOCKING  At this point, we have an active message
     * stored in <code>vertex_state[lvid].current_message</code>, and we
     * are waiting for lock acquisition to complete. When lock acquisition is
     * complete, we will be signalled by
     * \ref graphlab::distributed_chandy_misra via the lock callback
     * lock_ready() . Once a vertex is in a locking state, if the vertex
     * is ever signaled, the scheduler is bypassed and the signalled message
     * is merged into <code>vertex_state[lvid].current_message</code>
     * directly. When locks acquisition is complete, we broadcast a
     * remote_call to all mirrors to begin gathering
     * ( master_broadcast_gathering() ). The master vertex transits into
     * the GATHERING state, while all mirrored vertices transit into the
     * MIRROR_GATHERING state ( rpc_begin_gathering() ).
     * \li \c GATHERING/MIRROR_GATHERING Upon entry into the GATHERING state,
     * an internal task corresponding to this vertex is put on the intenal
     * task queue. When the task is popped and executed, it will complete
     * gathering on the local graph, and send the partial gathered result
     * back to the master. If the vertex is in MIRROR_GATHERING, it is a
     * mirror, and it transits back to NONE
     * ( rpc_gather_complete(), decrement_gather_counter() ) .
     * If the vertex is in GATHERING, it stays in GATHERING until all partial
     * gathers are received. After which it transits into APPLYING.
     * \li \c APPLYING Upon entry into the APPLYING state, an internal task
     *corresponding to this veretx is put on the internal task queue. When
     the task is popped and executed, the apply is performed on the master
     vertex. It then immediately transits into SCATTERING 
     * \li \c SCATTERING No internal task is generated for scattering. On
     * the master vertex, scattering is performed immediately after applying.
     * In SCATTERING, a partial scatter is performed on the local (master)
     * vertex, and all mirrors are called ( master_broadcast_scattering() )
     * to transit into MIRROR_SCATTERING. A partial lock release request is
     * also issued on the master's subgraph. The master vertex than
     * transits to the NONE state.
     * \li \c MIRROR_SCATTERING When a mirror enters the MIRROR_SCATTERING
     * state ( rpc_begin_scattering() ) it inserts a task into the internal
     * queue. When the task is popped and executed, it will perform
     * a partial scatter on the local (mirror) vertex, then perform
     * a partial lock release request on the mirror's subgraph.
     * \li \c MIRROR_SCATTERING_AND_NEXT_LOCKING Due to the distributed
     * nature, there is an unusual condition in which the master completes
     * scattering a new task has been scheduled on the same vertex, but
     * mirrors have not completed scattering yet. In other words, the master
     * vertex is in the LOCKING state while some mirrors are still in
     * MIRROR_SCATTERING. This state is unique on the mirror to handle
     * this case. If a lock request is received while the vertex is still in
     * MIRROR_SCATTERING, the vertex transitions into this state. This state
     * is equivalent to MIRROR_SCATTERING except that after completion
     * of the scatter, it performs the partial lock release, then
     * immediately re-requests the locks required on the local subgraph.
     */
    enum vertex_execution_state {
      NONE = 0,
      LOCKING,     // state on owner
      GATHERING,   // state on owner
      APPLYING,    // state on owner
      SCATTERING,  // state on owner
      MIRROR_GATHERING, // state on mirror
      MIRROR_SCATTERING, // state on mirror
      MIRROR_SCATTERING_AND_NEXT_LOCKING, // state on mirror
    }; // end of vertex execution state


    /**
     * \internal
     * The state machine + additional state maintained for each
     * vetex and mirrors. 
     */
    struct vertex_state {
      /**
       * This stores the active instance of the vertex program
       */
      vertex_program_type vertex_program;
      /**
       * This stores the current active message being executed.
       */
      message_type current_message;
      /**
       * This is temporary storage for the partial gathers on each mirror,
       * and the accumulated gather on each machine.
       */
      conditional_gather_type combined_gather;
      /**
       * This is used to count down the number of gathers completed.
       * It starts off at #mirrors + 1, and is decremented everytime
       * one partial gather is received. The decrement is performed
       * by decrement_gather_counter() . When the counter hits 0,
       * the master vertex transits into APPLYING
       */
      uint32_t apply_count_down;
      /** This condition happens when a vertex is scheduled (popped from
       * the scheduler) while the vertex is still running. In which case
       * the message is placed back in the scheduler (without scheduling it),
       * and the hasnext flag is set. When the vertex finishes execution,
       * the scheduler is signaled to schedule the blocked vertex.
       * This may only be set on a master vertex
       */
      bool hasnext;
      /// A 1 byte lock protecting the vertex state
      simple_spinlock slock;
      /// State of each vertex. \ref vertex_execution_state for details
      /// on the meaning of each state
      vertex_execution_state state;

      /// initializes the vertex state to default values.
      vertex_state(): current_message(message_type()),
                      apply_count_down(0),
                      hasnext(false),
                      state(NONE) { }
      
      /// Acquires a lock on the vertex state
      void lock() {
        slock.lock();
      }
      /// releases a lock on the vertex state
      void unlock() {
        slock.unlock();
      }
    }; // end of vertex_state
    
    /**
     * \internal
     * This stores the internal scheduler for each thread.
     * The representation is a single deque of pending vertices
     * which is continuually appended to.
     * If the vertex ID in the deque is < #local vertices, it is a task on
     * the vertex in question. The actual task will depend on the state
     * the vertex is currently in.
     * However, if the vertex ID is >= #local vertices, the task is an
     * aggregation. The aggregation key can be found in
     * aggregate_id_to_key[-id]
      */
    struct thread_local_data {
      mutex lock;
      size_t npending;
      std::deque<lvid_type> pending_vertices;
      thread_local_data() : npending(0) { }       
      void add_task(lvid_type v) {
        lock.lock();
        ++npending;
        pending_vertices.push_back(v);
        lock.unlock();
      }
      void add_task_priority(lvid_type v) {
        lock.lock();
        ++npending;
        pending_vertices.push_front(v);
        lock.unlock();
      }
      bool get_task(std::deque<lvid_type> &v) {
        v = std::deque<lvid_type>();
        lock.lock();
        if (npending == 0) { lock.unlock(); return false; }
        npending = 0;
        v.swap(pending_vertices);
        lock.unlock();
        return true;
      }
    }; // end of thread local data

    /// The RPC interface
    dc_dist_object<async_consistent_engine<VertexProgram> > rmi;

    /// A reference to the active graph
    graph_type& graph;

    /// A pointer to the lock implementation
    distributed_chandy_misra<graph_type>* cmlocks;

    /// Engine threads.
    thread_group thrgroup;
    
    //! The scheduler
    ischeduler_type* scheduler_ptr;

    /// vector of vertex states. vstate[lvid] contains the vertex state
    /// for local VID lvid.
    std::vector<vertex_state> vstate;
    /// Used if cache is enabled. ("use_cache=true"). Contains the result
    /// of partial gathers
    std::vector<conditional_gather_type> cache;
    
    typedef distributed_aggregator<graph_type, icontext_type> aggregator_type;
    aggregator_type aggregator;

    /// Number of engine threads
    size_t ncpus;
    /// set to true if engine is started
    bool started;
    /// A pointer to the distributed consensus object
    async_consensus* consensus;

    /// Thread local store. Used to hold the internal queus
    std::vector<thread_local_data> thrlocal;

    // Various counters.
    atomic<uint64_t> joined_messages;
    atomic<uint64_t> blocked_issues; // issued messages which either
                                     // 1: cannot start and have to be 
                                     //    reinjected into the 
                                     //    scheduler.
                                     // 2: issued but is combined into a
                                     //    message which is currently locking
    atomic<uint64_t> issued_messages;
    atomic<uint64_t> programs_executed;

    /// engine option. Maximum proportion of clean forks allowed
    float max_clean_fraction;
    /// Equals to #edges * max_clean_fraction
    /// Maximum number of clean forks allowed
    size_t max_clean_forks;
    /// Defaults to (-1), defines a timeout
    size_t timed_termination;
    /// engine option Sets to true if caching is enabled
    bool use_cache;
    /// Time when engine is started
    float engine_start_time;
    /// True when a force stop is triggered (possibly via a timeout)
    bool force_stop;
    
    graphlab_options opts_copy; // local copy of options to pass to 
                                // scheduler construction
    
    execution_status::status_enum termination_reason; 

    PERMANENT_DECLARE_DIST_EVENT_LOG(eventlog);
    DECLARE_TRACER(disteng_eval_sched_task);
    DECLARE_TRACER(disteng_chandy_misra);
    DECLARE_TRACER(disteng_init_gathering); 
    DECLARE_TRACER(disteng_init_scattering);
    DECLARE_TRACER(disteng_evalfac);
    DECLARE_TRACER(disteng_internal_task_queue);

 
    
    inline void ASSERT_I_AM_OWNER(const lvid_type lvid) const {
      ASSERT_EQ(graph.l_get_vertex_record(lvid).owner, rmi.procid());
    }
    inline void ASSERT_I_AM_NOT_OWNER(const lvid_type lvid) const {
      ASSERT_NE(graph.l_get_vertex_record(lvid).owner, rmi.procid());
    }

    enum {
      SCHEDULE_EVENT = 0,
      UPDATE_EVENT = 1,
      WORK_ISSUED_EVENT = 2,
      INTERNAL_TASK_EVENT = 3,
      NO_WORK_EVENT = 4,
      SCHEDULE_FROM_REMOTE_EVENT = 5,
      USER_OP_EVENT = 6,
      ENGINE_START_EVENT = 7,
      ENGINE_STOP_EVENT = 8
    };



  public:

    /**
     * Constructs an asynchronous consistent distributed engine.
     * The number of threads to create are read from 
     * opts::get_ncpus(). The scheduler to construct is read from
     * graphlab_options::get_scheduler_type(). The default scheduler
     * is the queued_fifo scheduler. For details on the scheduler types
     * \see scheduler_types
     *
     * Valid engine options (graphlab_options::get_engine_args()):
     * \arg \c max_clean_fraction The maximum proportion of edges which
     * can be locked at any one time. (This is a simplification of the actual
     * clean/dirty concept in the Chandy Misra locking algorithm used in this
     * implementation.
     * \arg \c timed_termination Maximum time in seconds the engine will run
     * for. The actual runtime may be marginally greater as the engine 
     * waits for all threads and processes to flush all active tasks before
     * returning.
     * \arg \c use_cache If set to true, partial gathers are cached.
     * See \ref gather_caching to understand the behavior of the
     * gather caching model and how it may be used to accelerate program
     * performance.
     * 
     * \param dc Distributed controller to associate with
     * \param graph The graph to schedule over. The graph must be fully
     *              constructed and finalized.
     * \param opts A graphlab_options object containing options and parameters
     *             for the scheduler and the engine.
     */
    async_consistent_engine(distributed_control &dc,
                            graph_type& graph, 
                            const graphlab_options& opts) : 
        rmi(dc, this), graph(graph), scheduler_ptr(NULL),
        aggregator(dc, graph, new context_type(*this, graph)), started(false),
        engine_start_time(timer::approx_time_seconds()), force_stop(false),
        vdata_exchange(dc),thread_barrier(opts.get_ncpus()) {
      rmi.barrier();

      // set default values
      max_clean_fraction = 1.0;
      max_clean_forks = (size_t)(-1);
      timed_termination = (size_t)(-1);
      use_cache = false;
      termination_reason = execution_status::UNSET;
      set_options(opts);
      
#ifdef USE_EVENT_LOG
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 3000, 
                                dist_event_log::RATE_BAR);
#else
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 3000, 
                                dist_event_log::LOG_FILE);
#endif

      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, SCHEDULE_EVENT, "Schedule");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, UPDATE_EVENT, "Updates");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, WORK_ISSUED_EVENT, "Work Issued");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, INTERNAL_TASK_EVENT, "Internal");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, NO_WORK_EVENT, "No Work");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, SCHEDULE_FROM_REMOTE_EVENT, "Remote Schedule");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, USER_OP_EVENT, "User Ops");
      PERMANENT_ADD_IMMEDIATE_DIST_EVENT_TYPE(eventlog, ENGINE_START_EVENT, "Engine Start");
      PERMANENT_ADD_IMMEDIATE_DIST_EVENT_TYPE(eventlog, ENGINE_STOP_EVENT, "Engine Stop");
      

      INITIALIZE_TRACER(disteng_eval_sched_task, 
                        "distributed_engine: Evaluate Scheduled Task");
      INITIALIZE_TRACER(disteng_init_gathering,
                        "distributed_engine: Initialize Gather");
      INITIALIZE_TRACER(disteng_init_scattering,
                        "distributed_engine: Initialize Scattering");
      INITIALIZE_TRACER(disteng_evalfac,
                      "distributed_engine: Time in Factorized Update user code");
      INITIALIZE_TRACER(disteng_internal_task_queue,
                      "distributed_engine: Time in Internal Task Queue");
      INITIALIZE_TRACER(disteng_chandy_misra,
                      "distributed_engine: Time in Chandy Misra");
      initialize();
      rmi.barrier();
    }
    
  private:

    /**
     * \internal
     * Configures the engine with the provided options.
     * The number of threads to create are read from 
     * opts::get_ncpus(). The scheduler to construct is read from
     * graphlab_options::get_scheduler_type(). The default scheduler
     * is the queued_fifo scheduler. For details on the scheduler types
     * \see scheduler_types
     *
     * Valid engine options (graphlab_options::get_engine_args()):
     * \arg \c max_clean_fraction The maximum proportion of edges which
     * can be locked at any one time. (This is a simplification of the actual
     * clean/dirty concept in the Chandy Misra locking algorithm used in this
     * implementation.
     * \arg \c timed_termination Maximum time in seconds the engine will run
     * for. The actual runtime may be marginally greater as the engine 
     * waits for all threads and processes to flush all active tasks before
     * returning.
     */
    void set_options(const graphlab_options& opts) {
      rmi.barrier();
      ncpus = opts.get_ncpus();
      ASSERT_GT(ncpus, 0);
      thread_barrier.resize_unsafe(opts.get_ncpus());
      std::vector<std::string> keys = opts.get_engine_args().get_option_keys();
      foreach(std::string opt, keys) {
        if (opt == "max_clean_fraction") {
          opts.get_engine_args().get_option("max_clean_fraction", max_clean_fraction);
        } else if (opt == "timed_termination") {
          opts.get_engine_args().get_option("timed_termination", timed_termination);
        } else if (opt == "use_cache") {
          opts.get_engine_args().get_option("use_cache", use_cache);
        } else {
          logstream(LOG_FATAL) << "Unexpected Engine Option: " << opt << std::endl;
        }
      }
      opts_copy = opts;
      // set a default scheduler if none
      if (opts_copy.get_scheduler_type() == "") {
        opts_copy.set_scheduler_type("queued_fifo");
      }
      rmi.barrier();
      
    }

    /**
     * \internal
     * Initializes the engine with respect to the associated graph.
     * This call will initialize all internal and scheduling datastructures.
     * This function must be called prior to any signal function. 
     */
    void initialize() {
      // construct all the required datastructures
      // deinitialize performs the reverse 
      graph.finalize();
      if (rmi.procid() == 0) memory_info::print_usage("Before Engine Initialization");
      logstream(LOG_INFO) 
        << rmi.procid() << ": Initializing..." << std::endl;
        
      // construct scheduler passing in the copy of the options from set_options
      scheduler_ptr = scheduler_factory<message_type>::
                    new_scheduler(graph.num_local_vertices(),
                                  opts_copy);

      // create initial fork arrangement based on the alternate vid mapping
      cmlocks = new distributed_chandy_misra<graph_type>(rmi.dc(), graph,
                                                        boost::bind(&engine_type::lock_ready, this, _1),
                                                        boost::bind(&engine_type::forward_cached_schedule, this, _1));

      cmlocks->compute_initial_fork_arrangement();
      
      // construct the vertex programs
      vstate.resize(graph.num_local_vertices());
      
      // construct the termination consensus object
      consensus = new async_consensus(rmi.dc(), ncpus);

      // if cache is enabled, allocate the cache
      if (use_cache) cache.resize(graph.num_local_vertices());
      
      // finally, the thread local queues
      thrlocal.resize(ncpus);
      if (rmi.procid() == 0) memory_info::print_usage("After Engine Initialization");
      rmi.barrier();
    }
  
  public:
    ~async_consistent_engine() {
      thrlocal.clear();
      delete consensus;
      vstate.clear();
      delete cmlocks;
      delete scheduler_ptr;
    }


    
    
    // documentation inherited from iengine
    size_t num_updates() const { 
      return programs_executed.value; 
    }





    // documentation inherited from iengine
    float elapsed_seconds() const { 
      return timer::approx_time_seconds() - engine_start_time; 
    }


    /**
     * \brief Not meaningful for the asynchronous engine. Returns -1.
     */
    int iteration() const { return -1; }


/**************************************************************************
 *                           Signaling Interface                          *
 **************************************************************************/

  private:

    void internal_post_delta(const vertex_type& vertex,
                             const gather_type& delta) {
      // only post deltas if caching is enabled
      if(use_cache) cache[vertex.local_id()] += delta;
    }

    void internal_clear_gather_cache(const vertex_type& vertex) {
      // only post clear cache if caching is enabled
      if(use_cache) cache[vertex.local_id()].clear();
    }



    /**
     * \internal
     * This is used to receive a message forwarded from another machine
     */
    void rpc_signal(vertex_id_type vid,
                            const message_type& message) {
      if (force_stop) return;
      const lvid_type local_vid = graph.local_vid(vid);
      BEGIN_TRACEPOINT(disteng_scheduler_task_queue);
      bool direct_injection = false;
      // if the vertex is still in the locking state
      // it has not received the message yet.
      // lets directly inject it into the vertex_state message cache
      if (vstate[local_vid].state == LOCKING) {
        vstate[local_vid].lock();
        if (vstate[local_vid].state == LOCKING) {
          vstate[local_vid].current_message += message;
          direct_injection = true;
          joined_messages.inc();
        }
        vstate[local_vid].unlock();
      }
      // if we cannot directly inject into the vertex, then we have no 
      // choice but to put the message into the scheduler
      if (direct_injection == false) {
        scheduler_ptr->schedule(local_vid, message);
      }
      END_TRACEPOINT(disteng_scheduler_task_queue);
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_FROM_REMOTE_EVENT, 1);
      consensus->cancel();
    }

    
    /**
     * \internal
     * This will inject a "placed" task back to be scheduled.
     */
    void signal_local_next(vertex_id_type local_vid) {
     /* This happens when a currently running vertex is activated.
     We cannot return the vertex back into the scheduler since that would
     cause the scheduler to loop indefinitely around a task which
     cannot run. Therefore, we "place" it back in the scheduler, without
     scheduling it, and set a "hasnext" flag in the vertex_state.
     Then when the vertex finishes execution, we must reschedule the vertex
     in the scheduler */
      if (force_stop) return;
      scheduler_ptr->schedule_from_execution_thread(thread::thread_id(),
                                                      local_vid);
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_EVENT, 1);
      consensus->cancel();
    }

    /**
     * \internal
     * \brief Signals a vertex with an optional message
     * 
     * Signals a vertex, and schedules it to be executed in the future.
     * must be called on a vertex accessible by the current machine.
     */
    void internal_signal(const vertex_type& vtx,
                         const message_type& message = message_type()) {
      if (force_stop) return;
      if (started) {
        BEGIN_TRACEPOINT(disteng_scheduler_task_queue);
        scheduler_ptr->schedule_from_execution_thread(thread::thread_id(),
                                                      vtx.local_id(), message);
        consensus->cancel();
        END_TRACEPOINT(disteng_scheduler_task_queue);
      }
      else {
        scheduler_ptr->schedule(vtx.local_id(), message);
      }
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, SCHEDULE_EVENT, 1);
    } // end of schedule


    /**
     * \internal
     * \brief Signals a vertex with an optional message
     * 
     * Signals a global vid, and schedules it to be executed in the future.
     * If current machine does not contain the vertex, it is ignored.
     */
    void internal_signal_gvid(vertex_id_type gvid,
                              const message_type& message = message_type()) {
      if (force_stop) return;
      if (graph.is_master(gvid)) {
        internal_signal(graph.vertex(gvid), message);
      }
    } // end of schedule


    void internal_signal_broadcast(vertex_id_type gvid,
                                   const message_type& message = message_type()) {
      for (size_t i = 0;i < rmi.numprocs(); ++i) {
        rmi.remote_call(i, &async_consistent_engine::internal_signal_gvid,
                        gvid, message);
      }
    } // end of signal_broadcast


    void rpc_internal_stop() { 
      force_stop = true;
      termination_reason = execution_status::FORCED_ABORT;
    }

    /**
     * \brief Force engine to terminate immediately.
     *
     * This function is used to stop the engine execution by forcing
     * immediate termination. 
     */
    void internal_stop() { 
      for (procid_t i = 0;i < rmi.numprocs(); ++i) {
        rmi.remote_call(i, &async_consistent_engine::rpc_internal_stop);
      }
    }

    
  public:



    void signal(vertex_id_type gvid,
                const message_type& message = message_type()) {
      rmi.barrier();
      internal_signal_gvid(gvid, message);
      rmi.barrier();
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
        scheduler_ptr->schedule(lvid, message);    
      }
      rmi.barrier();
    } // end of schedule all


/**************************************************************************
 *                         Computation Processing                         *
 * Internal vertex program scheduling. The functions are arranged roughly *
 * in the order they are used: locking, gathering, applying, scattering   *
 **************************************************************************/

  private:


    /**
     * \internal
     * This function is called after a message to vertex sched_lvid was
     * issued by the scheduler. The message must then be stored in the
     * vertex_state. Calling this function will begin requesting 
     * locks for this vertex on all mirrors.
     */
    void master_broadcast_locking(lvid_type sched_lvid) {
      local_vertex_type lvertex(graph.l_vertex(sched_lvid));

      logstream(LOG_DEBUG) << rmi.procid() << ": Broadcast Gathering: "
                           << lvertex.global_id() << std::endl;
//      ASSERT_I_AM_OWNER(sched_lvid);

      const unsigned char prevkey =
        rmi.dc().set_sequentialization_key(lvertex.global_id() % 254 + 1);
      rmi.remote_call(lvertex.mirrors().begin(), lvertex.mirrors().end(),
                      &engine_type::rpc_begin_locking, lvertex.global_id());
      rmi.dc().set_sequentialization_key(prevkey);
    }


    /**
     * \internal
     * RPC call issued by a master vertex requesting for locks
     * on sched_vid. This switches the vertex to a LOCKING state.
     * Note that, due to communication latencies, it is possible that
     * sched_vid is currently in a SCATTERING phase. If that is the cae,
     * sched_vid is switched to the MIRROR_SCATTERING_AND_NEXT_LOCKING
     * state
     */
    void rpc_begin_locking(vertex_id_type sched_vid) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Mirror Begin Locking: "
                           << sched_vid << std::endl;
      // immediately begin issuing the lock requests
      vertex_id_type sched_lvid = graph.local_vid(sched_vid);
      // set the vertex state
      vstate[sched_lvid].lock();
      if (vstate[sched_lvid].state == NONE) {
        vstate[sched_lvid].state = LOCKING;
        add_internal_task(sched_lvid);
      }
      else if (vstate[sched_lvid].state == MIRROR_SCATTERING) {
        vstate[sched_lvid].state = MIRROR_SCATTERING_AND_NEXT_LOCKING;
      }
      vstate[sched_lvid].unlock();
    }


    /**
     * \internal
     * When all distributed locks are acquired, this function is called
     * from the chandy misra implementation on the master vertex.
     * Here, we perform initialization
     * of the task and switch the vertex to a gathering state
     */
    void lock_ready(lvid_type lvid) {
      // locks for this vertex is ready. perform initialization and
      // and switch the vertex to the gathering state
      logstream(LOG_DEBUG) << "Lock ready on " << "L" << lvid << std::endl;
      vstate[lvid].lock();
      vstate[lvid].state = GATHERING;
      do_init_gather(lvid);
      vstate[lvid].unlock();
      master_broadcast_gathering(lvid, vstate[lvid].vertex_program);
    }


    /**
     * \internal
     * Broadcasts to all machines to begin gathering on sched_lvid.
     * The vertex program is transmitted as well
     */
    void master_broadcast_gathering(lvid_type sched_lvid,
                                    const vertex_program_type& prog) {
      local_vertex_type lvertex(graph.l_vertex(sched_lvid));
      BEGIN_TRACEPOINT(disteng_init_gathering);
      logstream(LOG_DEBUG) << rmi.procid() << ": Broadcast Gathering: "
                           << lvertex.global_id() << std::endl;
//      ASSERT_I_AM_OWNER(sched_lvid);
      // convert to local ID

      const unsigned char prevkey =
        rmi.dc().set_sequentialization_key(lvertex.global_id() % 254 + 1);
      rmi.remote_call(lvertex.mirrors().begin(), lvertex.mirrors().end(),
                      &engine_type::rpc_begin_gathering, lvertex.global_id() , prog);
      rmi.dc().set_sequentialization_key(prevkey);
      END_TRACEPOINT(disteng_init_gathering);
      add_internal_task(sched_lvid);
    }


    /**
     * \internal
     * Called remotely by master_broadcast_gathering.
     * This function stores the vertex program in the vertex
     * state and switches the vertex to a GATHERING state.
     */
    void rpc_begin_gathering(vertex_id_type sched_vid,
                             const vertex_program_type& prog) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Mirror Begin Gathering: "
                           << sched_vid << std::endl;
//      ASSERT_NE(graph.get_vertex_record(sched_vid).owner, rmi.procid());
      // immediately begin issuing the lock requests
      vertex_id_type sched_lvid = graph.local_vid(sched_vid);
      // set the vertex state
      vstate[sched_lvid].lock();
//      ASSERT_EQ(vstate[sched_lvid].state, LOCKING);

      vstate[sched_lvid].state = MIRROR_GATHERING;
      vstate[sched_lvid].vertex_program = prog;
      vstate[sched_lvid].combined_gather.clear();
      vstate[sched_lvid].unlock();
      // lets go
      add_internal_task(sched_lvid);
    }


   /**
     * \internal
     * When a mirror/master finishes its part of a gather, this function
     * is called on the master with the gathered value.
     */
    void rpc_gather_complete(vertex_id_type vid,
                             const conditional_gather_type& uf) {
      logstream(LOG_DEBUG) << rmi.procid() << ": Receiving Gather Complete of "
                           << vid << std::endl;
      lvid_type lvid = graph.local_vid(vid);
      vstate[lvid].lock();
      vstate[lvid].combined_gather += uf;
      decrement_gather_counter(lvid);
      vstate[lvid].unlock();
    }

    
    /**
     * \internal
     * when a machine finishes its part of the gather, it calls
     * decrement_gather_counter which will count the number of machines
     * which have completed the gather. When all machines have responded
     * this function switches the vertex to the APPLYING state.
     * Must be called with locks
     */
    void decrement_gather_counter(const lvid_type lvid) {
      vstate[lvid].apply_count_down--;
      logstream(LOG_DEBUG) << rmi.procid() << ": Partial Gather Complete: " 
                    << graph.global_vid(lvid) << "(" << vstate[lvid].apply_count_down << ")" << std::endl;
      if (vstate[lvid].apply_count_down == 0) {
        logstream(LOG_DEBUG) << rmi.procid() << ": Gather Complete " 
                             << graph.global_vid(lvid) << std::endl;
        vstate[lvid].state = APPLYING;
        add_internal_task(lvid);
      }
    }

    /**
     * \internal
     * Performs the recv_message operation on vertex lvid using
     * the message stored in the vertex_state. locks should be acquired.
     */
    void do_init_gather(lvid_type lvid) {
      context_type context(*this, graph);
      vstate[lvid].vertex_program.recv_message(context,
                                              vertex_type(graph.l_vertex(lvid)),
                                              vstate[lvid].current_message);
      vstate[lvid].current_message = message_type();
      vstate[lvid].combined_gather.clear();
    }
 


    /**
     * \internal
     * Performs the gather operation on vertex lvid. Locks should be acquired.
     */
    void do_gather(lvid_type lvid) { // Do gather
      BEGIN_TRACEPOINT(disteng_evalfac);
      local_vertex_type lvertex(graph.l_vertex(lvid));
      vertex_type vertex(lvertex);

      conditional_gather_type* gather_target = NULL;
      
      if (use_cache) {
        if (cache[lvid].not_empty()) {
          // there is something in the cache. Return that
          vstate[lvid].combined_gather += cache[lvid];
          return;
        }
        else {
          // there is nothing in the cache. Make the cache the
          // gather target, then later write back to the combined gather
          gather_target = &(cache[lvid]);
        }
      }
      else {
        // cache not enabled, gather directly to the combined gather
        gather_target = &(vstate[lvid].combined_gather);
      }
      
      context_type context(*this, graph);

      edge_dir_type gatherdir = vstate[lvid].vertex_program.gather_edges(context, vertex);

      if(gatherdir == graphlab::IN_EDGES ||
        gatherdir == graphlab::ALL_EDGES) {
        foreach(const local_edge_type& edge, lvertex.in_edges()) {
          edge_type e(edge);
          (*gather_target) +=
            vstate[lvid].vertex_program.gather(context, vertex, e);
        }
      }
      if(gatherdir == graphlab::OUT_EDGES ||
        gatherdir == graphlab::ALL_EDGES) {
        foreach(const local_edge_type& edge, lvertex.out_edges()) {
          edge_type e(edge);
          (*gather_target) +=
            vstate[lvid].vertex_program.gather(context, vertex, e);
        }
      }

      if (use_cache && cache[lvid].empty()) {
        // this is the condition where the gather target is the cache
        // now I need to combine it to the combined gather
        vstate[lvid].combined_gather += cache[lvid];
      }
      
      END_TRACEPOINT(disteng_evalfac);
    }

    /**
     * \internal
     * Main function to perform gathers on vertex lvid.
     * This function performs do_gather() on vertex lvid.
     * The resultant gathered value is then sent back to the master
     * for combining.
     */
    void process_gather(lvid_type lvid) {

      const vertex_id_type vid = graph.global_vid(lvid);
      logstream(LOG_DEBUG) << rmi.procid() << ": Gathering on " << vid
                           << std::endl;
      do_gather(lvid);

      const procid_t vowner = graph.l_get_vertex_record(lvid).owner;
      if (vowner == rmi.procid()) {
        decrement_gather_counter(lvid);
      } else {
        vstate[lvid].state = NONE;
        logstream(LOG_DEBUG) << rmi.procid() << ": Send Gather Complete of " << vid
                             << " to " << vowner << std::endl;

        rmi.remote_call(vowner,
                        &engine_type::rpc_gather_complete,
                        graph.global_vid(lvid),
                        vstate[lvid].combined_gather);

        vstate[lvid].combined_gather.clear();
      }
    }

    
    /**
     * \internal
     * Performs the apply operation on vertex lvid using
     * the gathered values stored in the vertex_state.
     * Locks should be acquired.
     */
    void do_apply(lvid_type lvid) { 
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, USER_OP_EVENT, 1);
      BEGIN_TRACEPOINT(disteng_evalfac);
      context_type context(*this, graph);
      
      vertex_type vertex(graph.l_vertex(lvid));
      
      logstream(LOG_DEBUG) << rmi.procid() << ": Apply On " << vertex.id() << std::endl;   
      vstate[lvid].vertex_program.apply(context, 
                                        vertex, 
                                        vstate[lvid].combined_gather.value);
      vstate[lvid].combined_gather.clear();
      END_TRACEPOINT(disteng_evalfac);
    }

    /**
     * \internal
     * Similar to master_broadcast_gathering. Sends the updated vertex program
     * and vertex data to all mirrors. And requests all mirrors to begin
     * scattering.
     */
    void master_broadcast_scattering(lvid_type sched_lvid,
                                     const vertex_program_type& prog,
                                     const vertex_data_type &central_vdata) {
      BEGIN_TRACEPOINT(disteng_init_scattering);
      local_vertex_type lvertex(graph.l_vertex(sched_lvid));
      logstream(LOG_DEBUG) << rmi.procid() << ": Broadcast Scattering: "
                           << lvertex.global_id() << std::endl;
//      ASSERT_I_AM_OWNER(sched_lvid);

      const unsigned char prevkey =
        rmi.dc().set_sequentialization_key(lvertex.global_id() % 254 + 1);
      rmi.remote_call(lvertex.mirrors().begin(), lvertex.mirrors().end(),
                      &engine_type::rpc_begin_scattering,
                      lvertex.global_id(), prog, central_vdata);
      rmi.dc().set_sequentialization_key(prevkey);
      END_TRACEPOINT(disteng_init_scattering);
    }


    /**
     * \internal
     * Called remotely by master_broadcast_scattering.
     * Stores the modified vertex program and vertex data and
     * switches the vertex to a SCATTERING state. 
     */
    void rpc_begin_scattering(vertex_id_type vid,
                              const vertex_program_type& prog,
                              const vertex_data_type &central_vdata) {
      vertex_id_type lvid = graph.local_vid(vid);
      vstate[lvid].lock();
//      ASSERT_I_AM_NOT_OWNER(lvid);
      vstate[lvid].state = MIRROR_SCATTERING;
      graph.get_local_graph().vertex_data(lvid) = central_vdata;
      vstate[lvid].vertex_program = prog;
      vstate[lvid].unlock();
      add_internal_task(lvid);
    }

    
    /**
     * \internal
     * Performs the scatter operation on vertex lvid. locks should be acquired.
     */
    void do_scatter(lvid_type lvid) {
      BEGIN_TRACEPOINT(disteng_evalfac);
      local_vertex_type lvertex(graph.l_vertex(lvid));
      vertex_type vertex(lvertex);

      context_type context(*this, graph);
      
      edge_dir_type scatterdir = vstate[lvid].vertex_program.scatter_edges(context, vertex);
      
      if(scatterdir == graphlab::IN_EDGES || 
         scatterdir == graphlab::ALL_EDGES) {
        foreach(const local_edge_type& edge, lvertex.in_edges()) {
          edge_type e(edge);
          vstate[lvid].vertex_program.scatter(context, vertex, e);
        }
      }
      if(scatterdir == graphlab::OUT_EDGES ||
         scatterdir == graphlab::ALL_EDGES) {
        foreach(const local_edge_type& edge, lvertex.out_edges()) {
          edge_type e(edge);
          vstate[lvid].vertex_program.scatter(context, vertex, e);
        }
      }
      END_TRACEPOINT(disteng_evalfac);
    } // end of do scatter


/**************************************************************************
 *                         Thread Management                              *
 * Functions which manage the thread queues, internal task queues,        *
 * termination, etc                                                       *
 **************************************************************************/
  private:
    /**
     * \internal
     * Key internal dispatch method.
     * Advances the state of a vertex based on its state
     * as decribed in the vertex_state.
     */
    void eval_internal_task(size_t threadid, lvid_type lvid) {
      // if lvid is >= #local vertices, this is an aggregator task
      if (lvid >= vstate.size()) {
        aggregator.tick_asynchronous_compute(threadid,
                                             aggregate_id_to_key[-lvid]);
        return;
      }
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, INTERNAL_TASK_EVENT, 1);
      vstate[lvid].lock();
      switch(vstate[lvid].state) {
      case NONE: 
        ASSERT_MSG(false, "Empty Internal Task");
      case LOCKING: {
          BEGIN_TRACEPOINT(disteng_chandy_misra);
          cmlocks->make_philosopher_hungry_per_replica(lvid);
          END_TRACEPOINT(disteng_chandy_misra);
          break;
      }
      case GATHERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Internal Task: " 
                              << graph.global_vid(lvid) << ": GATHERING(" << vstate[lvid].apply_count_down << ")" << std::endl;

          process_gather(lvid);
          break;
      }
      case MIRROR_GATHERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Internal Task: " 
                              << graph.global_vid(lvid) << ": MIRROR_GATHERING" << std::endl;
          process_gather(lvid);
          break;
        }
      case APPLYING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Internal Task: " 
                              << graph.global_vid(lvid) << ": APPLYING" << std::endl;

          do_apply(lvid);
          vstate[lvid].state = SCATTERING;
          master_broadcast_scattering(lvid,
                                      vstate[lvid].vertex_program,
                                      graph.get_local_graph().vertex_data(lvid));
          // fall through to scattering
        }
      case SCATTERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Scattering: " 
                              << graph.global_vid(lvid) << ": SCATTERING" << std::endl;

          do_scatter(lvid);
          programs_executed.inc();
          BEGIN_TRACEPOINT(disteng_chandy_misra);
          cmlocks->philosopher_stops_eating_per_replica(lvid);
          END_TRACEPOINT(disteng_chandy_misra);

          if (vstate[lvid].hasnext) {
            // stick next back into the scheduler
            signal_local_next(lvid);
            vstate[lvid].hasnext = false;
          } 
          vstate[lvid].state = NONE;
          break;
        }
      case MIRROR_SCATTERING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Scattering: " 
                              << graph.global_vid(lvid) << ": MIRROR_SCATTERING" << std::endl;
          do_scatter(lvid);
          vstate[lvid].state = NONE;
          cmlocks->philosopher_stops_eating_per_replica(lvid);
//          ASSERT_FALSE(vstate[lvid].hasnext);
          break;
        }
      case MIRROR_SCATTERING_AND_NEXT_LOCKING: {
          logstream(LOG_DEBUG) << rmi.procid() << ": Scattering: " 
                              << graph.global_vid(lvid) << ": MIRROR_SCATTERING_AND_NEXT_LOCKING" << std::endl;
          do_scatter(lvid);
          vstate[lvid].state = LOCKING;
//          ASSERT_FALSE(vstate[lvid].hasnext);          
          cmlocks->philosopher_stops_eating_per_replica(lvid);
          cmlocks->make_philosopher_hungry_per_replica(lvid);
          break;
        }
      }
      vstate[lvid].unlock();
    } // end of eval internal task


    /**
     * \internal
     * Picks up a task from the internal queue or the scheduler.
     * \param has_internal_task Set to true, if an internal task is returned
     * \param internal_lvid Contains a queue of internal tasks.
     *                      Set if has_internal_task = true
     * \param has_sched_msg Set to true if a scheduler entry is returned
     * \param sched_lvid A vertex to activate. Set if has_sched_msg = true
     * \param msg A message to send to sched_lvid. Set if has_sched_msg = true
     */
    void get_a_task(size_t threadid,
                    bool& has_internal_task,
                    std::deque<lvid_type>& internal_lvid,
                    bool& has_sched_msg,
                    lvid_type& sched_lvid,
                    message_type &msg) {
      has_internal_task = false;
      has_sched_msg = false;
      if (timer::approx_time_seconds() - engine_start_time > timed_termination) {
        return;
      }
      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      if (thrlocal[threadid].get_task(internal_lvid)) {
        has_internal_task = true;
        END_TRACEPOINT(disteng_internal_task_queue);
        return;
      }
      END_TRACEPOINT(disteng_internal_task_queue);

      if (cmlocks->num_clean_forks() >= max_clean_forks) {
          return;
        }

      sched_status::status_enum stat =
        scheduler_ptr->get_next(threadid, sched_lvid, msg);
      has_sched_msg = stat != sched_status::EMPTY;
    } // end of get a task


    /**
     * \internal
     * Called when get_a_task returns no internal task not a scheduler task.
     * This rechecks the status of the internal task queue and the scheduler
     * inside a consensus critical section.
     */
    bool try_to_quit(size_t threadid,
                     bool& has_internal_task,
                     std::deque<lvid_type>& internal_lvid,
                     bool& has_sched_msg,
                     lvid_type& sched_lvid,
                     message_type &msg) {
      static size_t ctr = 0;
      PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, NO_WORK_EVENT, 1);
      if (timer::approx_time_seconds() - engine_start_time > timed_termination) {
        termination_reason = execution_status::TIMEOUT;
        force_stop = true;
      }
      if (!force_stop &&
          issued_messages.value != programs_executed.value + blocked_issues.value) {
        ++ctr;
        if (ctr % 10 == 0) usleep(1);
        return false;
      }
      logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid << ": " << "Termination Attempt "
                           << programs_executed.value << "/" << issued_messages.value << std::endl;
      has_internal_task = false;
      has_sched_msg = false;
      consensus->begin_done_critical_section(threadid);

      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      if (thrlocal[threadid].get_task(internal_lvid)) {
        logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tCancelled by Internal Task"  << std::endl;
        has_internal_task = true;
        consensus->cancel_critical_section(threadid);
        END_TRACEPOINT(disteng_internal_task_queue);
        return false;
      }
      END_TRACEPOINT(disteng_internal_task_queue);

      sched_status::status_enum stat =
        scheduler_ptr->get_next(threadid, sched_lvid, msg);
      if (stat == sched_status::EMPTY) {
        logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tTermination Double Checked" << std::endl;
        bool ret = consensus->end_done_critical_section(threadid);
        if (ret == false) {
          logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tCancelled" << std::endl;
        }
        return ret;
      } else {
        logstream(LOG_DEBUG) << rmi.procid() << "-" << threadid <<  ": "
                             << "\tCancelled by Scheduler Task" << std::endl;
        consensus->cancel_critical_section(threadid);
        has_sched_msg = true;
        return false;
      }
    } // end of try to quit


    /**
     * \internal
     * Adds a vertex to the internal task queue
     */
    void add_internal_task(lvid_type lvid) {
      if (force_stop) return;
      BEGIN_TRACEPOINT(disteng_internal_task_queue);
      size_t i = lvid % ncpus;
      if (vstate[lvid].state == APPLYING || vstate[lvid].state == SCATTERING ||
          vstate[lvid].state == MIRROR_SCATTERING ||
          vstate[lvid].state == MIRROR_SCATTERING_AND_NEXT_LOCKING) {
        thrlocal[i].add_task_priority(lvid);
      }
      else {
        thrlocal[i].add_task(lvid); 
      }
      consensus->cancel_one(i);
      END_TRACEPOINT(disteng_internal_task_queue);
    }

    void add_internal_aggregation_task(const std::string& key) {
      ASSERT_GT(aggregate_key_to_id.count(key), 0);
      lvid_type id = -aggregate_key_to_id[key];
      // add to all internal task queues
      for (size_t i = 0;i < ncpus; ++i) {
        thrlocal[i].add_task_priority(id);
      }
      consensus->cancel();
    }
    
    
    /**
     * \internal
     * Callback from the Chandy misra implementation.
     * This is called when a vertex acquired all available local forks.
     * This is really an optimization. We extract all available tasks from the
     * on this vertex from the scheduler and forward it to the master.
     */
    void forward_cached_schedule(vertex_id_type lvid) {
      message_type msg;
      const typename graph_type::vertex_record& rec = graph.l_get_vertex_record(lvid);
      if (rec.owner != rmi.procid()) {
        if (scheduler_ptr->get_specific(lvid, msg) == sched_status::NEW_TASK) {
          rmi.remote_call(rec.owner, &engine_type::rpc_signal, rec.gvid, msg);
        }
      }
    }
    /**
     * \internal
     * Called when the scheduler returns a vertex to run.
     * If this function is called with vertex locks acquired, prelocked
     * should be true. Otherwise it should be false.
     */
    template <bool prelocked>
    void eval_sched_task(const lvid_type sched_lvid, 
                         const message_type& msg) {
      BEGIN_TRACEPOINT(disteng_eval_sched_task);
      logstream(LOG_DEBUG) << rmi.procid() << ": Schedule Task: "
                           << graph.global_vid(sched_lvid) << std::endl;
      // If I am not the owner just forward the task to the other
      // scheduler and return
      const typename graph_type::vertex_record& rec = graph.l_get_vertex_record(sched_lvid);
      const procid_t owner = rec.owner;
      bool acquirelock = false;
      if (owner != rmi.procid()) {
        const vertex_id_type vid = rec.gvid;
        rmi.remote_call(owner, &engine_type::rpc_signal, vid, msg);
        return;
      }
//      ASSERT_I_AM_OWNER(sched_lvid);
      // this is in local VIDs
      issued_messages.inc();
      if (prelocked == false) {
        vstate[sched_lvid].lock();
      }
      if (vstate[sched_lvid].state == NONE) {
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, UPDATE_EVENT, 1);
        PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, WORK_ISSUED_EVENT, rec.num_in_edges + rec.num_out_edges);

        // we start gather right here.
        // set up the state
        vstate[sched_lvid].state = LOCKING;
        vstate[sched_lvid].hasnext = false;
        vstate[sched_lvid].current_message = msg;
        vstate[sched_lvid].apply_count_down = graph.l_vertex(sched_lvid).num_mirrors() + 1;
        acquirelock = true;
        // we are going to broadcast after unlock
      } else if (vstate[sched_lvid].state == LOCKING) {
         blocked_issues.inc();
         vstate[sched_lvid].current_message += msg;
         joined_messages.inc();
      } else {
        blocked_issues.inc();
        if (vstate[sched_lvid].hasnext) {
          scheduler_ptr->place(sched_lvid, msg);
          joined_messages.inc();
        } else {
          vstate[sched_lvid].hasnext = true;
          scheduler_ptr->place(sched_lvid, msg);
        }
      }
      if (prelocked == false) vstate[sched_lvid].unlock();
      END_TRACEPOINT(disteng_eval_sched_task);
      if (acquirelock) {
        BEGIN_TRACEPOINT(disteng_chandy_misra);
        cmlocks->make_philosopher_hungry_per_replica(sched_lvid);
        END_TRACEPOINT(disteng_chandy_misra);
        master_broadcast_locking(sched_lvid);
      }
    }
    
    
    
    /**
     * \internal
     * Per thread main loop
     */
    void thread_start(size_t threadid) {
      bool has_internal_task = false;
      bool has_sched_msg = false;
      std::deque<vertex_id_type> internal_lvid;
      vertex_id_type sched_lvid;
      message_type msg;
//      size_t ctr = 0;
      float last_aggregator_check = timer::approx_time_seconds();
      while(1) {
        if (timer::approx_time_seconds() != last_aggregator_check) {
          last_aggregator_check = timer::approx_time_seconds();
          std::string key = aggregator.tick_asynchronous();
          if (key != "") add_internal_aggregation_task(key);
        }
/*        ++ctr;
        if (max_clean_forks != (size_t)(-1) && ctr % 10000 == 0) {
          std::cout << cmlocks->num_clean_forks() << "/" << max_clean_forks << "\n";
        }*/
       get_a_task(threadid, 
                   has_internal_task, internal_lvid,
                   has_sched_msg, sched_lvid, msg);
        // if we managed to get a task..
        if (has_internal_task) {
          while(!internal_lvid.empty()) {
            eval_internal_task(threadid, internal_lvid.front());
            internal_lvid.pop_front();
          }
        } else if (has_sched_msg) {
          eval_sched_task<false>(sched_lvid, msg);
        }
        /*
         * We failed to obtain a task, try to quit
         */
        else if (!try_to_quit(threadid,
                              has_internal_task, internal_lvid,
                              has_sched_msg, sched_lvid, msg)) {
          if (has_internal_task) {
            while(!internal_lvid.empty()) {
              eval_internal_task(threadid, internal_lvid.front());
              internal_lvid.pop_front();
            }
          } else if (has_sched_msg) {
            eval_sched_task<false>(sched_lvid, msg);
          }
        } else { break; }
      }
    } // end of thread start


/**************************************************************************
 *                         For init vertex program                        *
 * For convenience, the init() of vertex programs are managed by a        *
 * separate set of threads which are only activated at the start.         *
 * This is copied from the synchronous engine.                            *
***************************************************************************/
  private:
    /// \internal pair of vertex id and vertex data
    typedef std::pair<vertex_id_type, vertex_data_type> vid_vdata_pair_type;
    /// \internal Exchange type used to swap vertex data on init program
    typedef buffered_exchange<vid_vdata_pair_type> vdata_exchange_type;
    /// \internal Exchange structure used to swap vertex data on init program
    vdata_exchange_type vdata_exchange;

    /// \internal Synchronizes init vertex program threads
    barrier thread_barrier;

    /**
    * \internal
    * synchronizes the local vertex lvid
    */
    void sync_vertex_data(lvid_type lvid) {
      ASSERT_TRUE(graph.l_is_master(lvid));
      const vertex_id_type vid = graph.global_vid(lvid);
      local_vertex_type vertex = graph.l_vertex(lvid);
      foreach(const procid_t& mirror, vertex.mirrors()) {
        vdata_exchange.send(mirror, std::make_pair(vid, vertex.data()));
      }
    } // end of sync_vertex_data

    /**
    * \internal
    * Receives data from the enchange and saves it
    */
    void recv_vertex_data() {
      procid_t procid(-1);
      typename vdata_exchange_type::buffer_type buffer;
      while(vdata_exchange.recv(procid, buffer)) {
        foreach(const vid_vdata_pair_type& pair, buffer) {
          const lvid_type lvid = graph.local_vid(pair.first);
          ASSERT_FALSE(graph.l_is_master(lvid));
          graph.l_vertex(lvid).data() = pair.second;
        }
      }
    } // end of recv vertex data

    /**
    * \internal
    * Called simultaneously by all machines to initialize all vertex programs
    */
    void initialize_vertex_programs(size_t thread_id) {
      // For now we are using the engine as the context interface
      context_type context(*this, graph);
      for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices();
          lvid += ncpus) {
        if(graph.l_is_master(lvid)) {
          vertex_type vertex = local_vertex_type(graph.l_vertex(lvid));
          vstate[lvid].vertex_program.init(context, vertex);
          sync_vertex_data(lvid);
        }
        recv_vertex_data();
      }
      // Flush the buffer and finish receiving any remaining vertex
      // programs.
      thread_barrier.wait();
      if(thread_id == 0) { vdata_exchange.flush(); }
      thread_barrier.wait();
      recv_vertex_data();
    } // end of initialize_vertex_programs



/**************************************************************************
 *                         Main engine start()                            *
 **************************************************************************/

  public:
   
    /**
      * \brief Start the engine execution.
      *
      * This function starts the engine and does not
      * return until the scheduler has no tasks remaining.
      * If perform_init is set to true then the init
      * function is called once at the beginning of start().
      * 
      * \param perform_init_vertex_program If true, runs init on each
      * vertex program before any message processing happens.
      * Defaults to true.
      *
      * \return the reason for termination
      */
    execution_status::status_enum start(bool perform_init_vertex_program = true) {
      logstream(LOG_INFO) << "Spawning " << ncpus << " threads" << std::endl;
      ASSERT_TRUE(scheduler_ptr != NULL);
      consensus->reset();
      // start the scheduler
      scheduler_ptr->start();
      
      // start the aggregator
      aggregator.start(ncpus);
      aggregator.aggregate_all_periodic();
      generate_aggregate_keys();
      // now since I am overloading the internal task queue IDs to also
      // hold aggregator keys, this constraints the number of vertices
      // I can handle. Check for overflow
      {
        lvid_type lv = (lvid_type)graph.num_local_vertices();
        ASSERT_MSG(lv + aggregate_id_to_key.size() >= lv,
                   "Internal Queue IDs numeric overflow");
      }
      started = true;

      rmi.barrier();

      size_t allocatedmem = memory_info::allocated_bytes();
      rmi.all_reduce(allocatedmem);

      engine_start_time = timer::approx_time_seconds();
      force_stop = false;

      rmi.dc().flush_counters();

      termination_reason = execution_status::RUNNING;

      if (perform_init_vertex_program) {
        logstream(LOG_INFO) << "Initialize Vertex Programs: " << allocatedmem << std::endl;
        for (size_t i = 0; i < ncpus; ++i) {
          thrgroup.launch(boost::bind(&engine_type::initialize_vertex_programs, this, i));
        }
        thrgroup.join();
      }

      if (rmi.procid() == 0) {
        logstream(LOG_INFO) << "Total Allocated Bytes: " << allocatedmem << std::endl;
        PERMANENT_IMMEDIATE_DIST_EVENT(eventlog, ENGINE_START_EVENT);
      }
      for (size_t i = 0; i < ncpus; ++i) {
        thrgroup.launch(boost::bind(&engine_type::thread_start, this, i));
      }
      thrgroup.join();
      aggregator.stop();
      // if termination reason was not changed, then it must be depletion
      if (termination_reason == execution_status::RUNNING) {
        termination_reason = execution_status::TASK_DEPLETION;
      }

      rmi.dc().flush_counters();
      if (rmi.procid() == 0) {
        PERMANENT_IMMEDIATE_DIST_EVENT(eventlog, ENGINE_STOP_EVENT);
      }
      size_t ctasks = programs_executed.value;
      rmi.all_reduce(ctasks);
      programs_executed.value = ctasks;

      ctasks = issued_messages.value;
      rmi.all_reduce(ctasks);
      issued_messages.value = ctasks;

      ctasks = blocked_issues.value;
      rmi.all_reduce(ctasks);
      blocked_issues.value = ctasks;

      ctasks = joined_messages.value;
      ctasks += scheduler_ptr->num_joins();
      rmi.all_reduce(ctasks);
      joined_messages.value = ctasks;

      if (rmi.procid() == 0) {
        std::cout << "Completed Tasks: " << programs_executed.value << std::endl;
        std::cout << "Issued Tasks: " << issued_messages.value << std::endl;
        std::cout << "Blocked Issues: " << blocked_issues.value << std::endl;
        std::cout << "------------------" << std::endl;
        std::cout << "Joined Tasks: " << joined_messages.value << std::endl;
      }

      /*for (size_t i = 0;i < vstate.size(); ++i) {
          if(vstate[i].state != NONE) {
            std::cout << "Vertex: " << i << ": " << vstate[i].state << " " << (int)(cmlocks->philosopherset[i].state) << " " << cmlocks->philosopherset[i].num_edges << " " << cmlocks->philosopherset[i].forks_acquired << "\n";

            foreach(typename local_graph_type::edge_type edge, cmlocks->graph.in_edges(i)) {
              std::cout << (int)(cmlocks->forkset[cmlocks->graph.edge_id(edge)]) << " ";
            }
            std::cout << "\n";
            foreach(typename local_graph_type::edge_type edge, cmlocks->graph.out_edges(i)) {
              std::cout << (int)(cmlocks->forkset[cmlocks->graph.edge_id(edge)]) << " ";
            }
            std::cout << "\n";
            getchar();
          }
        }*/
      started = false;
      return termination_reason; 
    } // end of start

/************************************************************
 *                  Aggregators                             *
 ************************************************************/
  private:
    std::map<std::string, lvid_type> aggregate_key_to_id;
    std::vector<std::string> aggregate_id_to_key;

    /**
     * \internal
     * Assigns each periodic aggregation key a unique ID so we can use it in
     * internal task scheduling.
     */
    void generate_aggregate_keys() {
      aggregate_key_to_id.clear();
      aggregate_id_to_key.clear();
      // no ID 0. since we are using negatives for scheduling
      aggregate_id_to_key.push_back("");  
      std::set<std::string> keys = aggregator.get_all_periodic_keys();
      
      foreach(std::string key, keys) {
        aggregate_id_to_key.push_back(key);
        aggregate_key_to_id[key] = (lvid_type)(aggregate_id_to_key.size() - 1);
        
      }
    }
  public:
    // // Exposed aggregator functionality 
    // /**
    //  * \copydoc distributed_aggregator::add_vertex_aggregator()
    //  */
    // template <typename ReductionType>
    // bool add_vertex_aggregator(const std::string& key,
    //   boost::function<ReductionType(icontext_type&,
    //                                 vertex_type&)> map_function,
    //   boost::function<void(icontext_type&,
    //                        const ReductionType&)> finalize_function) {
    //   rmi.barrier();
    //   return aggregator.add_vertex_aggregator<ReductionType>(key,
    //                                                     map_function,
    //                                                     finalize_function);
    // }
    
    // /**
    //  * \copydoc distributed_aggregator::add_edge_aggregator()
    //  */
    // template <typename ReductionType>
    // bool add_edge_aggregator(const std::string& key,
    //   boost::function<ReductionType(icontext_type&,
    //                                 edge_type&)> map_function,
    //   boost::function<void(icontext_type&,
    //                        const ReductionType&)> finalize_function) {
    //   rmi.barrier();
    //   return aggregator.add_edge_aggregator<ReductionType>(key,
    //                                                        map_function,
    //                                                        finalize_function);
    // }
    
    // /**
    //  * \copydoc distributed_aggregator::aggregate_now()
    //  */
    // bool aggregate_now(const std::string& key) {
    //   rmi.barrier();
    //   return aggregator.aggregate_now(key);
    // }
    
    // /**
    //  * \copydoc distributed_aggregator::aggregate_periodic()
    //  */
    // bool aggregate_periodic(const std::string& key, float seconds) {
    //   rmi.barrier();
    //   return aggregator.aggregate_periodic(key, seconds);
    // }
  private:
    aggregator_type* get_aggregator() { return &aggregator; } 

  }; // end of class
} // namespace

#include <graphlab/macros_undef.hpp>

#endif // GRAPHLAB_DISTRIBUTED_ENGINE_HPP

