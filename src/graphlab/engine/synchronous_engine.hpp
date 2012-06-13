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



#ifndef GRAPHLAB_SYNCHRONOUS_ENGINE_HPP
#define GRAPHLAB_SYNCHRONOUS_ENGINE_HPP

#include <deque>
#include <boost/bind.hpp>

#include <graphlab/engine/iengine.hpp>

#include <graphlab/vertex_program/ivertex_program.hpp>
#include <graphlab/vertex_program/icontext.hpp>
#include <graphlab/vertex_program/context.hpp>

#include <graphlab/engine/execution_status.hpp>
#include <graphlab/options/graphlab_options.hpp>




#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic_add_vector.hpp>
#include <graphlab/util/tracepoint.hpp>
#include <graphlab/util/memory_info.hpp>

#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>






#include <graphlab/macros_def.hpp>

namespace graphlab {
  
  
  /**
   * \ingroup engines
   * 
   * \brief The synchronous engine executes all active vertex program
   * synchronously in a sequence of super-step (iterations) in both the
   * shared and distributed memory settings.
   * 
   * \tparam VertexProgram The user defined vertex program which
   * should implement the \ref graphlab::ivertex_program interface.   
   *
   *  
   * ### Execution Semantics
   * 
   * On start() the \ref graphlab::ivertex_program::init function is invoked
   * on all vertex programs in parallel to initialize the vertex program,
   * vertex data, and possibly signal vertices.
   * The engine then proceeds to execute a sequence of
   * super-steps (iterations) each of which is further decomposed into a 
   * sequence of minor-steps which are also executed synchronously:
   * \li Receive all incoming messages (signals) by invoking the 
   * \ref graphlab::ivertex_program::recv_message function on all
   * vertex-programs that have incoming messages.  If a
   * vertex-program does not have any incoming messages then it is
   * not active during this super-step.  
   * \li Execute all gathers for active vertex programs by invoking
   * the user defined \ref graphlab::ivertex_program::gather function
   * on the edge direction returned by the 
   * \ref graphlab::ivertex_program::gather_edges function.  The gather
   * functions can modify edge data but cannot modify the vertex
   * program or vertex data and therefore can be executed on multiple
   * edges in parallel.  The gather type is used to accumulate (sum)
   * the result of the gather function calls.
   * \li Execute all apply functions for active vertex-programs by
   * invoking the user defined \ref graphlab::ivertex_program::apply
   * function passing the sum of the gather functions.  If \ref
   * graphlab::ivertex_program::gather_edges returns no edges then
   * the default gather value is passed to apply.  The apply function
   * can modify the vertex program and vertex data.
   * \li Execute all scatters for active vertex programs by invoking
   * the user defined \ref graphlab::ivertex_program::scatter function
   * on the edge direction returned by the 
   * \ref graphlab::ivertex_program::scatter_edges function.  The scatter
   * functions can modify edge data but cannot modify the vertex
   * program or vertex data and therefore can be executed on multiple
   * edges in parallel.  
   * 
   * ### Construction
   *
   * The synchronous engine is constructed by passing in a 
   * \ref graphlab::distributed_control object which manages coordination 
   * between engine threads and a \ref graphlab::distributed_graph object 
   * which is the graph on which the engine should be run.  The graph should
   * already be populated and cannot change after the engine is constructed. 
   * In the distributed setting all program instances (running on each machine)
   * should construct an instance of the engine at the same time.
   * 
   * Computation is initiated by signaling vertices using either 
   * \ref graphlab::synchronous_engine::signal or 
   * \ref graphlab::synchronous_engine::signal_all.  In either case all
   * machines should invoke signal or signal all at the same time.  Finally,
   * computation is initiated by calling the 
   * \ref graphlab::synchronous_engine::start function. 
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
   *   graph.load_structure(graph_dir, "tsv");
   *   graph.finalize();
   *   std::cout << "#vertices: " << graph.num_vertices() 
   *             << " #edges:" << graph.num_edges() << std::endl;
   *   graphlab::synchronous_engine<pagerank_vprog> engine(dc, graph, clopts);
   *   engine.signal_all();
   *   engine.start();
   *   std::cout << "Runtime: " << engine.elapsed_time();
   *   graphlab::mpi_tools::finalize();
   * }
   * \endcode
   *
   * \see graphlab::async_consistent_engine
   */
  template<typename VertexProgram>
  class synchronous_engine : 
    public iengine<VertexProgram> {

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
    typedef typename VertexProgram::graph_type  graph_type;

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
    typedef icontext<graph_type, gather_type, message_type> icontext_type;
           
  private:

    /**
     * \brief Local vertex type used by the engine for fast indexing
     */
    typedef typename graph_type::local_vertex_type    local_vertex_type;
    
    /**
     * \brief Local edge type used by the engine for fast indexing
     */
    typedef typename graph_type::local_edge_type      local_edge_type;
    
    /**
     * \brief Local vertex id type used by the engine for fast indexing
     */
    typedef typename graph_type::lvid_type            lvid_type;

    /**
     * \brief The actual instance of the context type used by this engine.
     */
    typedef context<synchronous_engine> context_type;
    friend class context<synchronous_engine>;


    /**
     * \brief The type of the distributed aggregator inherited from iengine
     */
    typedef typename iengine<vertex_program_type>::aggregator_type aggregator_type;

    /**
     * \brief The object used to communicate with remote copies of the
     * synchronous engine.
     */
    dc_dist_object< synchronous_engine<VertexProgram> > rmi;

    /**
     * \brief A reference to the distributed graph on which this
     * synchronous engine is running.
     */
    graph_type& graph;

    /**
     * \brief The local worker threads used by this engine
     */
    thread_pool threads;

    /**
     * \brief A thread barrier that is used to control the threads in the
     * thread pool.
     */
    graphlab::barrier thread_barrier;

    /**
     * \brief The maximum number of super-steps (iterations) to run
     * before terminating.  If the max iterations is reached the
     * engine will terminate if their are no messages remaining.
     */
    size_t max_iterations;

    /**
     * \brief A counter that tracks the current iteration number since
     * start was last invoked.
     */
    size_t iteration_counter;

    /**
     * \brief The time in seconds at which the engine started.
     */
    float start_time;

    /**
     * \brief The timeout time in seconds
     */
    float timeout;
    
    /**
     * \brief Used to stop the engine prematurely
     */
    bool force_abort;

    /**
     * \brief The vertex locks protect access to vertex specific
     * data-structures including 
     * \ref graphlab::synchronous_engine::gather_accum
     * and \ref graphlab::synchronous_engine::messages.
     */
    std::vector<mutex> vlocks;

    /**
     * \brief The vertex programs associated with each vertex on this
     * machine.
     */
    std::vector<vertex_program_type> vertex_programs;

    /**
     * \brief Vector of messages associated with each vertex.
     */
    std::vector<message_type> messages;

    /**
     * \brief Bit indicating whether a message is present for each vertex.
     */
    dense_bitset has_message;
 

    /**
     * \brief Gather accumulator used for each master vertex to merge
     * the result of all the machine specific accumulators (or
     * caches).
     *
     * The gather accumulator can be accessed by multiple threads at
     * once and therefore must be guarded by a vertex locks in 
     * \ref graphlab::synchronous_engine::vlocks
     */
    std::vector<gather_type>  gather_accum;
    
    /**
     * \brief Bit indicating if the gather has accumulator contains any
     * values.  
     *
     * While dense bitsets are thread safe the value of this bit must
     * change concurrently with the 
     * \ref graphlab::synchronous_engine::gather_accum and therefore is
     * set while holding the lock in
     * \ref graphlab::synchronous_engine::vlocks.
     */
    dense_bitset has_gather_accum;


    /**
     * \brief This optional vector contains caches of previous gather
     * contributions for each machine.  
     *
     * Caching is done locally and therefore a high-degree vertex may
     * have multiple caches (one per machine).
     */
    std::vector<gather_type>  gather_cache;

    /**
     * \brief A bit indicating if the local gather for that vertex is
     * available.
     */
    dense_bitset has_cache;   

    /**
     * \brief A bit (for master vertices) indicating if that vertex is active
     * (received a message on this iteration).
     */
    dense_bitset active_superstep;

    /**
     * \brief  The number of local vertices (masters) that are active on this
     * iteration.
     */
    atomic<size_t> num_active_vertices;

    /**
     * \brief A bit indicating (for all vertices) whether to
     * participate in the current minor-step (gather or scatter).
     */
    dense_bitset active_minorstep;      

    /**
     * \brief A counter measuring the number of applys that have been completed
     */
    atomic<size_t> completed_applys;

    /**
     * \brief The pair type used to synchronize vertex programs across machines.
     */
    typedef std::pair<vertex_id_type, vertex_program_type> vid_prog_pair_type;

    /**
     * \brief The type of the exchange used to synchronize vertex programs
     */
    typedef buffered_exchange<vid_prog_pair_type> vprog_exchange_type;
   
    /**
     * \brief The distributed exchange used to synchronize changes to
     * vertex programs.
     */
    vprog_exchange_type vprog_exchange;

    /**
     * \brief The pair type used to synchronize vertex across across machines.
     */
    typedef std::pair<vertex_id_type, vertex_data_type> vid_vdata_pair_type;

    /**
     * \brief The type of the exchange used to synchronize vertex data
     */
    typedef buffered_exchange<vid_vdata_pair_type> vdata_exchange_type;

    /**
     * \brief The distributed exchange used to synchronize changes to
     * vertex programs.
     */
    vdata_exchange_type vdata_exchange;

    /**
     * \brief The pair type used to synchronize the results of the gather phase
     */
    typedef std::pair<vertex_id_type, gather_type> vid_gather_pair_type;

    /**
     * \brief The type of the exchange used to synchronize gather
     * accumulators
     */
    typedef buffered_exchange<vid_gather_pair_type> gather_exchange_type;
   
    /**
     * \brief The distributed exchange used to synchronize gather
     * accumulators.
     */
    gather_exchange_type gather_exchange;

    /**
     * \brief The pair type used to synchronize messages
     */
    typedef std::pair<vertex_id_type, message_type> vid_message_pair_type;
   
    /**
     * \brief The type of the exchange used to synchronize messages
     */
    typedef buffered_exchange<vid_message_pair_type> message_exchange_type;

    /**
     * \brief The distributed exchange used to synchronize messages
     */
    message_exchange_type message_exchange;


    /**
     * \brief The distributed aggregator used to manage background
     * aggregation.
     */
    aggregator_type aggregator;

  public:

    /**
     * \brief Construct a synchronous engine for a given graph and options.
     *
     * The synchronous engine should be constructed after the graph
     * has been loaded (e.g., \ref graphlab::distributed_graph::load)
     * and the graphlab options have been set 
     * (e.g., \ref graphlab::command_line_options).
     *
     * In the distributed engine the synchronous engine must be called
     * on all machines at the same time (in the same order) passing
     * the \ref graphlab::distributed_control object.  Upon
     * construction the synchronous engine allocates several
     * data-structures to store messages, gather accumulants, and
     * vertex programs and therefore may require considerable memory.
     *
     * The number of threads to create are read from 
     * \ref graphlab_options::get_ncpus "opts.get_ncpus()". 
     * Valid engine options (graphlab_options::get_engine_args()):
     * \arg \c use_cache If set to true, partial gathers are cached.
     * See \ref gather_caching to understand the behavior of the
     * gather caching model and how it may be used to accelerate program
     * performance.
     * \argc \c iterations Limit the number of iterations the engine 
     * may run.
     * 
     * @param [in] dc Distributed controller to associate with
     * @param [in,out] graph A reference to the graph object that this
     * engine will modify. The graph must be fully constructed and 
     * finalized.
     * @param [in] opts A graphlab::graphlab_options object specifying engine
     *                  parameters.  This is typically constructed using
     *                  \ref graphlab::command_line_options.
     */
    synchronous_engine(distributed_control& dc, graph_type& graph,
                       const graphlab_options& opts);


    /**
     * \brief Start execution of the synchronous engine.
     * 
     * The start function begins computation and does not return until
     * there are no remaining messages or until max_iterations has
     * been reached. 
     * 
     * The start() function modifies the data graph through the vertex
     * programs and so upon return the data graph should contain the
     * result of the computation.
     *
     * @return The reason for termination
     */
    execution_status::status_enum start();

    // documentation inherited from iengine
    size_t num_updates() const;

    // documentation inherited from iengine
    void signal(vertex_id_type vid,
                const message_type& message = message_type());
    
    // documentation inherited from iengine
    void signal_all(const message_type& message = message_type(),
                    const std::string& order = "sequential");

    // documentation inherited from iengine
    float elapsed_seconds() const;
    
    /**
     * \brief Get the current iteration number since start was last
     * invoked.
     *
     *  \return the current iteration
     */
    int iteration() const; 


    /**
     * \brief Compute the total memory used by the entire distributed
     * system.
     * 
     * @return The total memory used in bytes.
     */
    size_t total_memory_usage() const;

    /**
     * \brief Get a pointer to the distributed aggregator object. 
     *
     * This is currently used by the \ref graphlab::iengine interface to 
     * implement the calls to aggregation.
     *
     * @return a pointer to the local aggregator.
     */
    aggregator_type* get_aggregator();

  private:

    /**
     * \brief This internal stop function is called by the \ref graphlab::context to
     * terminate execution of the engine.
     */
    void internal_stop();    

    /**
     * \brief This function is called remote by the rpc to force the
     * engine to stop.
     */
    void rpc_stop();

    /**
     * \brief Signal a vertex.
     *
     * This function is called by the \ref graphlab::context.
     *
     * @param [in] vertex the vertex to signal
     * @param [in] message the message to send to that vertex.
     */
    void internal_signal(const vertex_type& vertex,
                         const message_type& message = message_type()); 

    /**
     * \brief Called by the context to signal an arbitrary vertex.
     * This must be done by finding the owner of that vertex. 
     *
     * @param [in] gvid the global vertex id of the vertex to signal
     * @param [in] message the message to send to that vertex.
     */
    void internal_signal_broadcast(vertex_id_type gvid,
                                   const message_type& message = message_type());
    
    /**
     * \brief This function tests if this machine is the master of
     * gvid and signals if successful.
     */
    void internal_signal_rpc(vertex_id_type gvid,
                              const message_type& message = message_type());


    /** 
     * \brief Post a to a previous gather for a give vertex.
     *
     * This function is called by the \ref graphlab::context.
     *
     * @param [in] vertex The vertex to which to post a change in the sum
     * @param [in] delta The change in that sum
     */
    void internal_post_delta(const vertex_type& vertex,
                             const gather_type& delta);

    /**
     * \brief Clear the cached gather for a vertex if one is
     * available.
     *
     * This function is called by the \ref graphlab::context.
     *
     * @param [in] vertex the vertex for which to clear the cache
     */
    void internal_clear_gather_cache(const vertex_type& vertex);


    // Program Steps ==========================================================
   
    /**
     * \brief Executes ncpus copies of a member function each with a
     * unique consecutive id (thread id). 
     *
     * This function is used by the main loop to execute each of the
     * stages in parallel.
     *
     * The member function must have the type:
     *
     * \code
     * void synchronous_engine::member_fun(size_t threadid);
     * \endcode
     *
     * This function runs an rmi barrier after termination
     *
     * @tparam the type of the member function.  
     * @param [in] member_fun the function to call.
     */
    template<typename MemberFunction>       
    void run_synchronous(MemberFunction member_fun) {
      // launch the initialization threads
      for(size_t i = 0; i < threads.size(); ++i) 
        threads.launch(boost::bind(member_fun, this, i));
      // Wait for all threads to finish
      threads.join();
      rmi.barrier();
    } // end of run_synchronous

    // /** 
    //  * \brief Initialize all vertex programs by invoking 
    //  * \ref graphlab::ivertex_program::init on all vertices.
    //  *
    //  * @param thread_id the thread to run this as which determines
    //  * which vertices to process.
    //  */
    // void initialize_vertex_programs(size_t thread_id);

    /** 
     * \brief Synchronize all message data.
     *
     * @param thread_id the thread to run this as which determines
     * which vertices to process.
     */
    void exchange_messages(size_t thread_id);


    /** 
     * \brief Invoke the \ref graphlab::ivertex_program::recv_message function 
     * on all vertex programs that have inbound messages.
     *
     * @param thread_id the thread to run this as which determines
     * which vertices to process.
     */
    void receive_messages(size_t thread_id);


    /** 
     * \brief Execute the \ref graphlab::ivertex_program::gather function on all 
     * vertices that received messages for the edges specified by the 
     * \ref graphlab::ivertex_program::gather_edges.
     *
     * @param thread_id the thread to run this as which determines
     * which vertices to process.
     */
    void execute_gathers(size_t thread_id);


    /** 
     * \brief Execute the \ref graphlab::ivertex_program::apply function on all
     * all vertices that received messages in this super-step (active).
     *
     * @param thread_id the thread to run this as which determines
     * which vertices to process.
     */
    void execute_applys(size_t thread_id);

    /** 
     * \brief Execute the \ref graphlab::ivertex_program::scatter function on all 
     * vertices that received messages for the edges specified by the 
     * \ref graphlab::ivertex_program::scatter_edges.
     *
     * @param thread_id the thread to run this as which determines
     * which vertices to process.
     */
    void execute_scatters(size_t thread_id);

    // Data Synchronization ===================================================
    /**
     * \brief Send the vertex program for the local vertex id to all
     * of its mirrors. 
     *
     * @param [in] lvid the vertex to sync.  This muster must be the
     * master of that vertex.
     */
    void sync_vertex_program(lvid_type lvid);

    /**
     * \brief Receive all incoming vertex programs and update the
     * local mirrors.
     *
     * This function returns when there are no more incoming vertex
     * programs and should be called after a flush of the vertex
     * program exchange.
     */
    void recv_vertex_programs();

    /**
     * \brief Send the vertex data for the local vertex id to all of
     * its mirrors.
     *
     * @param [in] lvid the vertex to sync.  This machine must be the master
     * of that vertex.
     */
    void sync_vertex_data(lvid_type lvid);
    
    /**
     * \brief Receive all incoming vertex data and update the local
     * mirrors.
     *
     * This function returns when there are no more incoming vertex
     * data and should be called after a flush of the vertex data
     * exchange.
     */
    void recv_vertex_data();

    /**
     * \brief Send the gather value for the vertex id to its master.
     *
     * @param [in] lvid the vertex to send the gather value to
     * @param [in] accum the locally computed gather value. 
     */
    void sync_gather(lvid_type lvid, const gather_type& accum);


    /**
     * \brief Receive the gather values from the buffered exchange.
     *
     * This function returns when there is nothing left in the
     * buffered exchange and should be called after the buffered
     * exchange has been flushed
     */
    void recv_gathers();

    /**
     * \brief Send the accumulated message for the local vertex to its
     * master.
     *
     * @param [in] lvid the vertex to send 
     */
    void sync_message(lvid_type lvid);

    /**
     * \brief Receive the messages from the buffered exchange.
     *
     * This function returns when there is nothing left in the
     * buffered exchange and should be called after the buffered
     * exchange has been flushed
     */
    void recv_messages();

  }; // end of class synchronous engine







  /**
   * Constructs an synchronous distributed engine.
   * The number of threads to create are read from
   * opts::get_ncpus().
   *
   * Valid engine options (graphlab_options::get_engine_args()):
   * \arg \c max_iterations Sets the maximum number of iterations the
   * engine will run for. 
   * \arg \c use_cache If set to true, partial gathers are cached.
   * See \ref gather_caching to understand the behavior of the
   * gather caching model and how it may be used to accelerate program
   * performance.
   *
   * \param dc Distributed controller to associate with
   * \param graph The graph to schedule over. The graph must be fully
   *              constructed and finalized.
   * \param opts A graphlab_options object containing options and parameters
   *             for the engine.
   */
  template<typename VertexProgram>
  synchronous_engine<VertexProgram>::
  synchronous_engine(distributed_control &dc, 
                     graph_type& graph,
                     const graphlab_options& opts) :
    rmi(dc, this), graph(graph), 
    threads(opts.get_ncpus()), 
    thread_barrier(opts.get_ncpus()),
    max_iterations(-1), iteration_counter(0),
    timeout(0),
    vprog_exchange(dc), vdata_exchange(dc), 
    gather_exchange(dc), message_exchange(dc),
    aggregator(dc, graph, new context_type(*this, graph)) {
    // Process any additional options
    std::vector<std::string> keys = opts.get_engine_args().get_option_keys();
    bool use_cache = false;
    foreach(std::string opt, keys) {
      if (opt == "max_iterations") {
        opts.get_engine_args().get_option("max_iterations", max_iterations);
      } else if (opt == "timeout") {
        opts.get_engine_args().get_option("timeout", timeout);
      } else if (opt == "use_cache") {
        opts.get_engine_args().get_option("use_cache", use_cache);
      } else {
        logstream(LOG_FATAL) << "Unexpected Engine Option: " << opt << std::endl;
      }
    }
    // Finalize the graph
    graph.finalize();
    memory_info::log_usage("Before Engine Initialization");
    // Allocate vertex locks and vertex programs
    vlocks.resize(graph.num_local_vertices());
    vertex_programs.resize(graph.num_local_vertices());
    // Allocate messages and message bitset
    messages.resize(graph.num_local_vertices(), message_type());
    has_message.resize(graph.num_local_vertices()); 
    has_message.clear();
    // Allocate gather accumulators and accumulator bitset
    gather_accum.resize(graph.num_local_vertices(), gather_type());
    has_gather_accum.resize(graph.num_local_vertices());
    has_gather_accum.clear();
    // If caching is used then allocate cache data-structures
    if (use_cache) {
      gather_cache.resize(graph.num_local_vertices(), gather_type());
      has_cache.resize(graph.num_local_vertices());
      has_cache.clear();
    }
    // Allocate bitset to track active vertices on each bitset.
    active_superstep.resize(graph.num_local_vertices());
    active_superstep.clear();
    active_minorstep.resize(graph.num_local_vertices());
    active_minorstep.clear();
    // Print memory usage after initialization
    memory_info::log_usage("After Engine Initialization");
    rmi.barrier();
  } // end of synchronous engine
  






  template<typename VertexProgram>
  typename synchronous_engine<VertexProgram>::aggregator_type*
  synchronous_engine<VertexProgram>::get_aggregator() {
    return &aggregator;
  } // end of get_aggregator



  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::internal_stop() {
    for (size_t i = 0; i < rmi.numprocs(); ++i) 
      rmi.remote_call(i, &synchronous_engine<VertexProgram>::rpc_stop);
  } // end of internal_stop

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::rpc_stop() {
    force_abort = true;
  } // end of rpc_stop


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  signal(vertex_id_type gvid, const message_type& message) {
    rmi.barrier();
    internal_signal_rpc(gvid, message);
    rmi.barrier();
  } // end of signal



  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  signal_all(const message_type& message, const std::string& order) {
    for(lvid_type lvid = 0; lvid < graph.num_local_vertices(); ++lvid) {
      if(graph.l_is_master(lvid)) 
        internal_signal(vertex_type(graph.l_vertex(lvid)), message);
    }
  } // end of signal all
  


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  internal_signal(const vertex_type& vertex,
                  const message_type& message) {
    const lvid_type lvid = vertex.local_id();
    vlocks[lvid].lock();
    if( has_message.get(lvid) ) {
      messages[lvid] += message;
    } else {
      messages[lvid] = message;
      has_message.set_bit(lvid);
    }
    vlocks[lvid].unlock();       
  } // end of internal_signal


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  internal_signal_broadcast(vertex_id_type gvid, const message_type& message) {
    for (size_t i = 0; i < rmi.numprocs(); ++i) {
      if(i == rmi.procid()) internal_signal_rpc(gvid, message);
      else rmi.remote_call(i, &synchronous_engine<VertexProgram>::internal_signal_rpc,
                          gvid, message);
    }
  } // end of internal_signal_broadcast

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  internal_signal_rpc(vertex_id_type gvid,
                      const message_type& message) {
    if (graph.is_master(gvid)) {
      internal_signal(graph.vertex(gvid), message);
    }
  } // end of internal_signal_rpc



  

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  internal_post_delta(const vertex_type& vertex, const gather_type& delta) {
    const bool caching_enabled = !gather_cache.empty();
    if(caching_enabled) {
      const lvid_type lvid = vertex.local_id();      
      vlocks[lvid].lock();
      if( has_cache.get(lvid) ) {
        gather_cache[lvid] += delta;
      } else {
        gather_cache[lvid] = delta;
        has_cache.set_bit(lvid);
      }
      vlocks[lvid].unlock();       
    }
  } // end of post_delta


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  internal_clear_gather_cache(const vertex_type& vertex) {
    const bool caching_enabled = !gather_cache.empty();
    const lvid_type lvid = vertex.local_id();
    if(caching_enabled && has_cache.get(lvid)) {
      vlocks[lvid].lock();
      gather_cache[lvid] = gather_type();
      has_cache.clear_bit(lvid);
      vlocks[lvid].unlock();
    }
  } // end of clear_gather_cache
  



  template<typename VertexProgram>
  size_t synchronous_engine<VertexProgram>::
  num_updates() const { return completed_applys.value; }

  template<typename VertexProgram>
  float synchronous_engine<VertexProgram>::
  elapsed_seconds() const { return timer::approx_time_seconds() - start_time; }

  template<typename VertexProgram>
  int synchronous_engine<VertexProgram>::
  iteration() const { return iteration_counter; }



  template<typename VertexProgram>
  size_t synchronous_engine<VertexProgram>::total_memory_usage() const {
    size_t allocated_memory = memory_info::allocated_bytes();
    rmi.all_reduce(allocated_memory);
    return allocated_memory;
  } // compute the total memory usage of the GraphLab system


  template<typename VertexProgram> execution_status::status_enum
  synchronous_engine<VertexProgram>::start() {
    rmi.barrier();
    graph.finalize();
    // Initialization code ==================================================     
    // Reset event log counters? 
    rmi.dc().flush_counters();
    // Start the timer
    graphlab::timer timer; timer.start();
    start_time = timer::approx_time_seconds();
    iteration_counter = 0;
    force_abort = false;
    execution_status::status_enum termination_reason = 
      execution_status::UNSET; 
    // if (perform_init_vtx_program) {
    //   // Initialize all vertex programs
    //   run_synchronous( &synchronous_engine::initialize_vertex_programs );
    // }
    aggregator.start();
    rmi.barrier();
    // Program Main loop ====================================================      
    while(iteration_counter < max_iterations && 
          !force_abort ) {

      // Check first to see if we are out of time
      if(timeout != 0 && timeout < elapsed_seconds()) {
        termination_reason = execution_status::TIMEOUT;
        break;
      }


      logstream(LOG_EMPH) 
        << rmi.procid() << ": Starting iteration: " << iteration_counter 
        << std::endl;
      // Reset Active vertices ---------------------------------------------- 
      // Clear the active super-step and minor-step bits which will
      // be set upon receiving messages
      active_superstep.clear(); active_minorstep.clear();
      has_gather_accum.clear(); 
      rmi.barrier();

      // Exchange Messages --------------------------------------------------
      // Exchange any messages in the local message vectors
      run_synchronous( &synchronous_engine::exchange_messages );
      /**
       * Post conditions:
       *   1) only master vertices have messages
       */

      // Receive Messages ---------------------------------------------------
      // Receive messages to master vertices and then synchronize
      // vertex programs with mirrors if gather is required
      num_active_vertices = 0;
      run_synchronous( &synchronous_engine::receive_messages );
      has_message.clear();
      /**
       * Post conditions:
       *   1) there are no messages remaining
       *   2) All masters that received messages have their
       *      active_superstep bit set
       *   3) All masters and mirrors that are to participate in the
       *      next gather phases have their active_minorstep bit
       *      set.
       *   4) num_active_vertices is the number of vertices that
       *      received messages.
       */
      
      // Check termination condition  ---------------------------------------
      size_t total_active_vertices = num_active_vertices; 
      rmi.all_reduce(total_active_vertices);
      if (rmi.procid() == 0) 
        logstream(LOG_EMPH)
          << "\tActive vertices: " << total_active_vertices << std::endl;
      if(total_active_vertices == 0 ) {
        termination_reason = execution_status::TASK_DEPLETION;
        break;
      }


      // Execute gather operations-------------------------------------------
      // Execute the gather operation for all vertices that are active
      // in this minor-step (active-minorstep bit set).  
      run_synchronous( &synchronous_engine::execute_gathers );
      // Clear the minor step bit since only super-step vertices
      // (only master vertices are required to participate in the
      // apply step)
      active_minorstep.clear(); // rmi.barrier();
      /**
       * Post conditions:
       *   1) gather_accum for all master vertices contains the
       *      result of all the gathers (even if they are drawn from
       *      cache)
       *   2) No minor-step bits are set
       */

      // Execute Apply Operations -------------------------------------------
      // Run the apply function on all active vertices
      run_synchronous( &synchronous_engine::execute_applys );
      /**
       * Post conditions:
       *   1) any changes to the vertex data have been synchronized
       *      with all mirrors.
       *   2) all gather accumulators have been cleared
       *   3) If a vertex program is participating in the scatter
       *      phase its minor-step bit has been set to active (both
       *      masters and mirrors) and the vertex program has been
       *      synchronized with the mirrors.         
       */


      // Execute Scatter Operations -----------------------------------------
      // Execute each of the scatters on all minor-step active vertices.
      run_synchronous( &synchronous_engine::execute_scatters );
      /**
       * Post conditions:
       *   1) NONE
       */
      logstream(LOG_EMPH) << "\t Running Aggregators" << std::endl;
      // probe the aggregator
      aggregator.tick_synchronous();

      ++iteration_counter;
    }
    // Final barrier to ensure that all engines terminate at the same time
    rmi.full_barrier();
    // Stop the aggregator
    aggregator.stop();
    // return the final reason for termination
    return termination_reason;
  } // end of start




  // template<typename VertexProgram>
  // void synchronous_engine<VertexProgram>::
  // initialize_vertex_programs(size_t thread_id) {
  //   // For now we are using the engine as the context interface
  //   context_type context(*this, graph);
  //   for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices(); 
  //       lvid += threads.size()) {
  //     if(graph.l_is_master(lvid)) {          
  //       vertex_type vertex = local_vertex_type(graph.l_vertex(lvid));
  //       vertex_programs[lvid].init(context, vertex);
  //       // send the vertex program and vertex data to all mirrors
  //       sync_vertex_program(lvid); sync_vertex_data(lvid);
  //     }
  //     // recv_vertex_programs(); recv_vertex_data();
  //   }
  //   // Flush the buffer and finish receiving any remaining vertex
  //   // programs.
  //   thread_barrier.wait();
  //   if(thread_id == 0) { vprog_exchange.flush(); vdata_exchange.flush(); }
  //   thread_barrier.wait();
  //   recv_vertex_programs();
  //   recv_vertex_data();
  // } // end of initialize_vertex_programs


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  exchange_messages(size_t thread_id) {
    context_type context(*this, graph);
    for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices(); 
        lvid += threads.size()) {
      // if the vertex is not local and has a message send the
      // message and clear the bit
      if(!graph.l_is_master(lvid) && has_message.get(lvid)) {
        sync_message(lvid); 
        has_message.clear_bit(lvid);
        // clear the message to save memory
        messages[lvid] = message_type();
      }
    } // end of loop over vertices to send messages
      // Finish sending and receiving all messages
    thread_barrier.wait();
    if(thread_id == 0) message_exchange.flush(); 
    thread_barrier.wait();
    recv_messages();
  } // end of exchange_messages

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  receive_messages(size_t thread_id) {
    context_type context(*this, graph);
    for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices(); 
        lvid += threads.size()) {
      // if this is the master of lvid and we have a message
      if(graph.l_is_master(lvid) && has_message.get(lvid)) {
        // The vertex becomes active for this superstep 
        active_superstep.set_bit(lvid);
        ++num_active_vertices;       
        // Pass the message to the vertex program
        vertex_type vertex = vertex_type(graph.l_vertex(lvid));
        vertex_programs[lvid].recv_message(context, vertex,
                                           messages[lvid]);
        // clear the message to save memory
        messages[lvid] = message_type();
        // Determine if the gather should be run
        const vertex_program_type& const_vprog = vertex_programs[lvid];
        const vertex_type const_vertex = vertex;
        if(const_vprog.gather_edges(context, const_vertex) != 
           graphlab::NO_EDGES) {
          active_minorstep.set_bit(lvid);
          sync_vertex_program(lvid);
        }  
      }
      // recv_vertex_programs();
    }
    // Flush the buffer and finish receiving any remaining vertex
    // programs.
    thread_barrier.wait();
    if(thread_id == 0) vprog_exchange.flush();
    thread_barrier.wait();
    recv_vertex_programs();
  } // end of receive messages


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  execute_gathers(size_t thread_id) {
    context_type context(*this, graph);
    const bool caching_enabled = !gather_cache.empty();
    for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices(); 
        lvid += threads.size()) {
      // If this vertex is active in the gather minorstep
      if(active_minorstep.get(lvid)) {
        bool accum_is_set = false;
        gather_type accum = gather_type();         
        // if caching is enabled and we have a cache entry then use
        // that as the accum
        if( caching_enabled && has_cache.get(lvid) ) {
          accum = gather_cache[lvid];
          accum_is_set = true;
        } else {
          // recompute the local contribution to the gather
          const vertex_program_type& vprog = vertex_programs[lvid];
          local_vertex_type local_vertex = graph.l_vertex(lvid);
          const vertex_type vertex(local_vertex);
          const edge_dir_type gather_dir = vprog.gather_edges(context, vertex);
          // Loop over in edges
          if(gather_dir == IN_EDGES || gather_dir == ALL_EDGES) {
            foreach(local_edge_type local_edge, local_vertex.in_edges()) {
              edge_type edge(local_edge);
              if(accum_is_set) { // \todo hint likely                
                accum += vprog.gather(context, vertex, edge);
              } else {
                accum = vprog.gather(context, vertex, edge); 
                accum_is_set = true;
              }
            }
          } // end of if in_edges/all_edges
            // Loop over out edges
          if(gather_dir == OUT_EDGES || gather_dir == ALL_EDGES) {
            foreach(local_edge_type local_edge, local_vertex.out_edges()) {
              edge_type edge(local_edge);              
              if(accum_is_set) { // \todo hint likely
                accum += vprog.gather(context, vertex, edge);              
              } else {
                accum = vprog.gather(context, vertex, edge);
                accum_is_set = true;
              }
            }
          } // end of if out_edges/all_edges
          // If caching is enabled then save the accumulator to the
          // cache for future iterations.  Note that it is possible
          // that the accumulator was never set in which case we are
          // effectively "zeroing out" the cache.
          if(caching_enabled && accum_is_set) {              
            gather_cache[lvid] = accum; has_cache.set_bit(lvid); 
          } // end of if caching enabled            
        }
        // If the accum contains a value for the local gather we put
        // that estimate in the gather exchange.
        if(accum_is_set) sync_gather(lvid, accum);  
        if(!graph.l_is_master(lvid)) {
          // if this is not the master clear the vertex program
          vertex_programs[lvid] = vertex_program_type();
        }
      } // end of if active
  
    } // end of loop over vertices to compute gather accumulators
      // Finish sending and receiving all gather operations
    thread_barrier.wait();
    if(thread_id == 0) gather_exchange.flush();
    thread_barrier.wait();
    recv_gathers();
  } // end of execute_gathers


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  execute_applys(size_t thread_id) {
    context_type context(*this, graph);
    for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices(); 
        lvid += threads.size()) {
      // If this vertex is active on this super-step 
      if( active_superstep.get(lvid) ) {          
        // Only master vertices can be active in a super-step
        ASSERT_TRUE(graph.l_is_master(lvid));
        vertex_type vertex(graph.l_vertex(lvid));
        // Get the local accumulator.  Note that it is possible that
        // the gather_accum was not set during the gather.
        const gather_type& accum = gather_accum[lvid];
        vertex_programs[lvid].apply(context, vertex, accum);
        // record an apply as a completed task
        ++completed_applys;
        // Clear the accumulator to save some memory
        gather_accum[lvid] = gather_type();
        // synchronize the changed vertex data with all mirrors
        sync_vertex_data(lvid);  
        // determine if a scatter operation is needed
        const vertex_program_type& const_vprog = vertex_programs[lvid];
        const vertex_type const_vertex = vertex;
        if(const_vprog.scatter_edges(context, const_vertex) != 
           graphlab::NO_EDGES) {
          active_minorstep.set_bit(lvid);
          sync_vertex_program(lvid);
        } else { // we are done so clear the vertex program
          vertex_programs[lvid] = vertex_program_type();
        }
      } // end of if apply
    } // end of loop over vertices to run apply
      // Finish sending and receiving all changes due to apply operations
    thread_barrier.wait();
    if(thread_id == 0) { vprog_exchange.flush(); vdata_exchange.flush(); }
    thread_barrier.wait();
    recv_vertex_programs();
    recv_vertex_data();
  } // end of execute_applys

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  execute_scatters(size_t thread_id) {
    context_type context(*this, graph);
    for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices(); 
        lvid += threads.size()) {
      // If this vertex is active in the scatter minorstep
      if(active_minorstep.get(lvid)) {
        const vertex_program_type& vprog = vertex_programs[lvid];
        local_vertex_type local_vertex = graph.l_vertex(lvid);
        const vertex_type vertex(local_vertex);
        const edge_dir_type scatter_dir = vprog.scatter_edges(context, vertex);
        // Loop over in edges
        if(scatter_dir == IN_EDGES || scatter_dir == ALL_EDGES) {
          foreach(local_edge_type local_edge, local_vertex.in_edges()) {
            edge_type edge(local_edge);
            vprog.scatter(context, vertex, edge);
          }
        } // end of if in_edges/all_edges
          // Loop over out edges
        if(scatter_dir == OUT_EDGES || scatter_dir == ALL_EDGES) {
          foreach(local_edge_type local_edge, local_vertex.out_edges()) {
            edge_type edge(local_edge);
            vprog.scatter(context, vertex, edge);
          }
        } // end of if out_edges/all_edges
        // Clear the vertex program
        vertex_programs[lvid] = vertex_program_type();
      } // end of if active on this minor step
    } // end of loop over vertices to complete scatter operation
  } // end of execute_scatters


  // Data Synchronization ===================================================
  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  sync_vertex_program(lvid_type lvid) {
    ASSERT_TRUE(graph.l_is_master(lvid));
    const vertex_id_type vid = graph.global_vid(lvid);
    local_vertex_type vertex = graph.l_vertex(lvid);
    foreach(const procid_t& mirror, vertex.mirrors()) {
      vprog_exchange.send(mirror, 
                          std::make_pair(vid, vertex_programs[lvid]));
    }
  } // end of sync_vertex_program
  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  recv_vertex_programs() {
    procid_t procid(-1);
    typename vprog_exchange_type::buffer_type buffer;
    while(vprog_exchange.recv(procid, buffer)) {
      foreach(const vid_prog_pair_type& pair, buffer) {
        const lvid_type lvid = graph.local_vid(pair.first);
        ASSERT_FALSE(graph.l_is_master(lvid));
        vertex_programs[lvid] = pair.second;
        active_minorstep.set_bit(lvid);
      }
    }
  } // end of recv vertex programs


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  sync_vertex_data(lvid_type lvid) {
    ASSERT_TRUE(graph.l_is_master(lvid));
    const vertex_id_type vid = graph.global_vid(lvid);
    local_vertex_type vertex = graph.l_vertex(lvid);
    foreach(const procid_t& mirror, vertex.mirrors()) {
      vdata_exchange.send(mirror, std::make_pair(vid, vertex.data()));
    }
  } // end of sync_vertex_data

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  recv_vertex_data() {
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


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  sync_gather(lvid_type lvid, const gather_type& accum) {
    if(graph.l_is_master(lvid)) {
      vlocks[lvid].lock();
      if(has_gather_accum.get(lvid)) {
        gather_accum[lvid] += accum;
      } else {
        gather_accum[lvid] = accum;
        has_gather_accum.set_bit(lvid);
      }
      vlocks[lvid].unlock();
    } else {
      const procid_t master = graph.l_master(lvid);
      const vertex_id_type vid = graph.global_vid(lvid);
      gather_exchange.send(master, std::make_pair(vid, accum));
    }
  } // end of sync_gather

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  recv_gathers() {
    procid_t procid(-1);
    typename gather_exchange_type::buffer_type buffer;
    while(gather_exchange.recv(procid, buffer)) {
      foreach(const vid_gather_pair_type& pair, buffer) {
        const lvid_type lvid = graph.local_vid(pair.first);
        const gather_type& accum = pair.second;
        ASSERT_TRUE(graph.l_is_master(lvid));
        vlocks[lvid].lock();
        if( has_gather_accum.get(lvid) ) {
          gather_accum[lvid] += accum;
        } else {
          gather_accum[lvid] = accum;
          has_gather_accum.set_bit(lvid);
        }
        vlocks[lvid].unlock();
      }
    }
  } // end of recv_gather


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  sync_message(lvid_type lvid) {
    ASSERT_FALSE(graph.l_is_master(lvid));
    const procid_t master = graph.l_master(lvid);
    const vertex_id_type vid = graph.global_vid(lvid);
    message_exchange.send(master, std::make_pair(vid, messages[lvid]));
  } // end of send_message




  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  recv_messages() {
    procid_t procid(-1);
    typename message_exchange_type::buffer_type buffer;
    while(message_exchange.recv(procid, buffer)) {
      foreach(const vid_message_pair_type& pair, buffer) {
        const lvid_type lvid = graph.local_vid(pair.first);
        ASSERT_TRUE(graph.l_is_master(lvid));
        vlocks[lvid].lock();
        if( has_message.get(lvid) ) {
          messages[lvid] += pair.second;
        } else {
          messages[lvid] = pair.second;
          has_message.set_bit(lvid);
        }
        vlocks[lvid].unlock();
      }
    }
  } // end of recv_messages











}; // namespace


#include <graphlab/macros_undef.hpp>

#endif 

