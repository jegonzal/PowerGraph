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
  
  

  template<typename VertexProgram>
  class synchronous_engine : 
    public iengine<VertexProgram> {

  public:

    typedef VertexProgram vertex_program_type;
    typedef typename VertexProgram::gather_type gather_type;
    typedef typename VertexProgram::message_type message_type;
    typedef typename VertexProgram::vertex_data_type vertex_data_type;
    typedef typename VertexProgram::edge_data_type edge_data_type;

    typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;

    typedef typename graph_type::vertex_type          vertex_type;
    typedef typename graph_type::edge_type            edge_type;

    typedef typename graph_type::local_vertex_type    local_vertex_type;
    typedef typename graph_type::local_edge_type      local_edge_type;
    typedef typename graph_type::lvid_type            lvid_type;




    typedef icontext<vertex_type, gather_type, message_type> icontext_type;
    typedef context<synchronous_engine> context_type;
       
  private:


    dc_dist_object< synchronous_engine<VertexProgram> > rmi;

    //! The local engine options
    graphlab_options opts; 
    //! the distributed graph
    graph_type& graph;
    //! local threads object
    thread_pool threads;
    graphlab::barrier thread_barrier;

    size_t max_iterations;
    size_t iteration_counter;

    //! Engine state
    bool started;

    graphlab::timer timer;
    float start_time;


    std::vector<mutex>    vlocks;
    std::vector<vertex_program_type> vertex_programs;

    std::vector<message_type> messages;
    dense_bitset              has_message;
 

    std::vector<gather_type>  gather_accum;
    dense_bitset              has_gather_accum;


    std::vector<gather_type>  gather_cache;
    dense_bitset              has_cache;   

    dense_bitset             active_superstep;
    atomic<size_t>           num_active_vertices;
    dense_bitset             active_minorstep;      

    atomic<size_t> completed_tasks;
   
    execution_status::status_enum termination_reason; 

    // Exchange used to swap vertex programs
    typedef std::pair<vertex_id_type, vertex_program_type> vid_prog_pair_type;
    typedef buffered_exchange<vid_prog_pair_type> vprog_exchange_type;
    vprog_exchange_type vprog_exchange;
    // Exchange used to swap vertex data between machines
    typedef std::pair<vertex_id_type, vertex_data_type> vid_vdata_pair_type;
    typedef buffered_exchange<vid_vdata_pair_type> vdata_exchange_type;
    vdata_exchange_type vdata_exchange;
    // Exchange used to transfer gather data
    typedef std::pair<vertex_id_type, gather_type> vid_gather_pair_type;
    typedef buffered_exchange<vid_gather_pair_type> gather_exchange_type;
    gather_exchange_type gather_exchange;
    // Exchange used to transfer message data
    typedef std::pair<vertex_id_type, gather_type> vid_message_pair_type;
    typedef buffered_exchange<vid_message_pair_type> message_exchange_type;
    message_exchange_type message_exchange;


    
  public:
    synchronous_engine(distributed_control &dc, 
                       graph_type& graph,
                       size_t ncpus);
    void initialize();
    size_t total_memory_usage() const;

    void start();
    void stop() { /* implement */ }
    execution_status::status_enum last_exec_status() const;
    size_t num_updates() const;
    void signal(const vertex_type& vertex,
                const message_type& message = message_type());    
    void signal_all(const message_type& message = message_type(),
                    const std::string& order = "sequential");

    float elapsed_seconds() const;
    size_t iteration() const; 

    void set_timeout(size_t timeout_secs) { /* implement */ }
    size_t elapsed_time() const { return 0; /* implement */ }
    void set_task_budget(size_t max_tasks) { /* implement */ }
    void set_options(const graphlab_options& opts) { /* implement */ }

    void post_delta(const vertex_type& vertex,
                    const gather_type& delta);
    void clear_gather_cache(const vertex_type& vertex);

  private:

    // Program Steps ==========================================================
    template<typename MemberFunction>       
    void run_synchronous(MemberFunction member_fun) {
      // launch the initialization threads
      for(size_t i = 0; i < threads.size(); ++i) 
        threads.launch(boost::bind(member_fun, this, i));
      // Wait for all threads to finish
      threads.join();
      rmi.barrier();
    } // end of run_synchronous

    void initialize_vertex_programs(size_t thread_id);
    void exchange_messages(size_t thread_id);
    void receive_messages(size_t thread_id);
    void execute_gathers(size_t thread_id);
    void execute_applys(size_t thread_id);
    void execute_scatters(size_t thread_id);

    // Data Synchronization ===================================================
    void sync_vertex_program(lvid_type lvid);
    void recv_vertex_programs();
    void sync_vertex_data(lvid_type lvid);
    void recv_vertex_data();
    void sync_gather(lvid_type lvid, const gather_type& accum);
    void recv_gathers();
    void sync_message(lvid_type lvid);
    void recv_messages();

  }; // end of class synchronous engine








  template<typename VertexProgram>
  synchronous_engine<VertexProgram>::
  synchronous_engine(distributed_control &dc, 
                     graph_type& graph,
                     size_t ncpus) : 
    rmi(dc, this), graph(graph), threads(ncpus), thread_barrier(ncpus), 
    max_iterations(-1), iteration_counter(0),
    vprog_exchange(dc), vdata_exchange(dc), 
    gather_exchange(dc), message_exchange(dc) {
    rmi.barrier();
  } // end of synchronous engine
  


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  signal(const vertex_type& vertex,
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
  } // end of send message


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  signal_all(const message_type& message, const std::string& order) {
    for(lvid_type lvid = 0; lvid < graph.num_local_vertices(); ++lvid) {
      if(graph.l_is_master(lvid)) 
        signal(vertex_type(graph.l_vertex(lvid)), message);
    }
  } // end of send message
  
  

  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  post_delta(const vertex_type& vertex, const gather_type& delta) {
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
  clear_gather_cache(const vertex_type& vertex) {
    const bool caching_enabled = !gather_cache.empty();
    if(caching_enabled) {
      const lvid_type lvid = vertex.local_id();      
      vlocks[lvid].lock();
      gather_cache[lvid] = gather_type();
      has_cache.clear_bit(lvid);
      vlocks[lvid].unlock();
    }
  } // end of clear_gather_cache
  




  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::initialize() {
    graph.finalize();
    if (rmi.procid() == 0) 
      memory_info::print_usage("Before Engine Initialization");
    logstream(LOG_INFO) 
      << rmi.procid() << ": Initializing..." << std::endl;
    vlocks.resize(graph.num_local_vertices());
    vertex_programs.resize(graph.num_local_vertices());

    messages.resize(graph.num_local_vertices());
    has_message.resize(graph.num_local_vertices());

    gather_accum.resize(graph.num_local_vertices());
    has_gather_accum.resize(graph.num_local_vertices());

    gather_cache.resize(graph.num_local_vertices());
    has_cache.resize(graph.num_local_vertices());

    active_superstep.resize(graph.num_local_vertices());
    active_superstep.clear();
    active_minorstep.resize(graph.num_local_vertices());
    active_minorstep.clear();
    if (rmi.procid() == 0) 
      memory_info::print_usage("After Engine Initialization");
    rmi.barrier();
  } // end of initialize



  template<typename VertexProgram>
  execution_status::status_enum synchronous_engine<VertexProgram>::
  last_exec_status() const { return termination_reason;  }



  template<typename VertexProgram>
  size_t synchronous_engine<VertexProgram>::
  num_updates() const { return completed_tasks.value; }

  template<typename VertexProgram>
  float synchronous_engine<VertexProgram>::
  elapsed_seconds() const { return lowres_time_seconds() - start_time; }

  template<typename VertexProgram>
  size_t synchronous_engine<VertexProgram>::
  iteration() const { return iteration_counter; }



  template<typename VertexProgram>
  size_t synchronous_engine<VertexProgram>::total_memory_usage() const {
    size_t allocated_memory = memory_info::allocated_bytes();
    rmi.all_reduce(allocated_memory);
    return allocated_memory;
  } // compute the total memory usage of the GraphLab system


  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::start() {
    // Initialization code ==================================================     
    // Reset event log counters? 
    rmi.dc().flush_counters();
    // Start the timer
    timer.start();
    start_time = lowres_time_seconds();
    iteration_counter = 0;

    // Initialize all vertex programs
    run_synchronous( &synchronous_engine::initialize_vertex_programs );

    // Program Main loop ====================================================      
    for (iteration_counter = 0; iteration_counter < max_iterations; 
         ++iteration_counter) {
      if (rmi.procid() == 0) 
        std::cout << "Starting iteration: " << iteration_counter << std::endl;
      // Reset Active vertices ----------------------------------------------

      // Clear the active super-step and minor-step bits which will
      // be set upon receiving messages
      active_superstep.clear(); active_minorstep.clear();
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
      has_gather_accum.clear();  // rmi.barrier();
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
    }
    // Final barrier to ensure that all engines terminate at the same time
    rmi.barrier();
  } // end of start




  template<typename VertexProgram>
  void synchronous_engine<VertexProgram>::
  initialize_vertex_programs(size_t thread_id) {
    // For now we are using the engine as the context interface
    context_type context(*this, graph);
    for(lvid_type lvid = thread_id; lvid < graph.num_local_vertices(); 
        lvid += threads.size()) {
      if(graph.l_is_master(lvid)) {          
        vertex_type vertex = local_vertex_type(graph.l_vertex(lvid));
        vertex_programs[lvid].init(context, vertex);
        // send the vertex program and vertex data to all mirrors
        sync_vertex_program(lvid); sync_vertex_data(lvid);
      }
      // recv_vertex_programs(); recv_vertex_data();
    }
    // Flush the buffer and finish receiving any remaining vertex
    // programs.
    thread_barrier.wait();
    if(thread_id == 0) { vprog_exchange.flush(); vdata_exchange.flush(); }
    thread_barrier.wait();
    recv_vertex_programs();
    recv_vertex_data();
  } // end of initialize_vertex_programs


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
  } // end of execute_applys

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
        gather_type accum; //! This machines contribution to the gather          
        // if caching is enabled and we have a cache entry then use
        // that as the accum
        if( caching_enabled || has_cache.get(lvid) ) {
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
          if(caching_enabled) {              
            gather_cache[lvid] = accum; has_cache.set_bit(lvid); 
          } // end of if caching enabled            
        }
        // If the accum contains a value for the local gather we put
        // that estimate in the gather exchange.
        if(accum_is_set) sync_gather(lvid, accum);
      }
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
        ++completed_tasks;
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
    gather_exchange.send(master, std::make_pair(vid, messages[lvid]));
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

