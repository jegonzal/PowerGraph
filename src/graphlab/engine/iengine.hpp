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



/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */

#ifndef GRAPHLAB_IENGINE_HPP
#define GRAPHLAB_IENGINE_HPP

#include <boost/bind.hpp>
#include <boost/functional.hpp>

#include <graphlab/vertex_program/icontext.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/aggregation/distributed_aggregator.hpp>



namespace graphlab {
  

  /**
   * \ingroup engine
   *
   * \brief The abstract interface of a GraphLab engine.  
   * 
   * A GraphLab engine is responsible for executing vertex programs in
   * parallel on one or more machines.  GraphLab has a collection of
   * different engines with different guarantees on how
   * vertex-programs are executed.  However each engine must implement
   * the iengine interface to allow them to be used "interchangeably."
   *
   * In addition to executing vertex programs GraphLab engines also
   * expose a synchronous aggregation framework. This allows users to
   * attach "map-reduce" style jobs that are run periodically on all
   * edges or vertices while GraphLab programs are actively running.
   *
   * Example Usage
   * =================
   *
   * One can use the iengine interface to select between different
   * engines at runtime:
   *
   * \code
   * iengine<pagerank>* engine_ptr = NULL;
   * if(cmdline_arg == "synchronous") {
   *   engine_ptr = new synchronous_engine<pagerank>(dc, graph, cmdopts);  
   * } else {
   *   engine_ptr = new async_consistent_engine<pagerank>(dc, graph, cmdopts);  
   * }
   * // Attach an aggregator
   * engine_ptr->add_edge_aggregator<float>("edge_map", 
   *                                        edge_map_fun, finalize_fun);
   * // Make it run every 3 seconds
   * engine_ptr->aggregate_periodic("edge_map");
   * // Signal all vertices
   * engine_ptr->signal_all();
   * // Run the engine
   * engine_ptr->start();
   * // do something interesting
   * delete engine_ptr; engine_ptr = NULL;
   * \endcode  
   *
   * @tparam VertexProgram The user defined vertex program which should extend the
   * \ref ivertex_program interface.
   */
  template<typename VertexProgram>
  class iengine {
  public:
    /**
     * \brief The user defined vertex program type which should extend
     * ivertex_program.
     */
    typedef VertexProgram vertex_program_type;

    /**
     * \brief The user defined message type which is defined in
     * ivertex_program::message_type. 
     *
     */
    typedef typename vertex_program_type::message_type message_type;

    /**
     * \brief The graph type which is defined in
     * ivertex_program::graph_type and will typically be
     * \ref distributed_graph.
     */
    typedef typename vertex_program_type::graph_type graph_type;

    /**
     * \brief The vertex identifier type defined in 
     * \ref graphlab::vertex_id_type.
     */
    typedef typename graph_type::vertex_id_type vertex_id_type;  

    /**
     * \brief the vertex object type which contains a reference to the
     * vertex data and is defined in the iengine::graph_type 
     * (see for example \ref distributed_graph::vertex_type).
     */
    typedef typename graph_type::vertex_type    vertex_type;

    /**
     * \brief the edge object type which contains a reference to the
     * edge data and is defined in the iengine::graph_type (see for
     * example \ref distributed_graph::edge_type).
     */
    typedef typename graph_type::edge_type      edge_type;

    /**
     * \brief The context type which is passed into vertex programs as
     * a callback to the engine.  
     *
     * Most engines use the \ref graphlab::context implementation.
     */
    typedef typename vertex_program_type::icontext_type icontext_type;

    /**
     * \brief The type of the distributed aggregator used by each engine to
     * implement distributed aggregation.
     */   
    typedef distributed_aggregator<graph_type, icontext_type> aggregator_type;


    /**
     * \brief Virtual destructor required for inheritance
     */ 
    virtual ~iengine() {};
    
    /**
     * \brief Start the engine execution.
     *
     * Behavior details depend on the engine implementation. See the
     * implementation documentation for specifics.
     * 
     * \param [in] perform_init_vertex_program If true, runs init on each
     * vertex program before any message processing happens.
     * Defaults to true.
     * 
     * @return the reason for termination
     */
    virtual execution_status::status_enum 
    start(bool perform_init_vtx_program = true) = 0;
   
    /**
     * \brief Compute the total number of updates (calls to apply)
     * executed since start was last invoked.
     *
     * \return Total number of updates
     */
    virtual size_t num_updates() const = 0;

    /**
     * \brief Get the elapsed time in seconds since start was last
     * called.
     * 
     * \return elapsed time in seconds
     */
    virtual float elapsed_seconds() const = 0;

    /**
     * \brief get the current iteration number.  This is not defined
     * for all engines in which case -1 is returned.
     *
     * \return the current iteration or -1 if not supported.
     */
    virtual int iteration() const { return -1; }

     
    /**
     * \brief Signals single a vertex with an optional message.
     * 
     * This function sends a message to particular vertex which will
     * receive that message on start. The signal function must be
     * invoked on all machines simultaneously.  For example:
     *
     * \code
     * graphlab::synchronous_engine<vprog> engine(dc, graph, opts);
     * engine.signal(0); // signal vertex zero
     * \endcode
     *
     * and _not_:
     *
     * \code
     * graphlab::synchronous_engine<vprog> engine(dc, graph, opts);
     * if(dc.procid() == 0) engine.signal(0); // signal vertex zero
     * \endcode
     *
     * Since signal is executed synchronously on all machines it
     * should only be used to schedule a small set of vertices. The
     * preferred method to signal a large set of vertices (e.g., all
     * vertices that are a certain type) is to use either the vertex
     * program init function or the aggregation framework.  For
     * example to signal all vertices that have a particular value one
     * could write:
     *
     * \code
     * struct bipartite_opt : 
     *   public graphlab::ivertex_program<graph_type, gather_type> {
     *   // The user defined init function
     *   void init(icontext_type& context, vertex_type& vertex) {
     *     // Signal myself if I am a certain type
     *     if(vertex.data().on_left) context.signal(vertex);
     *   }
     *   // other vastly more interesting code
     * };
     * \endcode
     *
     * @param [in] vid the vertex id to signal
     * @param [in] message the message to send to that vertex.  The
     * default message is sent if no message is provided. 
     * (See ivertex_program::message_type for details about the
     * message_type). 
     */
    virtual void signal(vertex_id_type vertex,
                        const message_type& message = message_type()) = 0;
    
    /**
     * \brief Signal all vertices with a particular message.
     * 
     * This function sends the same message to all vertices which will
     * receive that message on start. The signal_all function must be
     * invoked on all machines simultaneously.  For example:
     *
     * \code
     * graphlab::synchronous_engine<vprog> engine(dc, graph, opts);
     * engine.signal_all(); // signal all vertices
     * \endcode
     *
     * and _not_:
     *
     * \code
     * graphlab::synchronous_engine<vprog> engine(dc, graph, opts);
     * if(dc.procid() == 0) engine.signal_all(); // signal vertex zero
     * \endcode
     *
     * The signal_all function is the most common way to send messages
     * to the engine.  For example in the pagerank application we want
     * all vertices to be active on the first round.  Therefore we
     * would write:
     *
     * \code
     * graphlab::synchronous_engine<pagerank> engine(dc, graph, opts);
     * engine.signal_all();
     * engine.start();
     * \endcode
     *
     * @param [in] message the message to send to all vertices.  The
     * default message is sent if no message is provided
     * (See ivertex_program::message_type for details about the
     * message_type). 
     */
    virtual void signal_all(const message_type& message = message_type(),
                            const std::string& order = "shuffle") = 0;
   

    
    /** 
     * \brief Creates a vertex aggregator. Returns true on success.
     *        Returns false if an aggregator of the same name already
     *        exists.
     *
     * Creates a vertex aggregator associated to a particular key.
     * The map_function is called over every vertex in the graph, and the
     * return value of the map is summed. The finalize_function is then called
     * on the result of the reduction.
     *
     * \tparam ReductionType The output of the map function. Must be summable
     *                       and \ref Serializable.
     * \param [in] map_function The Map function to use. Must take an
     *                          icontext_type& as its first argument, and an
     *                          vertex_type& as its second argument. Returns a
     *                          ReductionType which must be summable and
     *                          \ref Serializable
     * \param [in] finalize_function The Finalize function to use. Must take
     *                               an icontext_type& as its first argument
     *                               and a const ReductionType& as its second
     *                               argument.
     *
     * \warning Pay attention to the types! A slightly erroneous type
     *          can produce screens of errors.
     */
    template <typename ReductionType,
              typename VertexMapType,
              typename FinalizerType>
    bool add_vertex_aggregator(const std::string& key,
                               const VertexMapType& map_function,
                               const FinalizerType& finalize_function) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->add_vertex_aggregator<ReductionType>(key, map_function, finalize_function);
    } // end of add vertex aggregator

#if defined(__cplusplus) && __cplusplus >= 201103L
    template <typename VertexMapType,
              typename FinalizerType>
    bool add_vertex_aggregator(const std::string& key,
                               const VertexMapType& map_function,
                               const FinalizerType& finalize_function) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->add_vertex_aggregator(key, map_function, finalize_function);
    } // end of add vertex aggregator

#endif
    
    /** \brief Creates a edge aggregator. Returns true on success.
     *         Returns false if an aggregator of the same name already exists
     *
     * Creates an edge aggregator associated to a particular key.
     * The map_function is called over every edge in the graph, and the return
     * value of the map is summed. The finalize_function is then called on
     * the result of the reduction.
     * 
     * \tparam ReductionType The output of the map function. Must be summable
     *                       and \ref Serializable.
     * \param [in] map_function The Map function to use. Must take an
     *                          icontext_type& as its first argument, and an
     *                          edge_type& as its second argument. Returns a
     *                          ReductionType which must be summable and
     *                          \ref Serializable
     * \param [in] finalize_function The Finalize function to use. Must take
     *                               an icontext_type& as its first argument
     *                               and a const ReductionType& as its second
     *                               argument.
     *
     * \warning Pay attention to the types! A slightly erroneous type
     *          can produce screens of errors
     */
    template <typename ReductionType,
              typename EdgeMapType,
              typename FinalizerType>
    bool add_edge_aggregator(const std::string& key,
                               const EdgeMapType& map_function,
                               const FinalizerType& finalize_function) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->add_edge_aggregator<ReductionType>(key, map_function, finalize_function);
    } // end of add edge aggregator


#if defined(__cplusplus) && __cplusplus >= 201103L
    template <typename EdgeMapType,
              typename FinalizerType>
    bool add_edge_aggregator(const std::string& key,
                               const EdgeMapType& map_function,
                               const FinalizerType& finalize_function) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->add_edge_aggregator(key, map_function, finalize_function);
    } // end of add edge aggregator
#endif

    /**
     * Performs an immediate aggregation on a key. All machines must
     * call this simultaneously. If the key is not found,
     * false is returned. Otherwise return true on success.
     *
     * \param[in] key Key to aggregate now
     * \return False if key not found, True on success.
     */
    bool aggregate_now(const std::string& key) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->aggregate_now(key);
    } // end of aggregate_now


    /**
     * Requests that the aggregator with a given key be aggregated
     * every certain number of seconds when the engine is running.
     * Note that the period is prescriptive: in practice the actual
     * period will be larger than the requested period. 
     * Seconds must be >= 0;
     *
     * \param [in] key Key to schedule
     * \param [in] seconds How frequently to schedule. Must be >= 0. In the
     *                     synchronous engine, seconds == 0 will ensure that
     *                     this key is recomputed every iteration.
     * 
     * All machines must call simultaneously.
     * \return Returns true if key is found and seconds >= 0,
     *         and false otherwise.
     */
    bool aggregate_periodic(const std::string& key, float seconds) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->aggregate_periodic(key, seconds);
    } // end of aggregate_periodic



    /**
     * \brief This is used by iengine to get the 
     * \ref distributed_aggregator from the derived class to support
     * the local templated aggregator interface. 
     *
     * \return a pointer to the distributed aggregator for that
     * engine. If no aggregator is available or aggregation is not
     * supported then return NULL.
     */
    virtual aggregator_type* get_aggregator() = 0;
    
  }; // end of iengine interface

} // end of namespace graphlab

#endif

