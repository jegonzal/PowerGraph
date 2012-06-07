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
    typedef typename graph_type::vertex_id_type vertex_id_type;  
    typedef typename graph_type::vertex_type    vertex_type;
    typedef typename graph_type::edge_type      edge_type;

    typedef typename vertex_program_type::icontext_type icontext_type;   
    typedef distributed_aggregator<graph_type, icontext_type> aggregator_type;


    //! Virtual destructor required for inheritance 
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
     */
    virtual int iteration() const {
      return -1;
    }

     
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
     * default message is sent if no message is provided.
     */
    virtual void signal_all(const message_type& message,
                            const std::string& order = "sequential") = 0;
   

    
    /**
     * \copydoc distributed_aggregator::add_vertex_aggregator
     */
    template <typename ReductionType>
    bool add_vertex_aggregator
    (const std::string& key,
     boost::function<ReductionType(icontext_type&, vertex_type&)> map_function,
     boost::function<void(icontext_type&, const ReductionType&)> finalize_function) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->add_vertex_aggregator(key, map_function, finalize_function);
    } // end of add vertex aggregator



    
    /**
     * \copydoc distributed_aggregator::add_edge_aggregator
     */
    template <typename ReductionType>
    bool add_edge_aggregator
    (const std::string& key,
     boost::function<ReductionType(icontext_type&, edge_type&)> map_function,
      boost::function<void(icontext_type&, const ReductionType&)> finalize_function) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->add_edge_aggregator(key, map_function, finalize_function);
    } // end of add edge aggregator


    /**
     * \copydoc distributed_aggregator::aggregate_now
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
     * \copydoc distributed_aggregator::aggregate_periodic()
     */
    bool aggregate_periodic(const std::string& key, float seconds) {
      aggregator_type* aggregator = get_aggregator();
      if(aggregator == NULL) {
        logstream(LOG_FATAL) << "Aggregation not supported by this engine!" << std::endl;
        return false; 
      }
      return aggregator->aggregate_periodic(key, seconds);
    } // end of aggregate_periodic


  private:
    virtual aggregator_type* get_aggregator() = 0;
    
  }; // end of iengine interface

}

#endif

