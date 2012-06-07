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


/* \file iengine.hpp
   \brief The file containing the iengine description
   
   This file contains the description of the engine interface.  All
   graphlab engines (single_threaded, multi_threaded, distributed, ...)
   should satisfy the functionality described below.
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
     \brief The abstract interface of a GraphLab engine.
     The graphlab engine interface describes the core functionality
     provided by all graphlab engines.  The engine is templatized over
     the type of graph.
     
     The GraphLab engines are a core element of the GraphLab
     framework.  The engines are responsible for applying a the update
     tasks and sync operations to a graph and shared data using the
     scheduler to determine the update schedule. This class provides a
     generic interface to interact with engines written to execute on
     different platforms.
     
     While users are free to directly instantiate the engine of their
     choice we highly recommend the use of the \ref core data
     structure to manage the creation of engines. Alternatively, users
     can use the 
     \ref gl_new_engine "graphlab::engine_factory::new_engine"
     static functions to create
     engines directly from configuration strings.
  */
  template<typename VertexProgram>
  class iengine {
  public:
    typedef VertexProgram vertex_program_type;

    typedef typename vertex_program_type::message_type message_type;
    typedef typename vertex_program_type::icontext_type icontext_type;   
    typedef typename vertex_program_type::graph_type graph_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;  
    typedef typename graph_type::vertex_type    vertex_type;
    typedef typename graph_type::edge_type      edge_type;

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

