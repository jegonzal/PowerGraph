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


#include <graphlab/vertex_program/icontext.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/options/graphlab_options.hpp>




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
    typedef typename vertex_program_type::graph_type graph_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;  
    typedef typename graph_type::vertex_type vertex_type;


    //! Virtual destructor required for inheritance 
    virtual ~iengine() {};
    
    /**
     * \brief Start the engine execution.
     *
     */
    virtual void start(bool perform_init_vtx_program = true) = 0;


    /**
     * \brief Force engine to terminate immediately.
     *
     */
    virtual void stop() = 0;

    
    /**
     * \brief Describe the reason for termination.
     *
     * Return the reason for the last termination.
     */
    virtual execution_status::status_enum last_exec_status() const = 0;
   
    /**
     * \brief Get the number of updates executed by the engine.
     *
     * This function returns the numbe of updates executed by the last
     * run of this engine.
     * 
     * \return the total number of updates
     */
    virtual size_t num_updates() const = 0;
           
    /**
     * \brief Signals a vertex with an optional message
     * 
     * Signals a vertex and schedules it to be executed in the future
     * Must be called on all machines simultaneously.
     */
    virtual void signal(vertex_id_type vertex,
                        const message_type& message = message_type()) = 0;

    /**
     * \brief Signals a vertex with an optional message
     * 
     * Signals a vertex, and schedules it to be executed in the future.
     * must be called on a vertex accessible by the current machine.
     */
    virtual void signal_internal(const vertex_type& vertex,
                                 const message_type& message = message_type()) = 0;


                                 
    /**
     * \brief Signals a global vid with an optional message. Ignored if current
     *        machine does not own it
     * 
     * Signals a global vid, and schedules it to be executed in the future.
     * If current machine does not contain the vertex, it is ignored.
     */
    virtual void signal_internal_gvid(vertex_id_type gvid,
                                 const message_type& message = message_type()) = 0;


    /**
     * \brief Signals a global vid with an optional message. 
     * Signals a global vid, and schedules it to be executed in the future.
     * Broadcast to all machines to guarantee gvid is signaled.
     */
    virtual void signal_broadcast(vertex_id_type gvid,
                                  const message_type& message = message_type()) = 0;

                                 
    /**
     * \brief Send a message to all vertices
     */
    virtual void signal_all(const message_type& message,
                            const std::string& order = "sequential") = 0;
    /**
     * Get the elapsed time since start was called in milliseconds
     */
    virtual size_t elapsed_time() const = 0;


    /** \brief get the current engine options. */
    //    virtual const graphlab_options& get_options() = 0;



    
    
  };

}

#endif

