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


#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/context/icontext.hpp>
#include <graphlab/update_functor/iupdate_functor.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
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
  template<typename Graph, typename UpdateFunctor>
  class iengine {
  public:

    //! The type of graph that the engine operates on
    typedef Graph graph_type;
    
    //! The type of the udpate functor
    typedef UpdateFunctor update_functor_type;

    //! The generic iupdate functor type
    typedef iupdate_functor<graph_type, update_functor_type> 
    iupdate_functor_type;


    //! The edge list type used by the graph
    typedef typename graph_type::edge_list_type  edge_list_type;

    //! The type of vertex color used by the graph
    typedef typename graph_type::vertex_color_type vertex_color_type;


    //! The type of scheduler
    typedef ischeduler<iengine> ischeduler_type;

    //! The type of context 
    typedef icontext<graph_type, update_functor_type> icontext_type;

    
    /**
     * The termination function is a function that reads the shared
     * data and returns true if the engine should terminate execution.
     * The termination function is called at fixed millisecond
     * intervals and therefore the engine may continue to execute even
     * after a termination function evaluates to true.  Because
     * termination functions are executed frequently and cannot
     * directly contribut to the computation, they should return
     * quickly.
     */
    typedef bool (*termination_function_type) ();
    

    //! Virtual destructor required for inheritance 
    virtual ~iengine() {};
    
    /**
     * \brief Start the engine execution.
     *
     * This \b blocking function starts the engine and does not
     * return until either one of the termination conditions evaluate
     * true or the scheduler has no tasks remaining.
     */
    virtual void start() = 0;


    /**
     * \brief Force engine to terminate immediately.
     *
     * This function is used to stop the engine execution by forcing
     * immediate termination.  Any existing update tasks will finish
     * but no new update tasks will be started and the call to start()
     * will return.
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
    virtual size_t last_update_count() const = 0;
           
    /**
     * \brief Adds an update task with a particular priority.
     * This function is forwarded to the scheduler.
     */
    virtual void schedule(vertex_id_type vid,
                          const update_functor_type& update_functor) = 0;


    /**
     * \brief Creates a collection of tasks on all the vertices in the
     * graph, with the same update function and priority This function
     * is forwarded to the scheduler.
     */
    virtual void schedule_all(const update_functor_type& update_functor) = 0;

    /**
     * Schedule an update on all the neighbors of a particular vertex
     */
    virtual void schedule_in_neighbors(const vertex_id_type& vertex, 
                                       const update_functor_type& update_fun) = 0;

    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    virtual void schedule_out_neighbors(const vertex_id_type& vertex, 
                                        const update_functor_type& update_fun) = 0;
                                                  
    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    virtual void schedule_neighbors(const vertex_id_type& vertex, 
                                    const update_functor_type& update_fun) = 0;


    /**
     * \brief associate a termination function with this engine.
     *
     * An engine can typically have many termination functions
     * associated with it. A termination function is a function which
     * takes a constant reference to the shared data and returns a
     * boolean which is true if the engine should terminate execution.
     *
     * A termination function has the following type:
     * \code
     * bool term_fun(const ishared_data_type* shared_data)
     * \endcode
     */
    virtual void add_termination_condition(termination_function_type term) = 0;

    //!  remove all associated termination functions
    virtual void clear_termination_conditions() = 0;
    
    /**
     *  \brief The timeout is the total
     *  ammount of time in seconds that the engine may run before
     *  exeuction is automatically terminated.
     */
    virtual void set_timeout(size_t timeout_secs) = 0;
    
    /**
     * Get the elapsed time since start was called in milliseconds
     */
    virtual size_t elapsed_time() const = 0;


    /**
     * \brief set a limit on the number of tasks that may be executed.
     * 
     * By once the engine has achived the max_task parameter execution
     * will be terminated. If max_tasks is set to zero then the
     * task_budget is ignored.  If max_tasks is greater than zero than
     * the value of max tasks is used.  Note that if max_task is
     * nonzero the engine encurs the cost of an additional atomic
     * operation in the main loop potentially reducing the overall
     * parallel performance.
     */
    virtual void set_task_budget(size_t max_tasks) = 0;


    /** \brief Update the engine options.  */
    virtual void set_options(const graphlab_options& opts) = 0;

    /** \brief get the current engine options. */
    virtual const graphlab_options& get_options() = 0;


    //! Get the global data and lock
    virtual void get_global(const std::string& key,      
                            graphlab::any_vector*& ret_values_ptr,
                            bool& ret_is_const) = 0; 

    //! Get the global data and lock
    virtual void acquire_global_lock(const std::string& key,      
                                     size_t index = 0) = 0;
    //! Release the global data lock
    virtual void release_global_lock(const std::string& key,      
                                     size_t index = 0) = 0;

    
    
  };

}

#endif

