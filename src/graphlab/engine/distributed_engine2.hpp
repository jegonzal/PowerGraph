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




/* \file iengine.hpp
   \brief The file containing the iengine description
   
   This file contains the description of the engine interface.  All
   graphlab engines (single_threaded, multi_threaded, distributed, ...)
   should satisfy the functionality described below.
*/

#ifndef GRAPHLAB_DISTRIBUTED_ENGINE_HPP
#define GRAPHLAB_DISTRIBUTED_ENGINE_HPP


#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/context/icontext.hpp>
#include <graphlab/update_functor/iupdate_functor.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>
#include <graphlab/util/chandy_misra.hpp>


namespace graphlab {
  

  template<typename Graph, typename UpdateFunctor>
  class distributed_engine: public iengine<Graph, UpdateFunctor> {
    
  public:
    // Include parent types
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef typename iengine_base::graph_type graph_type;
    typedef typename iengine_base::update_functor_type update_functor_type;
    
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef typename graph_type::edge_id_type   edge_id_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef typename graph_type::edge_type edge_type;

    typedef ischeduler<distributed_engine>      ischeduler_type;
    
    typedef typename iengine_base::icontext_type  icontext_type;
    typedef context<distributed_engine>         context_type;
    typedef context_manager<distributed_engine> context_manager_type;
   
   
    typedef typename iengine_base::termination_function_type termination_function_type;

    struct pipeline {
      p
    };
  private:
    dc_dist_object<distributed_engine<Graph, UpdateFunctor> > rmi;

    //! The local engine options
    graphlab_options opts; 

    graph_type& graph;
    chandy_misra locks;
    
    thread_group thrgroup;
    
    //! The scheduler
    ischeduler_type* scheduler_ptr;

  public:
    distributed_engine(distributed_control &dc, graph_type& graph): 
                                                        rmi(dc), graph(graph),
                                                        locks(graph){
      rmi.barrier();
    }

    //! Virtual destructor required for inheritance 
    ~distributed_engine() {
      rmi.barrier();
    }
    
    /**
     * \brief Force engine to terminate immediately.
     *
     * This function is used to stop the engine execution by forcing
     * immediate termination.  Any existing update tasks will finish
     * but no new update tasks will be started and the call to start()
     * will return.
     */
    void stop() { }

    
    /**
     * \brief Describe the reason for termination.
     *
     * Return the reason for the last termination.
     */
    execution_status::status_enum last_exec_status() const = 0;
   
    /**
     * \brief Get the number of updates executed by the engine.
     *
     * This function returns the numbe of updates executed by the last
     * run of this engine.
     * 
     * \return the total number of updates
     */
    size_t last_update_count() const = 0;
           
    /**
     * \brief Adds an update task with a particular priority.
     * This function is forwarded to the scheduler.
     */
    void schedule(vertex_id_type vid,
                          const update_functor_type& update_functor) {
      //.../
    }


    /**
     * \brief Creates a collection of tasks on all the vertices in the
     * graph, with the same update function and priority This function
     * is forwarded to the scheduler.
     */
    void schedule_all(const update_functor_type& update_functor) {
      scheduler_ptr->schedule_all(update_functor);
    }

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
    void add_termination_condition(termination_function_type term) { }

    //!  remove all associated termination functions
    void clear_termination_conditions() { };
    
    /**
     *  \brief The timeout is the total
     *  ammount of time in seconds that the engine may run before
     *  exeuction is automatically terminated.
     */
    void set_timeout(size_t timeout_secs) { };
    
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
    void set_task_budget(size_t max_tasks) { }


    /** \brief Update the engine options.  */
    void set_options(const graphlab_options& opts) {
      opts = new_opts;
    } 

    /** \brief get the current engine options. */
    const graphlab_options& get_options() {
      return opts;
    }
    
    /**
     * Unique to the distributed engine. The 
     */
    void initialize() {
      rmi.barrier();
      scheduler_ptr = scheduler_factory<distributed_engine>::
                          new_scheduler(opts.scheduler_type,
                          opts.scheduler_args,
                          graph,
                          opts.get_ncpus());

    }
    
    
    void thread_start() {
      
    }
    
    /**
     * \brief Start the engine execution.
     *
     * This \b blocking function starts the engine and does not
     * return until either one of the termination conditions evaluate
     * true or the scheduler has no tasks remaining.
     */
    void start() {
      rmi.barrier();
      thrgroup.launch(boost::function(...))
    }
  };

}

#endif // GRAPHLAB_DISTRIBUTED_ENGINE_HPP

