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


#ifndef GRAPHLAB_CORE_BASE_HPP
#define GRAPHLAB_CORE_BASE_HPP

#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/engine_options.hpp>

#include <graphlab/util/command_line_options.hpp>

#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/graph/graph.hpp>



#include <graphlab/metrics/metrics.hpp>
#include <graphlab/metrics/reporters/null_reporter.hpp>
#include <graphlab/metrics/reporters/basic_reporter.hpp>
#include <graphlab/metrics/reporters/file_reporter.hpp>
#include <graphlab/metrics/reporters/html_reporter.hpp>



namespace graphlab {


  /**
     \brief A GraphLab core is the base (or core) data structure in GraphLab.
     This provides a non type dependent base type for the core.
  */
  class core_base {
  public:
    /// default constructor does nothing
    core_base() { }
  private:
    //! Core is not copyable
    core_base(const core_base& other);
    //! Core is not copyable
    core_base& operator=(const core_base& other);

  public:


    virtual ~core_base() { }

    /**
     * \brief Set the type of scheduler.
     *
     * This will destroy the current engine and any tasks currently
     * associated with the scheduler.  See \ref Schedulers for the
     * list of supported schedulers.
     */
    virtual void set_scheduler_type(const std::string& scheduler_type) = 0;

    /**
     * \brief Set the scope consistency model used in this engine.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler.  The available scopes are:
     * 
     *  \li \b "full" This ensures full data consistency within the scope
     *  \li \b "edge" This ensures data consistency with just the
     *     vertex and edges
     *  \li \b "vertex" This ensures that a vertex cannot be updated
     *     by two processors simultaneously
     *  \li \b "none" This eliminates all locking 
     *
     * See \ref Scopes for details
     */
    virtual void set_scope_type(const std::string& scope_type) = 0;


    /**
     * \brief Set the engine type.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler. 
     *
     *  \li \b "async" This is the regular multithreaded engine
     *  \li \b "async_sim" This is a single threaded engine. But it can be 
     *                     be started with multiple "simulated threads".
     *                     The simulation is low-fidelity however, and should
     *                     be used with caution.
     */
    virtual void set_engine_type(const std::string& engine_type) = 0;
    
    /**
     * \brief Sets the output format of any recorded metrics
     *  \li \b "none" No reporting
     *  \li \b "basic" Outputs to screen
     *  \li \b "file" Outputs to a text file graphlab_metrics.txt
     *  \li \b "html" Outputs to a html file graphlab_metrics.html
     */
    virtual void set_metrics_type(const std::string& metrics_type) = 0;

    /**
       \brief Destroys a created engine (if any).
    */
    virtual void reset() = 0;
    
    /**
     * \brief Set the number of cpus that the engine will use.
     *
     * This will destroy the current engine and any tasks associated
     * with the current scheduler. 
     *
     */
    virtual void set_ncpus(size_t ncpus) = 0;


    /**
     * \brief Destroys and reconstructs the current engine,
     * reprocessing the engine arguments.  
     */
    virtual bool rebuild_engine() = 0;
    

    /**
     * \brief Set the engine options by passing in an engine options object.
     */
    virtual void set_engine_options(const engine_options& opts) = 0;

    virtual imetrics_reporter& get_reporter() = 0;
    

    /**
     * \brief Returns the engine options
     */
    virtual const engine_options& get_engine_options() const = 0;

    /**
     * \brief Returns a modifiable reference to the scheduler options
     */
    virtual scheduler_options& sched_options() = 0;
    

    /**
     * \brief Returns a constant reference to the scheduler options
     */
    virtual const scheduler_options& sched_options() const = 0;


    /**
     * \brief Set the engine options by simply parsing the command line
     * arguments. 
     */
    virtual bool parse_engine_options(int argc, char **argv) = 0;


    /**
     * \brief Run the engine until a termination condition is reached or
     * there are no more tasks remaining to execute.
     */
    virtual double start() = 0;
  

    /**
     * \brief Get the number of updates executed by the engine
     */
    virtual size_t last_update_count() = 0;
    
    
    virtual void fill_metrics() = 0;
    
    virtual void reset_metrics() = 0;
      
    /**
       \brief Outputs the recorded metrics
    */
    virtual void report_metrics() = 0;
    
    virtual void sync_now(glshared_base& shared) = 0;
  };
}

#endif

