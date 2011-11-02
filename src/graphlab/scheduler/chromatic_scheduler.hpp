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


#ifndef GRAPHLAB_CHROMATIC_SCHEDULER_HPP
#define GRAPHLAB_CHROMATIC_SCHEDULER_HPP

#include <vector>

#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/cache_line_pad.hpp>


#include <graphlab/parallel/atomic.hpp>



#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/engine/terminator/iterminator.hpp>
#include <graphlab/engine/terminator/controlled_termination.hpp>
#include <graphlab/options/options_map.hpp>






namespace graphlab {

  /**
   * \ingroup group_schedulers
   *
   * Chromatic Scheduler
   */
  template<typename Engine>
  class chromatic_scheduler : public ischeduler<Engine> {
  public:

    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;
    typedef typename base::vertex_color_type    vertex_color_type;

    typedef controlled_termination terminator_type;



  private:    

    //! The update functor
    std::vector<update_functor_type> functors;
       

    std::vector< std::vector< vertex_id_type> > color_blocks;
    std::vector< cache_line_pad<size_t> > cpu_index;
    std::vector< cache_line_pad<size_t> > cpu_color;
    std::vector< cache_line_pad<size_t> > cpu_waiting;

    

    //! The maximum number of iterations
    size_t max_iterations;

    atomic<size_t> color;
    atomic<size_t> waiting;

    terminator_type term;

  public:

    chromatic_scheduler(const graph_type& graph, 
                        size_t ncpus,
                        const options_map& opts) :
      functors(ncpus), cpu_index(ncpus), cpu_color(ncpus), 
      cpu_waiting(ncpus), max_iterations(0) {
      color.value = 0;
      // Verify the coloring
      //ASSERT_TRUE(graph.valid_coloring());
      
      // parse the options
      opts.get_int_option("max_iterations", max_iterations);
      // Initialize the chromatic blocks
      vertex_color_type ncolors(0);
      for(vertex_id_type i = 0; i < graph.num_vertices(); ++i) 
        ncolors = std::max(ncolors, vertex_color_type(graph.color(i) + 1));
      color_blocks.resize(ncolors);
      for(vertex_id_type i = 0; i < graph.num_vertices(); ++i) 
        color_blocks[graph.color(i)].push_back(i);
    } // end of constructor




    
    /** Called by engine before executing the schedule */
    void start() {
      color.value = 0;
      // Initialize the cpu indexs
      for(size_t i = 0; i < cpu_index.size(); ++i) {
        cpu_index[i] = i;
        cpu_color[i] = -1;
        cpu_waiting[i] = true;
      }
      // Set waiting to zero
      waiting.value = 0;
      term.reset();
    }

    void schedule(vertex_id_type vid, 
                  const update_functor_type& fun) {  
      // Does nothing
    }

   
    void schedule_all(const update_functor_type& fun) {
      for(size_t i = 0; i < functors.size(); ++i) 
        functors[i] = fun;
    }
    

    /**
     * This function is called by the engine to ask for new work to
     * do.  The update task to be executed is returned in ret_task.
     *
     *  \retval NEWTASK There is an update task in ret_task to be
     *   executed
     *  \retval EMPTY Scheduler is empty
     */
    sched_status::status_enum get_next(size_t cpuid, 
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {
      // See if we are waiting
      if(cpu_waiting[cpuid].value) {
        // Nothing has changed so we are still waiting
        if(cpu_color[cpuid].value == color.value) return sched_status::EMPTY;
        // Otherwise color has changed so we reset and leave waiting
        // state
        cpu_color[cpuid].value = color.value;
        cpu_index[cpuid].value = cpuid;
        cpu_waiting[cpuid].value = false; 
      } else {      
        // Increment the index
        cpu_index[cpuid].value += cpu_index.size();
      }

      const size_t current_color = 
        cpu_color[cpuid].value % color_blocks.size();

      // Check to see that we have not run the maximum number of iterations
      if(max_iterations > 0 && 
         cpu_color[cpuid].value / color_blocks.size() >= max_iterations) {
        term.complete();
        return sched_status::EMPTY;
      }
      
      // If the index is good then execute it
      if(cpu_index[cpuid].value < color_blocks[ current_color ].size()) {
        const vertex_id_type vertex =
          color_blocks[ current_color ][ cpu_index[cpuid].value ];
        ret_vid = vertex;
        ret_fun = functors[cpuid];
        return sched_status::NEW_TASK;
      }
      
      // We overran so switch to waiting and increment the waiting counter
      size_t current_waiting = waiting.inc();
      cpu_waiting[cpuid].value = true;
      // If everyone is waiting reset and try again
      if(current_waiting == cpu_index.size()) {
        waiting.value = 0;
        color.inc();
        // Let the engine callback again
        return sched_status::EMPTY;
      }
      return sched_status::EMPTY;
    } // end of get_next_task
    
    /**
     * This is called after a task has been executed
     */
    void completed(size_t cpuid, 
                   vertex_id_type vid, 
                   const update_functor_type& task) { } 

   

    iterminator& terminator() { return term; };


    static void print_options_help(std::ostream &out) {
      out << "max_iterations = [integer, default = 0]\n";
      out << "update_function = [update_function_type,"
        "default = set on add_task]\n";
    };

   
    
    

 


  }; // End of chromatic scheduler

}
#endif

