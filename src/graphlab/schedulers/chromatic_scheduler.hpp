/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_CHROMATIC_SCHEDULER_HPP
#define GRAPHLAB_CHROMATIC_SCHEDULER_HPP

#include <vector>

#include <graphlab/logger/assertions.hpp>

#include <graphlab/tasks/update_task.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/schedulers/icallback.hpp>
#include <graphlab/util/controlled_termination.hpp>

#include <graphlab/parallel/atomic.hpp>
#include <graphlab/schedulers/support/unused_scheduler_callback.hpp>


namespace graphlab {

  /**
   * \ingroup group_schedulers
   *
   * Chromatic Scheduler
   */
  template<typename Graph>
  class chromatic_scheduler : public ischeduler<Graph> {
  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef typename base::iengine_type         iengine_type;
    typedef typename base::update_task_type     update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type        callback_type;
    typedef typename base::monitor_type         monitor_type;
    typedef controlled_termination terminator_type;

    
    /**
     * Used to prevent false cache sharing by padding T
     */
    template <typename T> struct cache_line_pad  {
      T value;
      char pad[64 - (sizeof(T) % 64)];      
      cache_line_pad(const T& value = T()) : value(value) { }
      T& operator=(const T& other) { return value = other; }
      
    };

    
    chromatic_scheduler(iengine_type* engine,
                        Graph& graph, size_t ncpus) :
      graph(graph),
      callback(engine),
      cpu_index(ncpus), cpu_color(ncpus), cpu_waiting(ncpus),
      max_iterations(0),
      update_function(NULL) {
      color.value = 0;
      // Verify the coloring
      if (graph.valid_coloring() == false) {
        graph.compute_coloring();
      }
      // Initialize the chromatic blocks
      for(vertex_id_t i = 0; i < graph.num_vertices(); ++i) {
        graphlab::vertex_color_type color = graph.color(i);
        if( color >= color_blocks.size() ) color_blocks.resize(color + 1);
        color_blocks[color].push_back(i);        
      }
    }

    
    /** Called by engine before executing the schedule */
    void start() {
      if(update_function == NULL) {
        std::cout << "No update function provided!" << std::endl;
      }
      assert(update_function != NULL);
      color.value = 0;
      // Initialize the cpu indexs
      for(size_t i = 0; i < cpu_index.size(); ++i) {
        cpu_index[i] = i;
        cpu_color[i] = -1;
        cpu_waiting[i] = true;
      }
      // Set waiting to zero
      waiting.value = 0;
      // REset the control termination
      terminator.reset();
    }


    
    /// Adds an update task with a particular priority
    void add_task(update_task_type task, double priority) {
      update_function = task.function();
    }
    
    /** 
     * Creates a collection of tasks on all the vertices in
     * 'vertices', and all with the same update function and priority
     */
    void add_tasks(const std::vector<vertex_id_t>& vertices, 
                   update_function_type func, double priority) {
      update_function = func;
    }
    
    /** 
     * Creates a collection of tasks on all the vertices in the graph,
     * with the same update function and priority
     */
    void add_task_to_all(update_function_type func, 
                         double priority) {
      update_function = func;
    }
    
    /**
     * This function returns a reference to the scheduling callback to
     * be used for a particular cpu. This callback will be passed to
     * update functions, and is the main interface which allow the
     * update functions to create new tasks.
     */
    callback_type& get_callback(size_t cpuid) {
      return callback;
    }

    /**
     * This function is called by the engine to ask for new work to
     * do.  The update task to be executed is returned in ret_task.
     *
     *  \retval NEWTASK There is an update task in ret_task to be
     *   executed
     *  \retval EMPTY Scheduler is empty
     */
    sched_status::status_enum get_next_task(size_t cpuid, 
                                            update_task_type &ret_task) {
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

      size_t current_color = cpu_color[cpuid].value % color_blocks.size();

      // Check to see that we have not run the maximum number of iterations
      if(max_iterations > 0 && 
         cpu_color[cpuid].value / color_blocks.size() >= max_iterations) {
        terminator.complete();
        return sched_status::EMPTY;
      }
      
      // If the index is good then execute it
      if(cpu_index[cpuid].value < color_blocks[ current_color ].size()) {
        vertex_id_t vertex =
          color_blocks[ current_color ][ cpu_index[cpuid].value ];
        ret_task = update_task_type(vertex, update_function);
        return sched_status::NEWTASK;
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
    void completed_task(size_t cpuid, 
                        const update_task_type &task) {} 

   

    terminator_type& get_terminator() {
      return terminator;
    };


    void set_options(const scheduler_options &opts) {
      opts.get_int_option("max_iterations", max_iterations);
      any uf;
      if (opts.get_any_option("update_function", uf)) {
        update_function = uf.as<update_function_type>();
      }
    }

    static void print_options_help(std::ostream &out) {
      out << "max_iterations = [integer, default = 0]\n";
      out << "update_function = [update_function_type,"
        "default = set on add_task]\n";
    };
  private:
    Graph& graph;
    
    /// The callbacks pre-created for each cpuid
    unused_scheduler_callback<Graph> callback;

    
    std::vector< std::vector< vertex_id_t> > color_blocks;
    std::vector< cache_line_pad<size_t> > cpu_index;
    std::vector< cache_line_pad<size_t> > cpu_color;
    std::vector< cache_line_pad<size_t> > cpu_waiting;

    size_t max_iterations;

    update_function_type update_function;

    atomic<size_t> color;
    atomic<size_t> waiting;

    terminator_type terminator;

   
    
    

 


  }; // End of chromatic scheduler

}
#endif

