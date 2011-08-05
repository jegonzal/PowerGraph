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
 * This class defines a very simple scheduler that loops vertices that
 * are "dirty". Each cpu loops vertices of which id%num_cpus==cpuid. Also called
 * "partitioned" scheduler.
 **/

#ifndef GRAPHLAB_SWEEP_SCHEDULER_HPP
#define GRAPHLAB_SWEEP_SCHEDULER_HPP

#include <queue>
#include <cmath>
#include <cassert>

#include <graphlab/parallel/cache_line_pad.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/engine/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>
#include <graphlab/engine/terminator/shared_termination.hpp>
#include <graphlab/options/options_map.hpp>


// #include <graphlab/util/synchronized_queue.hpp>
// #include <graphlab/tasks/update_task.hpp>
// #include <graphlab/schedulers/ischeduler.hpp>
// #include <graphlab/parallel/pthread_tools.hpp>
// #include <graphlab/schedulers/support/direct_callback.hpp>
// #include <graphlab/schedulers/support/vertex_task_set.hpp>
// #include <graphlab/util/shared_termination.hpp>

// #include <graphlab/parallel/atomic.hpp>




#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** \ingroup group_schedulers
    */
  template<typename Engine>
  class sweep_scheduler : 
    public ischeduler<Engine> {
  public:


    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;
    typedef shared_termination terminator_type;

    
  private:

    inline size_t next_index(const size_t ncpus, 
                             const size_t nverts,
                             const size_t cpuid, 
                             const size_t index) const {
      const size_t nextidx = index + ncpus;
      if(__builtin_expect(nextidx < nverts, true)) 
        return nextidx; 
      else return cpuid;
    }// end of next index


    std::vector<vertex_id_type>     index2vid;
    std::vector<uint8_t>             vid2cpuid;
    std::vector< cache_line_pad<size_t> >   cpu2index;
    vertex_functor_set<engine_type> vfun_set;
    terminator_type                 term;



  public:
    sweep_scheduler(const graph_type& graph, 
                    size_t ncpus,
                    const options_map& opts) : 
      index2vid(graph.num_vertices()), 
      vid2cpuid(graph.num_vertices()),
      cpu2index(ncpus),
      vfun_set(graph.num_vertices()), term(ncpus) {
      // Construct the permutation
      for(size_t i = 0; i < graph.num_vertices(); ++i) index2vid[i] = i;
      std::string ordering = "shuffle";
      opts.get_string_option("ordering", ordering);
      if(ordering == "shuffle") random::shuffle(index2vid);     
      // construct the vid2cpuid map
      for(size_t i = 0; i < graph.num_vertices(); ++i) 
        vid2cpuid[index2vid[i]] = i % ncpus;
    }
        
    void start() { 
      for(size_t i = 0; i < cpu2index.size(); ++i) 
        cpu2index[i] = i;
      term.reset();
    }


    void schedule(vertex_id_type vid, 
                  const update_functor_type& fun) {      
      if(vfun_set.add(vid, fun)) term.new_job(vid2cpuid[vid]);       
    } // end of schedule

    void schedule_all(const update_functor_type& fun) {
      for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid)
        schedule(vid, fun);      
    } // end of schedule_all    
      

    
    /** Get next dirty vertex in queue. Each cpu checks vertices with
        modulo num_cpus = cpuid */
    sched_status::status_enum get_next(size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {         
      const size_t start_index = cpu2index[cpuid];
      // Loop through all vertices that are associated with this
      // processor searching for a vertex with an active task
      const size_t ncpus = cpu2index.size();
      const size_t nverts = index2vid.size();

      for(size_t i = next_index(ncpus,nverts,cpuid,start_index); 
          i != start_index; 
          i = next_index(ncpus,nverts,cpuid,i) ) {
        ASSERT_LT(i, nverts);
        const vertex_id_type vid = index2vid[i];
        const bool success = vfun_set.test_and_get(vid, ret_fun);
        if(success) {
          ret_vid = vid;
          cpu2index[cpuid] = i;
          return sched_status::NEW_TASK;
        }
      }
      return sched_status::EMPTY;
    } // end of get_next
    
    
    void completed(size_t cpuid,
                   vertex_id_type vid,
                   const update_functor_type& fun) {
      term.completed_job();
    } // end of completed


    iterminator& terminator() { return term; };


    static void print_options_help(std::ostream &out) {
      out << "ordering = [string: shuffle/linear, default=shuffle]\n";
    };

  }; 
  
  
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

