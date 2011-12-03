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
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>
#include <graphlab/scheduler/terminator/critical_termination.hpp>
#include <graphlab/options/options_map.hpp>
#include <graphlab/graph/graph_ops.hpp>



#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** \ingroup group_schedulers
    */
  template<typename Engine>
  class sweep_scheduler : public ischeduler<Engine> {
  public:


    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;
    typedef critical_termination terminator_type;


    
  private:

    inline size_t get_and_inc_index(const size_t cpuid) {
      const size_t nverts = index2vid.size();
      if (strict_round_robin) { 
        return rr_index++ % nverts; 
      } else {
        const size_t ncpus = cpu2index.size();
        const size_t index = cpu2index[cpuid];
        cpu2index[cpuid] += ncpus;
        // Address loop around
        if (__builtin_expect(cpu2index[cpuid] >= nverts, false)) 
          cpu2index[cpuid] = cpuid;
        return index;
      }
    }// end of next index

    bool strict_round_robin;
    atomic<size_t> rr_index;
    size_t max_iterations;

    std::vector<vertex_id_type>             index2vid;
    std::vector<vertex_id_type>             cpu2index;
    vertex_functor_set<update_functor_type> vfun_set;
    double                                  min_priority;
    terminator_type                         term;
   


  public:
    sweep_scheduler(const graph_type& graph, 
                    size_t ncpus,
                    const options_map& opts) :
      strict_round_robin(false), rr_index(0), 
      max_iterations(std::numeric_limits<size_t>::max()),
      index2vid(graph.num_vertices()), cpu2index(ncpus),
      vfun_set(graph.num_vertices()), 
      min_priority(-std::numeric_limits<double>::max()),
      term(ncpus) {
      // Determine whether a strict ordering is to be used
      opts.get_option("strict", strict_round_robin);
      if(strict_round_robin) {
        logstream(LOG_INFO) 
          << "Using a strict round robin schedule." << std::endl;
        // Max iterations only applies to strict round robin
        if(opts.get_option("niters", max_iterations) ) {
          logstream(LOG_INFO) 
            << "Using maximum iterations: " << max_iterations << std::endl;
        }
      } 
      
      // Determin the orering of updates
      std::string ordering = "random";
      opts.get_option("ordering", ordering) || 
        opts.get_option("order", ordering);
      if (ordering == "ascending") {
        logstream(LOG_INFO) 
          << "Using an ascending ordering of the vertices." << std::endl;
        for(size_t i = 0; i < graph.num_vertices(); ++i) index2vid[i] = i;
      } else if (ordering == "max_degree" || ordering == "min_degree") {
        logstream(LOG_INFO) 
          << "Constructing a " << ordering << " sweep ordering." << std::endl;
        const bool use_max_degree = (ordering == "max_degree");
        typedef std::pair<int, vertex_id_type> pair_type;
        std::vector< pair_type > vec(graph.num_vertices());      
        for(vertex_id_type i = 0; i < vec.size(); ++i) {
          typedef graph_ops<graph_type> graph_ops;
          const int degree = graph_ops::num_neighbors(graph, i);
          vec[i] = use_max_degree? pair_type(-degree,i) : pair_type(degree,i);
        }
        std::sort(vec.begin(), vec.end());
        for(size_t i = 0; i < vec.size(); ++i) index2vid[i] = vec[i].second;
      } else if (ordering == "color") {
        logstream(LOG_INFO) 
          << "Constructing a color based sweep ordering." << std::endl;
        typedef std::pair<vertex_color_type, vertex_id_type> pair_type;        
        std::vector< pair_type > vec(graph.num_vertices());      
        for(vertex_id_type i = 0; i < vec.size(); ++i) 
          vec[i] = pair_type(graph.color(i), i);
        std::sort(vec.begin(), vec.end());
        for(size_t i = 0; i < vec.size(); ++i) index2vid[i] = vec[i].second;
      } else { // Assume random ordering by default
        if(ordering != "random") {
          logstream(LOG_WARNING)
            << "The ordering \"" << ordering << "\" is not supported using default."
            << std::endl;
        }
        logstream(LOG_INFO) 
          << "Using a random ordering of the vertices." << std::endl;
        for(size_t i = 0; i < graph.num_vertices(); ++i) index2vid[i] = i;
        random::shuffle(index2vid);
      }
      // Initialize the cpu2index counters
      for(size_t i = 0; i < cpu2index.size(); ++i)  cpu2index[i] = i;

      // Get Min priority
      const bool is_set = opts.get_option("min_priority", min_priority);
      if(is_set) {
        logstream(LOG_INFO) 
          << "The minimum scheduling priority was set to " 
          << min_priority << std::endl;
      }
    }
        
   
    void start() { term.reset(); }

    void schedule(const size_t cpuid,
                  const vertex_id_type vid, 
                  const update_functor_type& fun) {      
      if(vfun_set.add(vid, fun)) term.new_job();
    } // end of schedule

    void schedule_all(const update_functor_type& fun) {
      for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid)
        schedule(0, vid, fun);      
    } // end of schedule_all    
      
    
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {         
      const size_t nverts    = index2vid.size();
      const size_t ncpus     = cpu2index.size();
      const size_t max_fails = (nverts/ncpus) + 1;
      // Check to see if max iterations have been achieved 
      if(strict_round_robin && (rr_index / nverts) >= max_iterations) 
        return sched_status::EMPTY;
      // Loop through all vertices that are associated with this
      // processor searching for a vertex with an active task
      for(size_t idx = get_and_inc_index(cpuid), fails = 0; 
          fails <= max_fails; // 
          idx = get_and_inc_index(cpuid), ++fails) {
        ASSERT_LT(idx, nverts);
        const vertex_id_type vid = index2vid[idx];
        const bool success = vfun_set.test_and_get(vid, ret_fun);
        if(success) { // Job found now decide whether to keep it
          if(ret_fun.priority() >= min_priority) {
            ret_vid = vid; return sched_status::NEW_TASK;
          } else {
            // Priority is insufficient so return to the schedule
            schedule(cpuid, vid, ret_fun);
          } 
        }
      } // end of for loop
      return sched_status::EMPTY;
    } // end of get_next
    
    
    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const update_functor_type& fun) {
      term.completed_job();
    } // end of completed


    iterminator& terminator() { return term; };


    static void print_options_help(std::ostream &out) {
      out << "ordering = [string: {random, ascending, max_degree, min_degree, color}, "
          << "\t vertex ordering, default=random]\n"
          << "strict = [bool, use strict round robin schedule, default=false]\n"
          << "min_priority = [double, minimum priority required to run an \n"
          << "\t update functor, default = -infinity]\n"
          << "niters = [integer, maximum number of iterations (requires strict=true) \n"
          << "\t default = -infinity]\n";

    } // end of print_options_help

  }; 
  
  
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

