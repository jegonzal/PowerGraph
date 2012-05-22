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

#ifndef GRAPHLAB_SWEEP_SCHEDULER_HPP
#define GRAPHLAB_SWEEP_SCHEDULER_HPP

#include <queue>
#include <cmath>
#include <cassert>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_map.hpp>
#include <graphlab/scheduler/terminator/critical_termination.hpp>
#include <graphlab/options/options_map.hpp>
#include <graphlab/graph/graph_ops.hpp>
#include <graphlab/graph/graph_basic_types.hpp>



#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** \ingroup group_schedulers
    */
  template<typename Message>
  class sweep_scheduler : public ischeduler<Message> {
  public:
    typedef Message message_type;
    typedef critical_termination terminator_type;
    
  private:

    const size_t ncpus;

    bool strict_round_robin;
    atomic<size_t> rr_index;
    size_t max_iterations;
   

    std::vector<vertex_id_type>             vids;
    std::vector<uint16_t>                   vid2cpu;
    std::vector<vertex_id_type>             cpu2index;

    vertex_map<message_type>                messages;
    double                                  min_priority;
    terminator_type                         term;   


  public:
    sweep_scheduler(size_t num_vertices,
                    size_t ncpus,
                    const options_map& opts) :
      ncpus(ncpus),
      strict_round_robin(false),
      max_iterations(std::numeric_limits<size_t>::max()),
      vids(num_vertices),
      messages(num_vertices), 
      min_priority(-std::numeric_limits<double>::max()),
      term(ncpus) {
      
      // Determin the orering of updates
      std::string ordering = "random";
      opts.get_option("ordering", ordering) || 
        opts.get_option("order", ordering);
      if (ordering == "ascending") {
        logstream(LOG_INFO) 
          << "Using an ascending ordering of the vertices." << std::endl;
        for(size_t i = 0; i < num_vertices; ++i) vids[i] = i;
      }  else { // Assume random ordering by default
        if(ordering != "random") {
          logstream(LOG_WARNING)
            << "The ordering \"" << ordering << "\" is not supported using default."
            << std::endl;
        }
        logstream(LOG_INFO) 
          << "Using a random ordering of the vertices." << std::endl;
        for(size_t i = 0; i < num_vertices; ++i) vids[i] = i;
        random::shuffle(vids);
      }


      // Determine whether a strict ordering is to be used
      opts.get_option("strict", strict_round_robin);
      if(strict_round_robin) {
        logstream(LOG_INFO) 
          << "Using a strict round robin schedule." << std::endl;
        // Max iterations only applies to strict round robin
        if(opts.get_option("niters", max_iterations) ) 
          logstream(LOG_INFO) 
            << "Using maximum iterations: " << max_iterations << std::endl;
        // Initialize the round robin index
        rr_index = 0;
      } else { // each cpu is responsible for its own subset of vertices
        // Initialize the cpu2index counters
        cpu2index.resize(ncpus);
        for(size_t i = 0; i < cpu2index.size(); ++i) cpu2index[i] = i;
        // Initialze the reverse map vid2cpu assignment
        vid2cpu.resize(vids.size());
        for(size_t i = 0; i < vids.size(); ++i) vid2cpu[vids[i]] = i % ncpus;
      }

      // Get Min priority
      const bool is_set = opts.get_option("min_priority", min_priority);
      if(is_set) {
        logstream(LOG_INFO) 
          << "The minimum scheduling priority was set to " 
          << min_priority << std::endl;
      }
    } // end of constructor
        
   
    void start() { term.reset(); }

    void schedule(const vertex_id_type vid, 
                  const message_type& msg) {      
      double ret_priority = 0;
      if(messages.add(vid, msg, ret_priority) && ret_priority >= min_priority) {
        if(!vid2cpu.empty()) term.new_job(vid2cpu[vid]);
        else term.new_job();
      } 
    } // end of schedule

    void schedule_all(const message_type& msg,
                      const std::string& order) {
      for (vertex_id_type vid = 0; vid < messages.size(); ++vid)
        schedule(vid, msg);      
    } // end of schedule_all    
      
    
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       message_type& ret_msg) {         
      const size_t nverts    = vids.size();
      const size_t max_fails = (nverts/ncpus) + 1;
      // Check to see if max iterations have been achieved 
      if(strict_round_robin && (rr_index / nverts) >= max_iterations) 
        return sched_status::EMPTY;
      // Loop through all vertices that are associated with this
      // processor searching for a vertex with an active task
      for(size_t idx = get_and_inc_index(cpuid), fails = 0; 
          fails <= max_fails; // 
          idx = get_and_inc_index(cpuid), ++fails) {
        // It is possible that the get_and_inc_index could return an
        // invalid index if the number of cpus exceeds the number of
        // vertices.  In This case we alwasy return empty
        if(__builtin_expect(idx >= nverts, false)) return sched_status::EMPTY;
        const vertex_id_type vid = vids[idx];
        bool success = messages.test_and_get(vid, ret_msg);
        while(success) { // Job found now decide whether to keep it
          if(ret_msg.priority() >= min_priority) {
            ret_vid = vid; return sched_status::NEW_TASK;
          } else {
            // Priority is insufficient so return to the schedule
            double ret_priority = 0;
            messages.add(vid, ret_msg, ret_priority);
            // when the job was added back it could boost the
            // priority.  If the priority is sufficiently high we have
            // to try and remove it again. Now it is possible that if
            // strict ordering is used it could be taken again so we
            // may need to repeat the process.
            if(ret_priority >= min_priority) 
              success = messages.test_and_get(vid, ret_msg);
            else success = false;
          } 
        }// end of while loop over success
      } // end of for loop
      return sched_status::EMPTY;
    } // end of get_next
    
    
    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const message_type& msg) {
      term.completed_job();
    } // end of completed


    iterminator& terminator() { return term; };

    size_t num_joins() const {
      return messages.num_joins();
    }


    static void print_options_help(std::ostream &out) {
      out << "ordering = [string: {random, ascending, max_degree, "
          << "min_degree,color},"
          << "\t vertex ordering, default=random]\n"
          << "strict = [bool, use strict round robin schedule, default=false]\n"
          << "min_priority = [double, minimum priority required to receive \n"
          << "\t a message, default = -infinity]\n"
          << "niters = [integer, maximum number of iterations "
          << " (requires strict=true) \n"
          << "\t default = -infinity]\n";
    } // end of print_options_help


  private:
    inline size_t get_and_inc_index(const size_t cpuid) {
      const size_t nverts = vids.size();
      if (strict_round_robin) { 
        return rr_index++ % nverts; 
      } else {
        const size_t index = cpu2index[cpuid];
        cpu2index[cpuid] += ncpus;
        // Address loop around
        if (__builtin_expect(cpu2index[cpuid] >= nverts, false)) 
          cpu2index[cpuid] = cpuid;
        return index;
      }
    }// end of next index

  }; 
  
  
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

