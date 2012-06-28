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
#include <graphlab/parallel/atomic_add_vector2.hpp>
#include <graphlab/graph/graph_basic_types.hpp>

#include <graphlab/scheduler/get_message_priority.hpp>
#include <graphlab/options/graphlab_options.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** \ingroup group_schedulers
    */
  template<typename Message>
  class sweep_scheduler : public ischeduler<Message> {
  public:
    typedef Message message_type;
    
  private:

    size_t ncpus;

    bool strict_round_robin;
    atomic<size_t> rr_index;
    size_t max_iterations;
   

    std::vector<lvid_type>             vids;
    std::vector<uint16_t>                   vid2cpu;
    std::vector<lvid_type>             cpu2index;

    atomic_add_vector2<message_type>         messages;
    double                                  min_priority;
    std::string                             ordering;

  public:
    sweep_scheduler(size_t num_vertices,
                    const graphlab_options& opts) :
          ncpus(opts.get_ncpus()),
          strict_round_robin(false),
          max_iterations(std::numeric_limits<size_t>::max()),
          vids(num_vertices),
          messages(num_vertices),
          min_priority(-std::numeric_limits<double>::max()) {
        // initialize defaults
      ordering = "random";
      set_options(opts);

      for(size_t i = 0; i < vids.size(); ++i) vids[i] = i;
      if (ordering == "ascending") {
        logstream(LOG_INFO) << "Using an ascending ordering of the vertices." << std::endl;
      } else if(ordering == "random") {
        logstream(LOG_INFO)  << "Using a random ordering of the vertices." << std::endl;
        random::shuffle(vids);
      }

      if(strict_round_robin) {
        logstream(LOG_INFO)
          << "Using a strict round robin schedule." << std::endl;
        // Max iterations only applies to strict round robin
        if(max_iterations != std::numeric_limits<size_t>::max()) {
          logstream(LOG_INFO)
            << "Using maximum iterations: " << max_iterations << std::endl;
        }
        rr_index = 0;
      } else {
        // each cpu is responsible for its own subset of vertices
        // Initialize the cpu2index counters
        cpu2index.resize(ncpus);
        for(size_t i = 0; i < cpu2index.size(); ++i) cpu2index[i] = i;
        // Initialze the reverse map vid2cpu assignment
        vid2cpu.resize(vids.size());
        for(size_t i = 0; i < vids.size(); ++i) vid2cpu[vids[i]] = i % ncpus;
      }
        
    } // end of constructor
        
   

    void set_options(const graphlab_options& opts) {
      size_t new_ncpus = opts.get_ncpus();
      if (new_ncpus != ncpus) {
        logstream(LOG_INFO) << "Changing ncpus from " << ncpus << " to " << new_ncpus << std::endl;
        ASSERT_GE(new_ncpus, 1);
        ncpus = new_ncpus;
      }
      std::vector<std::string> keys = opts.get_scheduler_args().get_option_keys();
      foreach(std::string opt, keys) {
        if (opt == "order") {
          opts.get_scheduler_args().get_option("order", ordering);
          ASSERT_TRUE(ordering == "random" || ordering == "ascending");
        } else if (opt == "strict") {
          opts.get_scheduler_args().get_option("strict", strict_round_robin);
        } else if (opt == "max_iterations") {
          opts.get_scheduler_args().get_option("max_iterations", max_iterations);
        } else if (opt == "min_priority") {
          opts.get_scheduler_args().get_option("min_priority", min_priority);
        } else {
          logstream(LOG_FATAL) << "Unexpected Scheduler Option: " << opt << std::endl;
        }
      }
    }
   
    void start() { 
    }

    void schedule(const lvid_type vid, 
                  const message_type& msg) {      
      messages.add(vid, msg);
    } // end of schedule


    void schedule_from_execution_thread(const size_t cpuid,
                                        const lvid_type vid) {      
    } // end of schedule


    void schedule_all(const message_type& msg,
                      const std::string& order) {
      for (lvid_type vid = 0; vid < messages.size(); ++vid)
        schedule(vid, msg);      
    } // end of schedule_all    
      
    
    sched_status::status_enum 
    get_specific(lvid_type vid,
                 message_type& ret_msg) {
      bool get_success = messages.test_and_get(vid, ret_msg); 
      if (get_success) return sched_status::NEW_TASK;
      else return sched_status::EMPTY;
    }

    void place(lvid_type vid,
                 const message_type& msg) {
      messages.add(vid, msg);
    }
    
    sched_status::status_enum get_next(const size_t cpuid,
                                       lvid_type& ret_vid,
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
        const lvid_type vid = vids[idx];
        bool success = messages.test_and_get(vid, ret_msg);
        while(success) { // Job found now decide whether to keep it
          if(scheduler_impl::get_message_priority(ret_msg) >= min_priority) {
            ret_vid = vid; return sched_status::NEW_TASK;
          } else {
            // Priority is insufficient so return to the schedule
            message_type combined_message;
            messages.add(vid, ret_msg, combined_message);
            double ret_priority = scheduler_impl::get_message_priority(combined_message);
            // when the job was added back it could boost the
            // priority.  If the priority is sufficiently high we have
            // to try and remove it again. Now it is possible that if
            // strict ordering is used it could be taken again so we
            // may need to repeat the process.
            if(ret_priority >= min_priority) {
              success = messages.test_and_get(vid, ret_msg);
            } else {
              success = false;
            }
          } 
        }// end of while loop over success
      } // end of for loop
      return sched_status::EMPTY;
    } // end of get_next
    
    
    void completed(const size_t cpuid,
                   const lvid_type vid,
                   const message_type& msg) {
    } // end of completed


    size_t num_joins() const {
      return messages.num_joins();
    }



    static void print_options_help(std::ostream &out) {
      out << "order = [string: {random, ascending} default=random]\n"
          << "strict = [bool, use strict round robin schedule, default=false]\n"
          << "min_priority = [double, minimum priority required to receive \n"
          << "\t a message, default = -inf]\n"
          << "max_iterations = [integer, maximum number of iterations "
          << " (requires strict=true) \n"
          << "\t default = inf]\n";
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

