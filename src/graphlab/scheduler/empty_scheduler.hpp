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

#ifndef GRAPHLAB_EMPTY_SCHEDULER_HPP
#define GRAPHLAB_EMPTY_SCHEDULER_HPP

#include <queue>
#include <cmath>
#include <cassert>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/parallel/cache_line_pad.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/options/options_map.hpp>
#include <graphlab/graph/graph_ops.hpp>



#include <graphlab/macros_def.hpp>

namespace graphlab {

   /** \ingroup group_schedulers 

    * This is the simplest of all schedulers.  It does absolutely
    * nothing but it requires very little memory to do nothing and
    * this makes some of our users happy. :-)
    */
  template<typename Message>
  class empty_scheduler : public ischeduler<Message> {
  public:

    typedef Message message_type;
  public:
    empty_scheduler(size_t num_vertices,
                    size_t ncpus,
                    const options_map& opts) { }

    void start() { }

    void schedule(const vertex_id_type vid, 
                  const message_type& msg) { }

    void schedule_all(const message_type& msg,
                      const std::string& order) { }
      
    
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       message_type& ret_msg) {         
      return sched_status::EMPTY;
    } // end of get_next
    
    
    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const message_type& msg) {  } // end of completed


    size_t num_joins() const { return 0; }

    static void print_options_help(std::ostream &out) {
      out << "This scheduler doesn't do anything." << std::endl;
    } // end of print_options_help
  }; 
  
  
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

