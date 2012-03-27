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

    * This is the simplest of all schedulers.  It does absolutely
    * nothing but it requires very little memory to do nothing and
    * this makes some of our users happy. :-)
    */
  template<typename Graph, typename UpdateFunctor>
  class empty_scheduler : public ischeduler<Graph, UpdateFunctor> {
  public:


    typedef ischeduler<Graph, UpdateFunctor> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;
    typedef critical_termination terminator_type;


    
  private:
    terminator_type                         term;   
  public:
    empty_scheduler(const graph_type& graph, 
                    size_t ncpus,
                    const options_map& opts) : term(ncpus) { }

    void start() { }

    void schedule(const vertex_id_type vid, 
                  const update_functor_type& fun) { }

    void schedule_all(const update_functor_type& fun,
                      const std::string& order) { }
      
    
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {         
      return sched_status::EMPTY;
    } // end of get_next
    
    
    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const update_functor_type& fun) {  } // end of completed


    iterminator& terminator() { return term; };

    size_t num_joins() const { return 0; }

    static void print_options_help(std::ostream &out) {
      out << "This scheduler doesn't do anything." << std::endl;
    } // end of print_options_help
  }; 
  
  
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

