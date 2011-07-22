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


#ifndef GRAPHLAB_SCHEDULER_FACTORY_HPP
#define GRAPHLAB_SCHEDULER_FACTORY_HPP

#include <string>
#include <sstream>



// Schedulers
#include <graphlab/options/options_map.hpp>
#include <graphlab/scheduler/scheduler_list.hpp>
#include <graphlab/engine/terminator/iterminator.hpp>

#include <boost/preprocessor.hpp>


namespace graphlab {
  
  
  /**
   *  helper for constructing graphlab engines.
   **/
  template<typename Engine>
  struct scheduler_factory {
    
    typedef Engine engine_type;
    typedef typename engine_type::graph_type graph_type;
    typedef typename engine_type::ischeduler_type ischeduler_type;

   
    /**
     * Construct the a scheduler
     */
    template<typename Scheduler>
    ischeduler_type* new_scheduler(const graph_type& graph,
                                   iterminator& terminator,
                                   const size_t& ncpus,
                                   const options_map& opts) {
      ischeduler_type* scheduler_ptr = 
        new Scheduler(graph, terminator, ncpus, opts);
      ASSERT_TRUE(scheduler_ptr != NULL);
      scheduler_ptr->set_scheduler_options(opts);
      return scheduler_ptr;
    }
  
   


    /**
     * This function returns a new scheduler for a particular engine
     */    
    ischeduler<Engine>* new_scheduler(const std::string& scheduler_str,
                                      const options_map& opts,
                                      const graph_type& graph,
                                      iterminator& terminator,
                                      const size_t& ncpus) {
#define __GENERATE_NEW_SCHEDULER__(r_unused, data_unused, i,  elem)     \
      BOOST_PP_EXPR_IF(i, else)                                         \
        if (scheduler_str == BOOST_PP_TUPLE_ELEM(3,0,elem)) {           \
          return new_scheduler<BOOST_PP_TUPLE_ELEM(3,1,elem)>           \
            ( graph, terminator, ncpus, opts );                         \
        }      
      // generate the construction calls
      BOOST_PP_SEQ_FOR_EACH_I(__GENERATE_NEW_SCHEDULER__, _, __SCHEDULER_LIST__);
#undef __GENERATE_NEW_SCHEDULER__        
      std::cout 
        << "Invalid scheduler type: " << scheduler_str << std::endl;
      return NULL;
    } // end of new_scheduler

  }; // end of class scheduler_factory

}; // End of namespace graphlab




#endif

