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


#ifndef GRAPHLAB_MAP_REDUCE_AGGREGATOR_HPP
#define	GRAPHLAB_MAP_REDUCE_AGGREGATOR_HPP



#include <graphlab/logger/logger.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/shared_data/glshared.hpp>
#include <graphlab/aggregation/iaggregate.hpp>



namespace graphlab {




  template<typename Graph, typename Accum>
  class map_reduce_aggregator : public iaggregator<Graph> {
  public:
    typedef Accum accumulator_type;
    /** 
     * The map function takes a scope and extracts the relevant
     * information into the result object.
     */
    typedef void(*map_function_type)(iscope_type& scope, Accum& result);

    /** 
     * The reduce function combines the right partial sum into the
     * left partial sum and behaves like:
     *
     *      partial_sum += rvalue
     *
     */
    typedef bool(*reduce_function_type)(Accum& partial_sum, 
                                        const Accum& rvalue);


    /**
     *  The apply function manipulates the partial sum and assigns it
     *  to the target gl shared object
     */
    typedef void(*apply_function_type)(const Accum& accum, graphlab::any& target);

  private:
    accumulator_type partial_sum;
    accumulator_type temporary;
    bool cleared;
    map_function_type map_function;
    reduce_function_type reduce_function;

  public:

    map_reduce_aggregator(map_function_type map_function,
                   reduce_function_type reduce_function,
                   apply_function_type apply_function) :
      cleared(true), 
      map_function(map_function),
      reduce_function(reduce_function),
      apply_function(apply_function) { }

    void clear() { cleared = true; }

    void add(iscope_type& scope) {
      if(cleared) {
        map_function(scope, partial_sum);
      } else {
        map_function(scope, temporary);
        const bool success = reduce_function(partial_sum, temporary);
        ASSERT_TRUE(success);
      }
    }
    
    void add(const iaggregate* other_ptr) {
      ASSERT_TRUE(other_ptr != NULL);
      const map_reduce_aggregator* agg_task = 
        dynamic_cast<const map_reduce_aggregator*>(other_ptr);
      if(!agg_task->cleared) { 
        const bool success = 
          reduce_function(partial_sum, agg_task->partial_sum); 
        return success;
      } else {
        partial_sum = agg_task->partial_sum;
      }
    }

    bool add_delta(const void* delta) {
      ASSERT_TRUE(accumulator_ptr != NULL);
      const accumulator_type* delta = dynamic_cast<const accumulator_type*>(delta);
      if(cleared) {
        partial_sum = *delta;
      } else {
        const bool success = 
          reduce_function(partial_sum, *delta);
        return success;
      }
    }

    void apply(any& base) {

    }

    

  }; // end of map_reduce_aggregator
  
}; // end of Namespace graphlab


#endif
