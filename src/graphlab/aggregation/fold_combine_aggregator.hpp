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


#ifndef GRAPHLAB_FOLD_COMBINE_AGGREGATOR_HPP
#define	GRAPHLAB_FOLD_COMBINE_AGGREGATOR_HPP



#include <graphlab/logger/logger.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/aggregation/iaggregator.hpp>



namespace graphlab {


  template<typename Graph, typename Accum>
  void default_fold(const iscope<Graph>& scope, Accum& result) {
    result += scope.const_vertex_data();
  }
  
  template<typename Target, typename Accum>
  void default_apply(Target& target, Accum& accum) {
    target = accum;
  }

  template<typename Accum>
  void default_combine(Accum& acc, const Accum& rvalue) {
    acc += rvalue;
  }
  


  template<typename Graph, typename Target, typename Accum>
  class fold_combine_aggregator : public iaggregator<Graph> {
  public:

    typedef iaggregator<Graph> iaggregator_type;

    typedef typename iaggregator_type::iscope_type iscope_type;

    //! The target should supply = operation
    typedef Target target_type;
    typedef Accum  accumulator_type;
    typedef typename Target::contained_type contained_type;
    

    /** 
     * The map function takes a scope and extracts the relevant
     * information into the result object.
     */
    typedef void(*fold_function_type)(iscope_type& scope, Accum& result);

    /**
     *  The apply function manipulates the partial sum and assigns it
     *  to the target gl shared object
     */
    typedef void(*apply_function_type)(contained_type& target, Accum& accum);


    /** 
     * The reduce function combines the right partial sum into the
     * left partial sum and behaves like:
     *
     *      partial_sum += rvalue
     *
     */
    typedef void(*combine_function_type)(Accum& partial_sum, 
                                         const Accum& rvalue);







  private:
    target_type&     target;
    accumulator_type zero;
    accumulator_type acc;


    fold_function_type    fold_function;
    apply_function_type   apply_function;
    combine_function_type combine_function;

  public:

    fold_combine_aggregator(target_type& target,
                            const accumulator_type& zero,
                            fold_function_type fold_function,
                            apply_function_type apply_function,
                            combine_function_type combine_function) :
      target(target),
      zero(zero),
      acc(zero),
      fold_function(fold_function),
      apply_function(apply_function),
      combine_function(combine_function) { }

    
    iaggregator_type* clone() { return new fold_combine_aggregator(*this); }
    void clear() { acc = zero; }
    void operator+=(const iscope_type& scope) { fold_function(scope, acc); }
    void operator+=(const iaggregator_type& iother) {
      const fold_combine_aggregator& other = 
        *dynamic_cast<const fold_combine_aggregator*>(&iother);
      combine_function(acc, other.acc);
    }
    void apply() { target.apply(apply_function, acc); }
  }; // end of fold_combine_aggregator
  
}; // end of Namespace graphlab


#endif
