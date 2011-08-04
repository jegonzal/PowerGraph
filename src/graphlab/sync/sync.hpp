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


#ifndef GRAPHLAB_SYNC_HPP
#define	GRAPHLAB_SYNC_HPP



#include <graphlab/logger/logger.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/sync/isync.hpp>



namespace graphlab {



  namespace sync_defaults {
    template<typename T, typename Accum>
    void apply(T& lvalue, const Accum& accum) {
      lvalue = accum;
    } // end of default_apply
  }; // end of sync_defaults namespace

  template<typename Graph, typename T, typename Accum >
  class fold_sync : public isync<Graph> {
  public:


    typedef isync<Graph> isync_type;
    typedef typename isync_type::iscope_type iscope_type;

    //! The target should supply = operation
    typedef glshared<T> glshared_type;
    typedef T      contained_type;
    typedef Accum  accumulator_type;

    

    /** 
     * The map function takes a scope and extracts the relevant
     * information into the result object.
     */
    typedef void(*fold_function_type)(iscope_type& scope, Accum& result);

    /**
     *  The apply function manipulates the partial sum and assigns it
     *  to the target gl shared object
     */
    typedef void(*apply_function_type)(T& lvalue, const Accum& accum);


    // /** 
    //  * The reduce function combines the right partial sum into the
    //  * left partial sum and behaves like:
    //  *
    //  *      partial_sum += rvalue
    //  *
    //  */
    // typedef void(*combine_function_type)(Accum& partial_sum, 
    //                                      const Accum& rvalue);







  private:
    glshared_type& target;
    accumulator_type zero;
    accumulator_type acc;

    fold_function_type fold_function;
    apply_function_type apply_function;

  public:

    fold_sync(const fold_sync& other) : 
      target(other.target), zero(other.zero), acc(other.acc) { }

    fold_sync(glshared_type& target,
              fold_function_type fold_function,
              const accumulator_type& zero = accumulator_type(0),
              apply_function_type apply_function = (sync_defaults::apply<T, Accum>) ) :
      target(target), zero(zero), acc(zero) { }


    
    isync_type* clone() const { return new fold_sync(*this); }
    void clear() { acc = zero; }
    void operator+=(iscope_type& scope) { fold_function(scope, acc); }
    void operator+=(const isync_type& iother) {
      const fold_sync& other = 
        *dynamic_cast<const fold_sync*>(&iother);
      acc += other.acc;
    }
    void apply() { 
      target.apply(apply_function, acc); 
    }
  }; // end of fold_sync
    

  
}; // end of Namespace graphlab


#endif
