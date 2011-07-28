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

 
  namespace sync_impl {
    template<typename T, typename Accum>
    void default_apply(T& lvalue, const Accum& accum) {
      lvalue = accum;
    }; // end of default_apply
  }; // end of namespace sync_impl;


  template<typename Graph>
  struct sync {

    typedef isync<Graph> isync_type;
    typedef iscop<Graph> iscope_type;

    typename<typename T, typename Accum,
              void(*fold_function)(iscope_type& scope, Accum& result),
              void(*apply_function)(T& lvalue, const Accum& accum) = 
              sync_impl::default_apply<T, Accum>  >
    class fold : public isync_type {
    public:


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


    public:

      fold(glshared_type& target,
           const accumulator_type& zero = accumulator_type(0) ) :
        target(target),
        zero(zero),
        acc(zero) { }


    
      isync_type* clone() { return new fold_sync(*this); }
      void clear() { acc = zero; }
      void operator+=(iscope_type& scope) { fold_function(scope, acc); }
      void operator+=(const isync_type& iother) {
        const fold& other = 
          *dynamic_cast<const fold_sync*>(&iother);
        acc += other.acc;
      }
      void apply() { target.apply(apply_function, acc); }
    }; // end of fold_sync
    
  }; // end of struct sync  
  
}; // end of Namespace graphlab


#endif
