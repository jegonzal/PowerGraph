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


#ifndef GRAPHLAB_SHAREDSUM_HPP
#define GRAPHLAB_SHAREDSUM_HPP
// #include <boost/shared_ptr.hpp>
// #include <boost/function.hpp>
// #include <boost/type_traits/function_traits.hpp>
// #include <boost/type_traits/remove_reference.hpp>

//#include <graphlab/parallel/atomic.hpp>

#include <graphlab/shared_data/isharedsum.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/logger/assertions.hpp>



namespace graphlab {

  namespace sharedsum_impl {
    template<typename T, typename Acc> 
    struct identity {
      void operator()(T& ret, const Acc& acc) const { ret = acc; }
    };
    struct ipartial_sum { virtual ~ipartial_sum() { } };
    ipartial_sum* get_tl_partial_sum(void* sharedsum_ptr);
    void set_tl_partial_sum(void* sharedsum_ptr, ipartial_sum* ps_ptr);
  };

  template <typename T, typename Acc = T, 
            typename ApplyFunctor = sharedsum_impl::identity<T,Acc> >
  class sharedsum {

  public:
    typedef T   contained_type;
    typedef Acc accumulator_type;
    typedef ApplyFunctor apply_functor_type;
    typedef sharedsum_impl::ipartial_sum ipartial_sum;

    struct partial_sum : public ipartial_sum {
      lock lock;
      accumlator_type acc;
      partial_sum(const accumulator_type& acc = accumulator_type()) :
        acc(acc) { }
      void zero() { acc = accumulator_type(); }
    };

  private:
    rwlock lock;
    contained_type contents;

    lock ps_lock;
    std::vector<partial_sum*> partial_sums;

    partial_sum& get_tl_partial_sum() {
      ipartial_sum* ips_ptr = sharedsum_impl::get_tl_partial_sum(this);
      if(ips_ptr == NULL) {
        ips_ptr = new partial_sum();
        shared_sum_impl::set_tl_partial_sum(this, ips_ptr);
        ps_lock.lock();
        partial_sums.push_back(
      }
      ASSERT_NE(ips_ptr, NULL);
      partial_sum* ps_ptr = dynamic_cast<partial_sum*>(ips_ptr);
      return *ps_ptr;
    }


  public:

    sharedsum(const T& val = T()) : contents(val) { }
  

    //! Assign a new value 
    void operator=(const T& val) {
      lock.writelock();
      contents = val;
      lock.unlock();
    }

    //! Add a delta function:
    void operator+=(const Acc& acc) { 
      
    }


    /// Returns a copy of the data
    inline T get_val() const {
      lock.readlock();
      const T copy = contents;
      lock.unlock();
      return copy;
    }


  };


 
} 
#endif

