/**  
 *  Software submitted by 
 *  Systems & Technology Research / Vision Systems Inc., 2013
 *
 *  Approved for public release; distribution is unlimited. [DISTAR Case #21428]
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


#ifndef TABLE_BASE_HPP
#define TABLE_BASE_HPP

/**
 * This file defines the root of the table hierarchy for
 * the various table types.
 *
 *  \author Scott Richardson     10/2012
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <cassert>
#include <iostream>

#include <graphlab/serialization/serialization_includes.hpp>

#include "discrete_variable.hpp"


namespace graphlab {


template<size_t MAX_DIM>
class table_base {
public:
  typedef table_base const *const const_ptr;

  virtual ~table_base() { }

  static inline double APPROX_LOG_ZERO() { 
    //return -std::numeric_limits<double>::max();
    return -1e6; 
  }

  virtual table_base& deep_copy(const table_base& base) = 0;
  virtual table_base& copy_onto(const table_base& base) = 0;

  //virtual table_base&   plus_equals(const table_base& base) = 0;
  virtual table_base&  times_equals(const table_base& base) = 0;
  virtual table_base& divide_equals(const table_base& base) = 0;
  //virtual void    plus(const table_base& base, table_base& out) const = 0;
  virtual void   times(const table_base& base, table_base& out) const = 0;
  virtual void  divide(const table_base& base, table_base& out) const = 0;

  virtual const discrete_variable& var(size_t index) const = 0;
  
  virtual void     MAP(table_base& msg) const = 0;
  virtual void    zero() = 0;
  virtual size_t numel() const = 0; // REVIEW might not be necessary
  virtual size_t ndims() const = 0;
  virtual void    load(graphlab::iarchive& arc) = 0;
  virtual void    save(graphlab::oarchive& arc) const = 0;

  friend std::ostream& operator<<(std::ostream& out,
                            const table_base<MAX_DIM>& factor) {
    factor.print(out);
    return out;
  }

private:
  virtual std::ostream& print(std::ostream& out = std::cout) const = 0;
}; // end of table_base

} // end of namespace graphlab 

#endif // TABLE_BASE_HPP
