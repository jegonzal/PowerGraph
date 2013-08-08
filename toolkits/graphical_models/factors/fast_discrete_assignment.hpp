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


#ifndef FAST_DISCRETE_ASSIGNMENT_HPP
#define FAST_DISCRETE_ASSIGNMENT_HPP

#include "discrete_variable.hpp"
#include "discrete_domain.hpp"
#include "discrete_assignment.hpp"




//#include <graphlab/macros_def.hpp>
namespace graphlab {


  /**
   * A limited version of discrete_assignment which supports
   * a smaller (and different) set of operations
   */
  template<size_t MAX_DIM>
  class fast_discrete_assignment {
      
  public:

    fast_discrete_assignment(const discrete_domain<MAX_DIM>& args) :
        _num_vars(args.num_vars()), _index(0)
    { 
      transposed = false;
      for (size_t i = 0;i < _num_vars; ++i) {
        _vars[i] = args.var(i);
        _asgs[i] = 0;
      }
      size_t multiple = 1;
      for (size_t i = 0;i < _num_vars; ++i) {
        _increment_step[i] = multiple;
        multiple *= _vars[i].size();
      }
    }

    //! Construct a fast_discrete_assignment from a discrete_assignment
    fast_discrete_assignment(const discrete_assignment<MAX_DIM>& asg) : 
      _num_vars(asg.args().num_vars()), _index(asg.linear_index()) { 
      transposed = false;
      typename discrete_assignment<MAX_DIM>::const_iterator asg_it = asg.begin();
      for (size_t i = 0; i < _num_vars; ++i) {
        _vars[i] = asg.args().var(i);
        _asgs[i] = asg_it[i];
      }
        
      size_t multiple = 1;
      for (size_t i = 0;i < _num_vars; ++i) {
        _increment_step[i] = multiple;
        multiple *= _vars[i].size();
      }
    }
    
    //! get the number of variables
    size_t num_vars() const { return _num_vars; }

    //! get the number of variables
    discrete_variable var(size_t i) const { return _vars[i]; }


    size_t linear_index() const { return _index; }
      
    //! Get the next fast_discrete_assignment
    fast_discrete_assignment& operator++() {
      // Update the discrete_assignments
      for(size_t i = 0; i < num_vars(); ++i) {
        if (_asgs[i] < (_vars[i].size() - 1)) {
          _asgs[i] = (_asgs[i] + 1);
          _index += _increment_step[i];
          return *this;
        }
        else {
          _index -= _asgs[i] * _increment_step[i];
          _asgs[i] = 0;
        }
      }
      // Reached end
      make_end();
      return *this;
    }

    void set_index(size_t index) {
      _index = index;
      recompute_asgs();
    }

    size_t asg(size_t var_id) const {
      size_t idx = var_location(var_id);
      assert(idx < num_vars());
      return _asgs[idx];
    }

    void set_asg(size_t var_id, size_t value) {
      size_t idx = var_location(var_id);
      assert(idx < num_vars());
      assert(value < var(idx).size());
      _asgs[idx] = value;
      recompute_linear_index();
    }

    //! Tests whether two fast_discrete_assignments are equal
    bool operator==(const fast_discrete_assignment& other) const {
      return _index == other._index;
    }
    //! Tests whether two fast_discrete_assignments are not equal
    bool operator!=(const fast_discrete_assignment& other) const {
      return _index != other._index;
    }

    //! Make this an ending fast_discrete_assignment
    void make_end() {
      _index = -1;
    }

    /** Makes the sub_domain the first set of variables to be incremented over
     * Can only be called once
     */
    void transpose_to_start(const discrete_domain<MAX_DIM>& sub_domain) {
      ASSERT_FALSE(transposed);
      transposed = true;
        
      size_t reorder_map[MAX_DIM];
      size_t cursubdomain_idx = 0;
      size_t remainder_idx = sub_domain.num_vars();
      for (size_t i = 0;i < num_vars(); ++i) {
        if (cursubdomain_idx < sub_domain.num_vars() && 
            _vars[i].id() == sub_domain.var(cursubdomain_idx).id()) {
          reorder_map[cursubdomain_idx] = i;
          ++cursubdomain_idx;
        }
        else {
          reorder_map[remainder_idx] = i;
          ++remainder_idx;
        }
      }
      //move the asg around

        
      uint16_t newasgs[MAX_DIM];
      size_t newincrement_step[MAX_DIM]; 
      discrete_variable newvars[MAX_DIM];
      for (size_t i = 0;i < num_vars() ; ++i) {
        newincrement_step[i] = _increment_step[reorder_map[i]];
        newasgs[i] = _asgs[reorder_map[i]];
        newvars[i] = _vars[reorder_map[i]];
      }
      // copyback
      for (size_t i = 0;i < num_vars(); ++i) {
        _asgs[i] = newasgs[i];
        _vars[i] = newvars[i];
        _increment_step[i] = newincrement_step[i];
      }
    } // end of transpose_to_start



  private:
    //! Recompute the index from the discrete_assignment
    void recompute_linear_index() {
      size_t multiple = 1;
      // Clear the index
      _index = 0;
      for(size_t i = 0; i < num_vars(); ++i) {
        _index += multiple * _asgs[i];
        //        assert(_args.var(i).nasgs > 0);
        multiple *= _vars[i].size();
      }
    }
      
    //! Recompute the discrete_assignments from the index
    void recompute_asgs() {
      size_t quotient = _index;
      for(size_t i = 0; i < num_vars(); ++i) {
        _asgs[i] = quotient % _vars[i].size();
        quotient /= _vars[i].size();
        // assert(_asgs[i] < _args.var(i).size());
      }
    }

    /** get the index of the variable or returns number of variables
        if the index is not found */
    size_t var_location(size_t var_id) const {
      size_t location = num_vars();
      for(size_t i = 0; i < num_vars() && !(location < num_vars()); ++i) {
        if(_vars[i].id() == var_id) location = i;
      }
      return location;
    }


    size_t _num_vars;
    uint32_t _index;    
    discrete_variable _vars[MAX_DIM]; // actual ordering of the assignments
    size_t _increment_step[MAX_DIM];  //increment ordering according to _vars
    uint16_t _asgs[MAX_DIM];          // assignments with respect to _vars
    bool transposed;
  };
}; // end of namespace graphlab





template<size_t MAX_DIM>
std::ostream& operator<<(std::ostream& out,
                         const graphlab::fast_discrete_assignment<MAX_DIM>& asg) {
  out << "{";
  for(size_t i = 0; i < asg.num_vars(); ++i) {
    out << "v_" << asg.var(i).id();
    if(i < asg.num_vars() - 1) out << ", ";
  }
  out << "}=" << asg.linear_index();
  return out;
}


//#include <graphlab/macros_undef.hpp>
#endif // FAST_DISCRETE_ASSIGNMENT_HPP

