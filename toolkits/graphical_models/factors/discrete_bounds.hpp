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


#ifndef DISCRETE_BOUNDS_HPP
#define DISCRETE_BOUNDS_HPP

#include <graphlab/logger/assertions.hpp>

#include "discrete_variable.hpp"


#include <graphlab/macros_def.hpp>
namespace graphlab {


/**
 * This class respresents a discrete domain over a set of variables.
 */
template<size_t MAX_DIM>
class discrete_bounds {
public:
  //! Make an empty domain
  discrete_bounds() : _num_vars(0) { }

  //! Make a single variable discrete_bounds
  discrete_bounds(const discrete_variable& v1) :
      _num_vars(1) 
  {
    DCHECK_LE(_num_vars, MAX_DIM);
    _vars[0] = v1;
  }

  //! Make a two variable discrete_bounds
  discrete_bounds(const discrete_variable& v1, const discrete_variable& v2) :
      _num_vars(2) 
  {
    DCHECK_LE(_num_vars, MAX_DIM);
    assert(v1 != v2);
    if(v1 < v2) {
      _vars[0] = v1;
      _vars[1] = v2;
    } else {
      _vars[0] = v2;
      _vars[1] = v1;
    }
  }

  //! Make a three variable discrete_bounds
  discrete_bounds(const discrete_variable& v1,
                  const discrete_variable& v2,
                  const discrete_variable& v3) :
      _num_vars(3) 
  {
    DCHECK_LE(_num_vars, MAX_DIM);
    DCHECK_NE(v1, v2);
    DCHECK_NE(v2, v3);
    DCHECK_NE(v1, v3);
      
    if(v1 < v2 && v2 < v3) {
      _vars[0] = v1;
      _vars[1] = v2;
      _vars[2] = v3;
    } else if( v1 < v3 && v3 < v2) {
      _vars[0] = v1;
      _vars[1] = v3;
      _vars[2] = v2;
    } else if( v2 < v1 && v1 < v3) {
      _vars[0] = v2;
      _vars[1] = v1;
      _vars[2] = v3;
    } else if( v2 < v3 && v3 < v1) {
      _vars[0] = v2;
      _vars[1] = v3;
      _vars[2] = v1;
    } else if( v3 < v1 && v1 < v2) {
      _vars[0] = v3;
      _vars[1] = v1;
      _vars[2] = v2;
    } else if( v3 < v2 && v2 < v1) {
      _vars[0] = v3;
      _vars[1] = v2;
      _vars[2] = v1;
    } else { throw("Invalid Case!"); }
  }

  //! Make a discrete_bounds from a vector of variables
  explicit discrete_bounds(const std::vector<discrete_variable>& vars) :
      _num_vars(vars.size()) 
  {
    DCHECK_LE(_num_vars, MAX_DIM);
    for(size_t i = 0; i < _num_vars; ++i) {
      _vars[i] = vars[i];
    }
    std::sort(_vars, _vars + std::min(MAX_DIM, _num_vars) );
  }

  //! Make a discrete_bounds from a set of variables
  explicit discrete_bounds(const std::set<discrete_variable>& vars) :
      _num_vars(vars.size()) 
  {
    DCHECK_LE(_num_vars, MAX_DIM); 
    size_t i = 0; 
    foreach(const discrete_variable& var, vars) _vars[i++] = var;
  }
  
  discrete_bounds(const discrete_bounds& other) :
      _num_vars(other._num_vars) 
  {
    *this = other;
  }

  virtual ~discrete_bounds() { }

  discrete_bounds& operator=(const discrete_bounds& other) {
    if(this == &other) 
      return *this;
    
    _num_vars = other._num_vars;
    DCHECK_LE(_num_vars, MAX_DIM);
    for(size_t i = 0; i < _num_vars; ++i) {
      _vars[i] = other.var(i);
    }
    return *this;
  }

  //! test whether two discrete_boundss are equal
  bool operator==(const discrete_bounds& other) const {
    if( num_vars() != other.num_vars() ) return false;  
    for(size_t i = 0; i < num_vars(); ++i) {
      if(var(i) != other.var(i)) return false;
    }
    return true;
  }
    
  //! test whether two discrete_boundss are not equal
  bool operator!=(const discrete_bounds& other) const {
    return !(*this == other);
  }

  //! add the other discrete_bounds to this discrete_bounds
  discrete_bounds& operator+=(const discrete_variable& var) {
    if(_vars[_num_vars - 1] < var) {
      _vars[_num_vars] = var;
      _num_vars++;
      return *this;
    }
    return operator+=(discrete_bounds(var));
  }

  //! add the discrete_bounds to this discrete_bounds
  discrete_bounds& operator+=(const discrete_bounds& other) {
    if(other.num_vars() == 0) return *this;
    discrete_bounds backup = *this;
    _num_vars = 0;
    for(size_t i = 0, j = 0; 
        i < backup.num_vars() || j < other.num_vars(); ) {
      DCHECK_LE(_num_vars, MAX_DIM);
      // Both 
      if(i < backup.num_vars() && j < other.num_vars() 
         && _num_vars < MAX_DIM) {
        if(backup.var(i) < other.var(j))  
          _vars[_num_vars++] = backup.var(i++);
        else if(other.var(j) < backup.var(i))  
          _vars[_num_vars++] = other.var(j++);
        else { _vars[_num_vars++] = backup.var(i++); j++; }
      } else if(i < backup.num_vars() && _num_vars < MAX_DIM) {
        _vars[_num_vars++] = backup.var(i++);
      } else if(j < other.num_vars() && _num_vars < MAX_DIM) {
        _vars[_num_vars++] = other.var(j++);
      } else {
        *this = backup;
        // Unreachable
        throw("Unreachable case in domain operator+=");
      }
    }
    return *this;
  }
  
  //! add two discrete_boundss together
  discrete_bounds operator+(const discrete_variable& var) const {
    discrete_bounds dom = *this;
    return dom += var;
  }

  //! add the other discrete_bounds to this discrete_bounds
  discrete_bounds operator+(const discrete_bounds& other) const {
    discrete_bounds dom = *this;
    return dom += other;
  }

  
  //! subtract the other discrete_bounds from this discrete_bounds
  discrete_bounds& operator-=(const discrete_bounds& other) {
    if(other.num_vars() == 0) return *this;
      
    size_t tmp_num_vars = 0;
    for(size_t i = 0, j = 0; i < _num_vars; ++i ) {
      // advance the other index
      for( ; j < other._num_vars && _vars[i].id() > other._vars[j].id(); ++j) { }

      if(!(j < other._num_vars && _vars[i].id() == other._vars[j].id())) {
        _vars[tmp_num_vars++] = _vars[i];
      }
    }
    _num_vars = tmp_num_vars;
    return *this;
  }

  //! subtract the other discrete_bounds from this discrete_bounds
  discrete_bounds operator-(const discrete_bounds& other) const {
    discrete_bounds dom(*this);
    return dom -= other;
  }

  discrete_bounds intersect(const discrete_bounds& other) const {
    discrete_bounds new_dom;
    new_dom._num_vars = 0;
    for(size_t i = 0, j = 0;
        i < num_vars() && j < other.num_vars(); ) {
      if(_vars[i] == other._vars[j]) {
        // new discrete_bounds gets the variable
        new_dom._vars[new_dom.num_vars()] = _vars[i];
        // Everyone advances
        new_dom._num_vars++;  i++; j++;
      } else {
        // otherwise increment one of the variables          
        if(_vars[i] < other._vars[j]) i++; else j++;
      }
    }
    return new_dom;
  }

  //! Get the number of variables
  inline size_t num_vars() const { return _num_vars; }

  //! Get the ith variable
  inline const discrete_variable& var(size_t index) const {
    DCHECK_LT(index, _num_vars);
    return _vars[index];
  }

  /** get the index of the variable or returns number of variables
      if the index is not found */
  size_t var_location(discrete_variable var) const {
    return var_location(var.id());
  }

  size_t var_location(size_t var_id) const {
    size_t location = _num_vars;
    for(size_t i = 0; i < _num_vars && !(location < _num_vars); ++i) {
      if(_vars[i].id() == var_id) location = i;
    }
    return location;
  }

  // get the index within our domain of each variable in other.
  // other must be a subset 
  std::vector<size_t> vars_location(const discrete_bounds& other) const {
    // ensure that the other domain is a subset of our domain
    DCHECK_EQ((*this + other).num_vars(), num_vars());
    
    std::vector<size_t> locations(other.num_vars());
    // NOTE this depends on the list of variables in both domains being sorted
    size_t subdomain_idx = 0;
    for(size_t i=0; i<num_vars(); ++i) {
      if(subdomain_idx < other.num_vars() && 
          var(i) == other.var(subdomain_idx)) {
        locations.at(subdomain_idx) = i;
        ++subdomain_idx;
      }
    }
    
    return locations;
  }

  // get the index within our domain of each variable in other.
  // other must be a subset 
  std::vector<size_t> vars_location(const std::vector<discrete_variable>& other) const {
    // ensure that the other domain is a subset of our domain
    DCHECK_EQ((*this + discrete_bounds(other)).num_vars(), num_vars());
    
    std::vector<size_t> locations(other.size());
    for(size_t i=0; i<other.size(); ++i) {
      locations.at(i) = var_location(other[i].id());
    }
    
    return locations;
  }

  //! determine the number of assignments
  // NOTE recomputing this every time it is needed is burdensome. the size could 
  // be cached (in a mutable var) and only recomputed when _num_vars changes. 
  size_t size() const { 
    size_t sum = 0;
    if(num_vars() > 0) {
      sum = 1;
      for(size_t i = 0; i < num_vars(); ++i) {
        // ensure variables to be sorted order
        if(i > 0) { DCHECK_LT( _vars[ i-1], _vars[i] ); }
        // and have positive arity
        DCHECK_GT(_vars[i].size(), 0);
        sum *= _vars[i].size();
      }
    }
    return sum;
  }

  void load(graphlab::iarchive& arc) {
    arc >> _num_vars;
    DCHECK_LE(_num_vars, MAX_DIM);
    for(size_t i = 0; i < _num_vars; ++i) arc >> _vars[i];
  }
    
  void save(graphlab::oarchive& arc) const {
    arc << _num_vars;
    for(size_t i = 0; i < _num_vars; ++i) arc << _vars[i];
  }

private:
  //mutable size_t _size_cached;
  //mutable size_t _last_num_vars;
  size_t _num_vars;
  // REVIEW C style array is space inefficient. a vector might be 
  // better. it would have to be resized if the domain was modified; 
  // however, it would remove the need for the template parameter.
  //   i tried this (37fef6f16b6d). creating the vector on the heap is 
  //   far too slow 
  discrete_variable _vars[MAX_DIM];
};


template<size_t MAX_DIM>
std::ostream& operator<<(std::ostream& out,
                   const discrete_bounds<MAX_DIM>& dom) {
  out << "{";
  for(size_t i = 0; i < dom.num_vars(); ++i) {
    out << dom.var(i) << "[0:" << dom.var(i).size()-1 << "]";
    if( i < dom.num_vars()-1 ) out << ", ";
  }
  return out << "} ";  
}


}; // end of namespace graphlab



#include <graphlab/macros_undef.hpp>

#endif // DISCRETE_BOUNDS_HPP
