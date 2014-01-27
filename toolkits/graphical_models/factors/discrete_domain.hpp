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


#ifndef DISCRETE_DOMAIN_HPP
#define DISCRETE_DOMAIN_HPP

#include <graphlab/logger/assertions.hpp>

#include "discrete_variable.hpp"
#include "discrete_bounds.hpp"
#include "discrete_assignment.hpp"

#include <graphlab/macros_def.hpp>
namespace graphlab {


/**
 * This class respresents a discrete domain over a set of variables.
 *
 * \author Scott Richardson     4/2013
 */
template<size_t MAX_DIM>
class discrete_domain : public discrete_bounds<MAX_DIM> {
  typedef discrete_bounds<MAX_DIM>     bounds_type;
  typedef discrete_assignment<MAX_DIM> assignment_type;

public:
  //! Make an empty domain
  discrete_domain() : 
      bounds_type() { }

  //! Make a single variable discrete_domain
  discrete_domain(const discrete_variable& v1) :
      bounds_type(v1) { }

  //! Make a two variable discrete_domain
  discrete_domain(const discrete_variable& v1, const discrete_variable& v2) :
      bounds_type(v1, v2) { }

  //! Make a three variable discrete_domain
  discrete_domain(const discrete_variable& v1,
                  const discrete_variable& v2,
                  const discrete_variable& v3) :
      bounds_type(v1, v2, v3) { }

  //! Make a discrete_domain from a vector of variables
  explicit discrete_domain(const std::vector<discrete_variable>& variables) :
      bounds_type(variables) { }

  //! Make a discrete_domain from a set of variables
  explicit discrete_domain(const std::set<discrete_variable>& variables) :
      bounds_type(variables) { }
  
  discrete_domain(const discrete_domain& other) : 
      bounds_type(other) { }

  discrete_domain(const bounds_type& other) : 
      bounds_type(other) { }

  virtual ~discrete_domain() { } 

  /** Standard assignment operator */
  discrete_domain& operator=(const discrete_domain& other) {
    if(this == &other) 
      return *this;
    
    bounds_type::operator=(other);
    return *this;
  }

  class ConstIterator;

  // Iterators 
  // from http://www.oreillynet.com/pub/a/network/2005/11/21/what-is-iterator-in-c-plus-plus-part2.html
  // and http://www.cs.helsinki.fi/u/tpkarkka/alglib/k06/lectures/Iterators.html
  // although i think i should use this reference: 
  // http://www.drdobbs.com/the-standard-librarian-defining-iterato/184401331?pgno=3
  class Iterator : 
      public std::iterator<std::forward_iterator_tag, assignment_type >
  {
  private:
    friend class ConstIterator;
    typedef discrete_domain<MAX_DIM> domain_type;

    typedef std::iterator<std::forward_iterator_tag, assignment_type > iterator_t;
    // "using typename" doesn't actually work in GCC < 4.7, which we don't have installed
    // everywhere. The "typedef typename" construct seems to, so stick with that for now.
    //using typename iterator_t::value_type;
    //using typename iterator_t::reference;
    //using typename iterator_t::pointer;
    typedef typename iterator_t::value_type   value_type_custom;
    typedef typename iterator_t::reference    reference_custom;
    typedef typename iterator_t::pointer      pointer_custom;

    assignment_type _asg;

  public:
    explicit Iterator(const domain_type& dom) {
      // initilize a new assignment
      _asg = assignment_type(dom);
    }
    explicit Iterator(const assignment_type& asg) : _asg(asg) { }
    Iterator& operator=(const Iterator& other) {
      if(this == &other) return *this;

      _asg = other._asg;
      return *this;
    }
    // implicit copy constructor, copy assignment and destructor
    bool operator==(const Iterator& other) const {
      return _asg == other._asg;
    }
    bool operator!=(const Iterator& other) const {
      return !(*this == other);
    }
    reference_custom operator*() {
      return _asg;
    }
    pointer_custom operator->() {
      // this may be more correct, but less clear
      //return &*(domain_type::Iterator)*this;
      return &_asg;
    }
    Iterator& operator++() {
      ++_asg;
      return *this;
    }
    Iterator operator++(int) {
      Iterator orig = *this; 
      ++(*this); 
      return orig;
    }
  };

  class ConstIterator : 
      public std::iterator<std::forward_iterator_tag, const assignment_type > 
  {
    typedef discrete_domain<MAX_DIM> domain_type;

    typedef std::iterator<std::forward_iterator_tag, const assignment_type > const_iterator_t;
    //using typename const_iterator_t::value_type;
    //using typename const_iterator_t::reference;
    //using typename const_iterator_t::pointer;
    typedef typename const_iterator_t::value_type   value_type_custom;
    typedef typename const_iterator_t::reference    reference_custom;
    typedef typename const_iterator_t::pointer      pointer_custom;

    assignment_type _asg;

  public:
    explicit ConstIterator(const domain_type& dom) {
      // initilize a new assignment
      _asg = assignment_type(dom);
    }
    explicit ConstIterator(const assignment_type& asg) : _asg(asg) { }
    ConstIterator(const Iterator& other) : _asg(other._asg) { }
    // implicit copy constructor and destructor
    ConstIterator& operator=(const ConstIterator& other) {
      if(this == &other) return *this;

      _asg = other._asg;
      return *this;
    }
    bool operator==(const ConstIterator& other) const {
      return _asg == other._asg;
    }
    bool operator!=(const ConstIterator& other) const {
      return !(*this == other);
    }
    reference_custom operator*() const {
      return _asg;
    }
    pointer_custom operator->() const {
      return &_asg;
    }
    ConstIterator& operator++() {
      ++_asg;
      return *this;
    }
    ConstIterator operator++(int) {
      ConstIterator orig = *this; 
      ++(*this); 
      return orig;
    }
  };

  Iterator begin() const { 
    return Iterator(*this); 
  }

  Iterator end() const { 
    Iterator ret(*this);
    ret->make_end();
    return ret;
  }

public:
  typedef Iterator       iterator;
  typedef ConstIterator  const_iterator;
};

}; // end of namespace graphlab



#include <graphlab/macros_undef.hpp>

#endif // DISCRETE_DOMAIN_HPP
