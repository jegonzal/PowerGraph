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


#ifndef SPARSE_TABLE_HPP
#define SPARSE_TABLE_HPP

#include <stdint.h>
#include <assert.h>

#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <map>

#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/logger/assertions.hpp>

#include "table_base.hpp"
#include "dense_table.hpp"
#include "sparse_index.hpp"


namespace graphlab {


/**
 * An n-D sparse table up to max_dim dimensions. 
 * SEE dense_table.hpp for more detail
 *
 * \author Scott Richardson     10/2012
 */
template<size_t MAX_DIM>
class sparse_table : public table_base<MAX_DIM> {
private:
  typedef sparse_table const *const    const_ptr;
  typedef table_base<MAX_DIM>          table_base_t;

  typedef discrete_variable            variable_t;
  typedef discrete_domain<MAX_DIM>     domain_t;
  typedef discrete_assignment<MAX_DIM> assignment_t;
  typedef dense_table<MAX_DIM>         dense_table_t;

public:
  typedef sparse_index<MAX_DIM> sparse_index_t;
  typedef std::vector<std::pair<sparse_index_t, double> > compact_data_t;
private:
  typedef std::map<sparse_index_t, double> sparse_data_t;
  typedef std::vector<std::pair<const sparse_index_t*, double*> >       compact_view_t;
  typedef std::vector<std::pair<const sparse_index_t*, const double*> > compact_const_view_t;

public:
  /** Construct an empty table */
  sparse_table() { } 

  /** Construct a table over the given domain */
  sparse_table(const domain_t& dom) {
    set_domain(dom);
  }

  /** Construct a table over the given domain 
   * dom  : the domain over which the table is defined 
   * data : a vector of assignment-value pairs. the assignment must
   *        be sorted according to dom; that is, such that the 
   *        variable with the smallest id iterates fastest 
   */
  sparse_table(const domain_t& dom, const sparse_data_t& data) {
    set_domain(dom);
    _dataAtAsg = data;
  }
  
  /** Construct a table over the given domain 
   * vars : a vector of variables that compose the domain
   * data : a vector of values serialized such that the first 
   *        variable in vars iterates the fastest
   * NOTE this is a convenience constructor. the entries in the 
   * vector are re-sorted such that the variable with the smallest
   * id iterates fastest
   */
  sparse_table(const std::vector<variable_t>& vars, 
      const std::vector<std::pair<size_t, double> >& data) 
  {
    // Construct the arguments (which will remap the domain)
    set_domain(domain_t(vars));

    // create a faux domain with the size of the dimensions ordered correctly. this
    // is essentially a permute operation.
    domain_t dom;
    for(size_t i=0; i<vars.size(); ++i) {
      domain_t d1(variable_t(i, vars[i].size()));
      dom += d1;
    }

    for(size_t i=0; i < data.size(); ++i) { 
      size_t sparse_idx = data[i].first;
      assignment_t asg(dom, sparse_idx);

      // permute the assignment
      std::vector<size_t> asgs(asg.begin(), asg.end());
      assignment_t fast_asg(vars, asgs);
      set_logP(fast_asg, data[i].second);
    }
  }

  /** Construct an empty table over the given variable */
  sparse_table(const variable_t& args) { 
    // Construct the arguments (which will remap the domain)
    set_domain(domain_t(args));
  }
  
  /** Construct an empty table over the given domain */
  sparse_table(const std::vector<variable_t>& args) { 
    // Construct the arguments (which will remap the domain)
    set_domain(domain_t(args));
  }

// NOTE currently, implementing the (big) three isnt strictly necessary
  /** Construct a copy */
  sparse_table(const sparse_table& other) : 
      _args(other._args), _dataAtAsg(other._dataAtAsg) { }

  /** Destructor */
  virtual ~sparse_table() { }

  // REVIEW currently, this isnt necessary
  /** Standard assignment operator */
  sparse_table& operator=(const sparse_table& other) {
    if(this == &other) 
      return *this;
  
    _args = other._args;
    //_dataAtAsg.insert(other._dataAtAsg.begin(), other._dataAtAsg.end());
    _dataAtAsg = other._dataAtAsg;
    return *this;
  }

public: 
  using table_base_t::APPROX_LOG_ZERO;

  // if the data structures between the two tables is equivilent, this is faster
  sparse_table& copy_onto(const sparse_table& other) {
    if(this == &other) 
      return *this;
  
    // ensure the domains are the same 
    DCHECK_EQ(args(), other.args());
    // ensure the number of non-zero entries are the same (sanity check)
    DCHECK_EQ(_dataAtAsg.size(), other._dataAtAsg.size());

    typename sparse_data_t::iterator entry     = _dataAtAsg.begin();
    typename sparse_data_t::const_iterator end = _dataAtAsg.end();
    typename sparse_data_t::const_iterator oentry = other._dataAtAsg.begin();
    for( ; entry != end; ++entry, ++oentry) {
      // ensure the two assignments are equivilent (std::map should sort them similarly)
      DCHECK_EQ(entry->first, oentry->first);
      __set_logP(entry->second, oentry->second);
    }
    // slower
    //_dataAtAsg.insert(other._dataAtAsg.begin(), other._dataAtAsg.end());

    return *this;
  }

  /** 
   * Reset the domain for the table. A domain is defined by a vector of
   * variables, and an assignment is defined over that domain. 
   */
  void set_domain(const domain_t& args) {
    _args = args;
    _dataAtAsg.clear();
  }
  const domain_t& domain() const {
    return args();
  }

  bool operator==(const sparse_table& other) {
    // are the two domains equal
    if(args() != other.args()) return false;
    // are there the same number of non-zero elements in the two tables
    if(_dataAtAsg.size() != other._dataAtAsg.size()) return false;

    typename sparse_data_t::iterator entry     = _dataAtAsg.begin();
    typename sparse_data_t::const_iterator end = _dataAtAsg.end();
    typename sparse_data_t::const_iterator oentry = other._dataAtAsg.begin();
    for( ; entry != end; ++entry, ++oentry) {
      // is the assignment the same (std::map should sort them similarly)
      if(entry->first != oentry->first) return false;
      // is the value the same
      if(entry->second != oentry->second) return false;
    }
    return true;
  }
  bool operator!=(const sparse_table& other) {
    return !this->operator==(other);
  }

  /** 
   * Return the variable at the given index within the domain. (var(i) 
   * specifies the dimension associated with the i'th element of an 
   * assignment, i.e., sparse_index::_asg[i])
   */
  virtual const variable_t& var(const size_t index) const {
    return args().var(index);
  }
  /** 
   * Return the index for a given variable within the domain (as well as
   * into sparse_index::_asg[]).
   */
  size_t var_location(const variable_t& var) {
    return args().var_location(var);
  }

  /** Return the number of dimensions in the domain */
  virtual size_t ndims() const { return args().num_vars(); }
  /** Return the number of elements in the domain: prod(size(table)) */
  virtual size_t numel() const { return args().size(); } 
  /** Return the number of non-zero elements in the table */
  size_t nnz() const { return _dataAtAsg.size(); }

  /** Zero existing entries in the table */
  virtual void zero() {
    typename sparse_data_t::iterator entry     = _dataAtAsg.begin();
    typename sparse_data_t::const_iterator end = _dataAtAsg.end();
    for( ; entry != end; ++entry) {
      entry->second = 0.0;
    }
  }
  
private:

  inline void remove_logP(const sparse_index_t& asg) {
    // the assignment must be within the domain
    DASSERT_TRUE(validate_asg(asg)); 

    if(_dataAtAsg.count(asg) == 0) return;

    _dataAtAsg.erase(asg);
  }

  // REVIEW i dont love this method, but it does afford me some bounds checking. 
  // set_logP(sparse_index_t&, double), seems cleaner, but it has to 
  // re-lookup the pointer. very slow. 
  inline void __set_logP(double& tbl_ref, const double& val) {
    tbl_ref = std::max(val, APPROX_LOG_ZERO());
  }

  inline void set_logP(const size_t linear_index, const double& val) {
    set_logP(compute_asg(linear_index), val);
  }
  inline void set_logP(const std::vector<size_t>& asg, const double& val) {
    set_logP(sparse_index_t(asg), val); 
  }
  inline void set_logP(const sparse_index_t& asg, const double& val) {
    // the assignment must be within the domain
    DASSERT_TRUE(validate_asg(asg)); 

    _dataAtAsg[asg] = std::max(val, APPROX_LOG_ZERO());
  }

  inline double logP(const sparse_index_t& asg) const {
    DASSERT_TRUE(validate_asg(asg)); 

    // O(log(n))
    typename sparse_data_t::const_iterator val = _dataAtAsg.find(asg);
    return val == _dataAtAsg.end() ? APPROX_LOG_ZERO() : val->second;
  }

public:
  /** 
   *  Add an entry to the sparse table (indexed by its sparse_index).
   *  Clip values to be greater than or equal to APPROX_LOG_ZERO.
   */
  // NOTE the assignment is not removed from the domain if val is APPROX_LOG_ZERO. 
  // in the future, if these values are removed, it could invalidate any iterator
  // over the list of sparse assignments.
  inline void set_logP(const assignment_t& asg, const double& val) {
    DCHECK_EQ(asg.args(), args());
    set_logP(as_sparse_index(asg), val);
  }
  /** Remove an entry from the sparse table and its corresponding sparse_index) */ 
  inline void remove_logP(const assignment_t& asg) {
    DCHECK_EQ(asg.args(), args());
    remove_logP(as_sparse_index(asg));
  }
  // NOTE index is serialized according to the linear indexing of the domain
  // TODO can i make this private? 
  inline double logP(const size_t linear_index) const {
    return logP(compute_asg(linear_index));
  }
  /** Return an entry from the sparse table (indexed by its sparse_index) */ 
  inline double logP(const assignment_t& asg) const {
    DCHECK_EQ(asg.args(), args());
    return logP(as_sparse_index(asg));
  }

  //! this(x) /= other(x);
  // supports broadcasting of a sub-domain across the full domain 
  sparse_table& operator/=(const dense_table_t& other) {
    return for_each_assignment(other, divides());
    return *this;
  }

  //! this(x) *= other(x);
  // supports broadcasting of a sub-domain across the full domain 
  sparse_table& operator*=(const dense_table_t& other) {
    return for_each_assignment(other, multiplies());
  }

  //! this(x) /= other(x);
  // supports broadcasting of a sub-domain across the full domain 
  sparse_table& operator/=(const sparse_table& other) {
    return for_each_assignment(other, divides());
  }
    //! this(x) *= other(x);
  // supports broadcasting of a sub-domain across the full domain 
  sparse_table& operator*=(const sparse_table& other) {
    return for_each_assignment(other, multiplies());
  }

private:
  struct divides {
    inline double operator()(const double& a, const double& b) const {
      return a - b; 
    }
  };
  struct multiplies {
    inline double operator()(const double& a, const double& b) const {
      return a + b; 
    }
  };

  template<class Func>
  inline sparse_table& for_each_assignment(const dense_table_t& other, const Func& f) {
    // other domain must be a subset of this domain
    DCHECK_EQ((args() + other.args()).num_vars(), args().num_vars());

    assignment_t dense_asg(args()); 
    // only need to operate on the the assignments in the sparse table
    // (equivalently, the intersection of the sparse and dense assignments)
    typename sparse_data_t::iterator it        = _dataAtAsg.begin();
    typename sparse_data_t::const_iterator end = _dataAtAsg.end();
    for( ; it != end; ++it) {
      dense_asg.set_index(linear_index(it->first));
      //double val = it->second + other.logP(dense_asg));
      //it->second = val;
      double val = f(it->second, other.logP(dense_asg));
      __set_logP(it->second, val);
    }
    return *this;
  }

  template<class Func>
  sparse_table& for_each_assignment(const sparse_table& other, const Func& f) {

    // if the tables span the same domain
    if(args() == other.args()) {
      DCHECK_EQ(numel(), other.numel());
      
      // NOTE the assignments NOT in the intersection of the two sparse tables 
      // will be APPROX_LOG_ZERO() and are removed in *this
      intersect(other);
      typename sparse_data_t::iterator it        = _dataAtAsg.begin();
      typename sparse_data_t::const_iterator end = _dataAtAsg.end();
      typename sparse_data_t::const_iterator other_it = other._dataAtAsg.begin();
      for( ; it != end; ++it, ++other_it) {
        //double val = it->second + other_it->second);
        //it->second = val;
        double val = f(it->second, other_it->second);
        __set_logP(it->second, val);
      }
    }
    // else, broadcast the sub-domain across the full domain
    else { 
      // other domain must be a subset of this domain
      DCHECK_EQ((args() + other.args()).num_vars(), args().num_vars());
      
      compact_view_t compact_view = as_vector_view();
      compact_const_view_t other_compact_view = other.as_vector_view();
      // define the one-to-one mapping from other's domain to our's
      std::vector<size_t> sorting_inds = args().vars_location(other.args());
      // reorder the assignments so they can be quickly iterated over
      permute(sorting_inds, compact_view);
      other.permute(other_compact_view);
      
      // Loop over x
      // NOTE the assignments are sorted the same. ie. our assignments share the same 
      // ordering over the sub-domain spaned by msg as the assignments in msg.
      typename compact_const_view_t::const_iterator x_fastasg = other_compact_view.begin();
      typename compact_const_view_t::const_iterator x_end     = other_compact_view.end();
      typename compact_view_t::iterator       y_fastasg = compact_view.begin();
      typename compact_view_t::const_iterator y_end     = compact_view.end();
      sparse_index_t yasg;
      for( ; x_fastasg < x_end; ++x_fastasg) {
        while(y_fastasg != y_end) {
          yasg = restrict(*(y_fastasg->first), other.args());
          if(*(x_fastasg->first) > yasg) { 
            ++y_fastasg;
            continue;
          } 
          else if(*(x_fastasg->first) < yasg) {
            ++x_fastasg;
            break;
          }
          // else the sub-assignments are equal
          else {
            //double val = *(y_fastasg->second) + *(x_fastasg->second);
            //*(y_fastasg->second) = val;
            double val = f(*(y_fastasg->second), *(x_fastasg->second));
            __set_logP(*(y_fastasg->second), val);
            ++y_fastasg;
          }
        }
      }
    }
    return *this;
  }

public:
  using table_base_t::MAP;

  // since the message is always a unary distribution, this is basically 
  // >>> max(reshape(permute(
  //     cavity, circshift(1:ndims(cavity), [1 -msg.dim])), [], msg.numel),
  //     [], 1)
  // or more generally, 
  // >>> max(reshape(permute(
  //     cavity, [setdiff(1:ndims(cavity), msg.dims), msg.dims]), [], msg.numel), 
  //     [], 1)
  void MAP(dense_table_t& msg) const {
    // No need to marginalize
    if(args() == msg.args()) {
      // Just copy and return
      as_dense_table(msg);
      return;
    }
    
    // the domains cannot be disjoint
    DCHECK_GT((args() - msg.args()).num_vars(), 0);
  
    compact_const_view_t fast_view = as_vector_view();
    // define the one-to-one mapping from the msg's domain to our's
    std::vector<size_t> sorting_inds = args().vars_location(msg.args());
    // reorder the assignments so they can be quickly iterated over
    permute(sorting_inds, fast_view); 

    assignment_t yasg(args());
    // Loop over x
    // NOTE our assignments have been reordered so we can index assignments in 
    // the two domains consecutively. e.g., if the domain of msg, {v1,v2}, is 
    // sorted in ascending order, then our assignments must also be sorted in 
    // assending order over {v1,v2} (although these sub-domains need not be 
    // be sorted the same.)
    typename compact_const_view_t::const_iterator fastyasg = fast_view.begin();
    typename compact_const_view_t::const_iterator yend     = fast_view.end();
    typename domain_t::const_iterator xasg = msg.args().begin();
    typename domain_t::const_iterator xend = msg.args().end();
    for( ; xasg != xend; ++xasg) {
      double maxval = APPROX_LOG_ZERO();

      // loop over y
      while(fastyasg != yend) {
        yasg.set_index(linear_index(*(fastyasg->first))); 
        if(*xasg != yasg.restrict(xasg->args())) break;
        
        //maxval = std::max(maxval, _dataAtAsg[*fastyasg]);
        maxval = std::max(maxval, *(fastyasg->second));
        ++fastyasg;
      }
      msg.set_logP( *xasg, maxval );
    }
  }

  void intersect(const sparse_table& other) {
    map_left_intersection(_dataAtAsg, other._dataAtAsg);
  }

private:
  //! Compute the index from the sparse_index
  // NOTE index is serialized according to the linear indexing of the domain
  size_t linear_index(const sparse_index_t& asg) const {
    size_t multiple = 1;
    // Clear the index
    size_t index = 0;
    for(size_t i = 0; i < args().num_vars(); ++i) {
      index += multiple * asg.asg_at(i);
      // assert(args().var(i).nasgs > 0);
      multiple *= args().var(i).size();
    }
    return index;
  }

  /** Ensure that an asg falls within the domain */
  bool validate_asg(const sparse_index_t& asg) const {
    // no index can be larger than the number of labels in that dimension
    //return asg <= end_asg();
    for(size_t i=0; i<args().num_vars(); ++i)
      if(asg.asg_at(i) >= args().var(i).size()) return false;
    return true;
  }

  //! Compute the sparse_index from the index
  // NOTE index is serialized according to the linear indexing of the domain
  sparse_index_t compute_asg(const size_t index) const {
    DCHECK_LT(index, args().size());
  
    sparse_index_t asg(args().num_vars());
    size_t quotient = index;
    for(size_t i = 0; i < args().num_vars(); ++i) {
      asg.set_asg_at(i, quotient % args().var(i).size());
      quotient /= args().var(i).size();
      // assert(asg.asg_at(i) < args().var(i).size());
    }
    return asg;
  }

public:
  //! Compute the largest assignment possible
  assignment_t end_asg() {
    sparse_index_t asg;
    for(size_t i=0; i<args().num_vars(); ++i)
      asg.set_asg_at(i, args().var(i).size() - 1);
    return as_assignment(asg);
  }

  //! WARNING this could lead to a very large table
  void as_dense_table(dense_table_t& other) const {
    other.set_domain(args());

    typename domain_t::const_iterator asg = other.args().begin();
    typename domain_t::const_iterator end = other.args().end();
    for( ; asg != end; ++asg) {
      other.set_logP( *asg, logP(*asg) );
    }
  }

  compact_data_t as_vector() const {
    compact_data_t compact_data;
    compact_data.resize(_dataAtAsg.size());

    typename sparse_data_t::const_iterator it = _dataAtAsg.begin();
    for(size_t i = 0; i < _dataAtAsg.size(); ++i, ++it) {
      compact_data.at(i) = *it;
    }
    //std::copy(_dataAtAsg.begin(), _dataAtAsg.end(), compact_data.begin());
    return compact_data;
  }

  std::vector<sparse_index_t> keyset() const {
    std::vector<sparse_index_t> keys;
  
    typename sparse_data_t::const_iterator it  = _dataAtAsg.begin();
    typename sparse_data_t::const_iterator end = _dataAtAsg.end();
    for( ; it != end; ++it) {
      keys.push_back(it->first);
    }
    return keys;
  }

//: virtual methods
public:
  virtual sparse_table& deep_copy(const table_base_t& base) {
    if(this == &base) return *this;
  
    // ensure we are dealing with a sparse_table
    const_ptr other = dynamic_cast<const_ptr>(&base);
    if(other == NULL) {
      std::cout << "ERROR: std::bad_cast" << std::endl;
      // REVIEW should probably raise an exception
      ASSERT_TRUE(false);
    }
  
    *this = *other;
    return *this;
  }
  virtual sparse_table& copy_onto(const table_base_t& base) {
    if(this == &base) return *this;
  
    // ensure we are dealing with a sparse_table
    const_ptr other = dynamic_cast<const_ptr>(&base);
    if(other == NULL) {
      std::cout << "ERROR: std::bad_cast" << std::endl;
      // REVIEW should probably raise an exception
      ASSERT_TRUE(false);
    }
  
    this->copy_onto(*other);
    return *this;
  }
/*  
  // NOTE this operation would turn a sparse table into a dense table
  //! this(x) += other(x);
  virtual sparse_table& plus_equals(const table_base_t& base) {
    // ensure we are dealing with a sparse_table
    const_ptr other = dynamic_cast<const_ptr>(&base);
    if(other == NULL) {
      std::cout << "ERROR: std::bad_cast" << std::endl;
      // REVIEW should probably raise an exception
      ASSERT_TRUE(false);
    }
  
    // TODO implement operator
    *this += *other;
  
    return *this;
  }
*/
  //! this(x) *= other(x);
  virtual sparse_table& times_equals(const table_base_t& base) {
    // ensure we are dealing with a sparse_table
    {
      sparse_table const* const other = 
        dynamic_cast<sparse_table const* const>(&base);
      if( NULL !=  other) {
        *this *= *other;
        return *this;
      }
    } {
      dense_table_t const* const other = 
        dynamic_cast<dense_table_t const* const>(&base);
      if( NULL != other ) { 
        *this *= *other;
        return *this;
      }
    }
    std::cout << "ERROR: std::bad_cast" << std::endl;
    // REVIEW should probably raise an exception
    ASSERT_TRUE(false);
  }
  
  //! this(x) /= other(x);
  virtual sparse_table& divide_equals(const table_base_t& base) {
    // ensure we are dealing with a sparse_table
    {
      sparse_table const* const other = 
          dynamic_cast<sparse_table const* const>(&base);
      if( NULL !=  other) {
        *this /= *other;
        return *this;
      }
    } {
      dense_table_t const* const other = 
          dynamic_cast<dense_table_t const* const>(&base);
      if( NULL != other ) { 
        *this /= *other;
        return *this;
      }
   }
    std::cout << "ERROR: std::bad_cast" << std::endl;
    // REVIEW should probably raise an exception
    ASSERT_TRUE(false);
  }
  
  //! (out(x) = this(x)) * other(x);
  virtual void times(const table_base_t& base, 
      table_base_t& out_base) const {
  
    // ensure we are dealing with a sparse_table
    sparse_table *const out = 
        dynamic_cast<sparse_table *const>(&out_base);
    if(out == NULL) {
      std::cout << "ERROR: std::bad_cast" << std::endl;
      // REVIEW should probably raise an exception
      ASSERT_TRUE(false);
    }
  
    *out = *this; // deep copy
    out->times_equals(base);
  }
  
  //! (out(x) = this(x)) / other(x);
  virtual void divide(const table_base_t& base, 
      table_base_t& out_base) const {
  
    // ensure we are dealing with a sparse_table
    sparse_table *const out = 
        dynamic_cast<sparse_table *const>(&out_base);
    if(out == NULL) {
      std::cout << "ERROR: std::bad_cast" << std::endl;
      // REVIEW should probably raise an exception
      ASSERT_TRUE(false);
    }
  
    *out = *this; // deep copy
    out->divide_equals(base);
  }
  
  virtual void MAP(table_base_t& base) const {
    // ensure we are dealing with a dense_table
    dense_table_t* msg = dynamic_cast<dense_table_t*>(&base);
    if(msg == NULL) {
      std::cout << "ERROR: std::bad_cast" << std::endl;
      // REVIEW should probably raise an exception
      ASSERT_TRUE(false);
    }
  
    MAP(*msg);
  }

  virtual std::ostream& print(std::ostream& out = std::cout) const {
    // ensure we are dealing with a sparse_table
    const_ptr tbl = dynamic_cast<const_ptr>(this);
    if(tbl == NULL) {
      std::cout << "ERROR: std::bad_cast" << std::endl;
      // REVIEW should probably raise an exception
      ASSERT_TRUE(false);
    }
  
    out << *tbl;
    return out;
  }

  friend std::ostream& operator<<(std::ostream& out,
      const sparse_table<MAX_DIM>& tbl) 
  {
    out << "Sparse Table: " << tbl.args() << "{" << std::endl;
    typename sparse_data_t::const_iterator val = tbl._dataAtAsg.begin();
    typename sparse_data_t::const_iterator end = tbl._dataAtAsg.end();
    for( ; val != end; ++val) {
      out << "\tLogP({" << val->first << "}=" << tbl.linear_index(val->first) << ")=" << val->second << std::endl;
    }
    out << "}";
    return out;
  }
  
  virtual void load(graphlab::iarchive& arc) {
    arc >> _args;
    arc >> _dataAtAsg; // uses graphlab serialization operator
  }
  virtual void save(graphlab::oarchive& arc) const {
    arc << _args;
    arc << _dataAtAsg; // uses graphlab serialization operator
  }

private:
  //! Tests whether one sparse_index's i'th index is < another's 
  struct less_than_by_index {
    less_than_by_index(size_t ind) : _sorting_ind(ind) { } 
  
    inline bool operator()(
        const std::pair<sparse_index_t, double>& a, 
        const std::pair<sparse_index_t, double>& b) const
    {
      return less_than(a.first, b.first);
    }
    inline bool operator()(
        const std::pair<const sparse_index_t*, double*>& a, 
        const std::pair<const sparse_index_t*, double*>& b) const
    {
      return less_than(*(a.first), *(b.first));
    }
    inline bool operator()(
        const std::pair<const sparse_index_t*, const double*>& a, 
        const std::pair<const sparse_index_t*, const double*>& b) const
    {
      return less_than(*(a.first), *(b.first));
    }
  private: 
    inline bool less_than(const sparse_index_t& a, 
        const sparse_index_t& b) const 
    {
      return a.asg_at(_sorting_ind) < b.asg_at(_sorting_ind);
    }
    size_t _sorting_ind;
  };
  template<class T>
  void permute(T& data) const {
    std::vector<size_t> sorting_inds;
    for(size_t i=0; i<args().num_vars(); ++i)
      sorting_inds.push_back(i);
    permute(sorting_inds, data);
  }
  template<class T>
  void permute( const size_t sorting_ind, T& data ) const {
    std::vector<size_t> sorting_inds;
    sorting_inds.push_back(sorting_ind);
    permute(sorting_inds, data);
  }
  // REVIEW i might be able to copy _dataAtAsg and sort it with a new predicate
  template<class T>
  void permute( const std::vector<size_t>& sorting_inds, T& data ) const {
    std::vector<size_t>::const_reverse_iterator s    = sorting_inds.rbegin();
    std::vector<size_t>::const_reverse_iterator rend = sorting_inds.rend();
    for( ; s != rend; ++s) {
      DCHECK_LT(*s, args().num_vars());
  
      std::stable_sort(data.begin(), data.end(), less_than_by_index(*s));
    }
  }

  //! Restrict the sparse_index to a sparse_index over the subdomain
  sparse_index_t restrict(const sparse_index_t& asg, 
      const domain_t& sub_domain) const 
  {
    sparse_index_t other_asg(sub_domain.num_vars());
    size_t index = 0;
    // Map the variables 
    for(size_t i = 0; i < args().num_vars() && 
        index < sub_domain.num_vars(); ++i) {
      if(sub_domain.var(index) == args().var(i)) {
        other_asg.set_asg_at(index, asg.asg_at(i));
        index++;
      }
    }
    DCHECK_EQ(index, sub_domain.num_vars());
    return other_asg;
  } // end of restrict
  
  // O(n)
  // from http://stackoverflow.com/questions/3772664/intersection-of-two-stl-maps 
  // and http://stackoverflow.com/questions/1773526/in-place-c-set-intersection
  // REVIEW is there any problem with invalidated iterators? 
  template<typename KeyType, typename LeftValue, typename RightValue>
  void map_left_intersection(
      std::map<KeyType, LeftValue>& left, 
      const std::map<KeyType, RightValue>& right) const
  {
    typename std::map<KeyType, LeftValue>::iterator    il = left.begin();
    typename std::map<KeyType, LeftValue>::iterator l_end = left.end();
    typename std::map<KeyType, RightValue>::const_iterator ir    = right.begin();
    typename std::map<KeyType, RightValue>::const_iterator r_end = right.end();
    while (il != l_end && ir != r_end) {
      if (il->first < ir->first) {
        left.erase(il);
        ++il;
      }
      else if (ir->first < il->first) {
        ++ir;
      }
      else {
        ++il;
        ++ir;
      }
    }
    left.erase(il, l_end);
  }

  compact_view_t as_vector_view() {
    compact_view_t compact_view;
    compact_view.resize(_dataAtAsg.size());

    typename sparse_data_t::iterator it = _dataAtAsg.begin();
    for(size_t i = 0; i < _dataAtAsg.size(); ++i, ++it) {
      compact_view.at(i) = std::make_pair(&(it->first), &(it->second));
    }
    //std::copy(_dataAtAsg.begin(), _dataAtAsg.end(), compact_view.begin());
    return compact_view;
  }
  compact_const_view_t as_vector_view() const {
    compact_const_view_t compact_view;
    compact_view.resize(_dataAtAsg.size());

    typename sparse_data_t::const_iterator it = _dataAtAsg.begin();
    for(size_t i = 0; i < _dataAtAsg.size(); ++i, ++it) {
      compact_view.at(i) = std::make_pair(&(it->first), &(it->second));
    }
    //std::copy(_dataAtAsg.begin(), _dataAtAsg.end(), compact_view.begin());
    return compact_view;
  }

private:
  inline const domain_t& args() const { return _args; }

  inline assignment_t as_assignment(const sparse_index_t &asg) const {
    std::vector<size_t> asgs(asg.begin(), asg.end());
    return assignment_t(args(), asgs);
  }

  inline sparse_index_t as_sparse_index(const assignment_t &asg) const {
    std::vector<size_t> asgs(asg.begin(), asg.end());
    return sparse_index_t(asgs);
  }

private:
  //! The indicies in an assignment are mapped (one-to-one and in-order) to the 
  //  variables in a domain.
  domain_t _args;

  //! Map between the sparse assignment in the domain and its value. Sorted by
  //  assignment.
  sparse_data_t  _dataAtAsg;
};


} // end of namespace graphlab

#endif // SPARSE_TABLE_HPP
