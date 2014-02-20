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


#ifndef DENSE_TABLE_HPP
#define DENSE_TABLE_HPP


/**
 * This file contains the definitions of some of the basic factor
 * types needed for loopy belief propagation. 
 *
 *  \author Joseph Gonzalez
 *  \author Scott Richardson     09/2012
 *          
 */



// INCLUDES ===================================================================>

// Including Standard Libraries
#include <cassert>
#include <cmath>

#include <iostream>
#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>
#include <set>


// Random number generation
#include <graphlab/util/random.hpp>

#include <graphlab/serialization/serialization_includes.hpp>

#include "discrete_variable.hpp"
#include "discrete_domain.hpp"
#include "discrete_assignment.hpp"
#include "fast_discrete_assignment.hpp"
#include "table_base.hpp"


// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * An n-D table up to max_dim dimensions. 
   * NOTE this table stores the data in log-space, although this 
   * implementation detail is abstracted away. E.g., operator* actually
   * adds values in log-space.
   * NOTE you can use begin()/end() to iterate over all assignments
   * in the domain. for example, a domain {var1[0,4), var2[0,5)} has two 
   * dimensions, var1 and var2, with four and five labels respectively.
   * logP()/set_logP() provide access to the underlying data via an 
   * assignment, e.g., [0,3]. the data in the table is serialized 
   * according to the linear indexing of the domain, which is ordered 
   * such that the variable with the lowest id iterates fastest.
   */
  template<size_t MAX_DIM>
  class dense_table_impl : public table_base<MAX_DIM> {
  public:
    typedef discrete_variable            variable_type;
    typedef discrete_domain<MAX_DIM>     domain_type;
    typedef discrete_assignment<MAX_DIM> assignment_type;
    typedef table_base<MAX_DIM>          table_base_t;


    /** Construct an empty table */
    dense_table_impl() { }
    
    /** Construct a table over the given domain */
    // dom : the domain over which the table is defined 
    dense_table_impl(const domain_type& dom) :
      _args(dom), _data(dom.size()) {  }

    /** Construct a dense table over the given domain */
    dense_table_impl(const std::vector<variable_type> &args) { 
      // Construct the arguments (which will remap the domain)
      set_domain(domain_type(args));
    }

    /** Construct a dense table over the given domain 
     * dom  : the domain over which the table is defined 
     * data : a vector of values serialized according to dom; that is, such
     *        that the variable with the smallest id iterates fastest 
     */
    dense_table_impl(const domain_type& dom, const std::vector<double> &data) :
      _args(dom), _data(dom.size()) 
    {
      set_data(data);
    }

    /** Construct a dense table over the given domain 
     * vars : a vector of variables that compose the domain
     * data : a vector of values serialized such that the first 
     *        variable in vars iterates the fastest
     * NOTE this is a convenience constructor. the entries in the 
     * vector are re-sorted such that the variable with the smallest
     * id iterates fastest
     * REVIEW make these static factory methods:
     * e.g., static dense_table<MAX_DIM>& table_from_serialized_data(...)
     */
    dense_table_impl(const std::vector<variable_type>& vars, 
        const std::vector<double>& data) 
    { 
      // Construct the arguments (which will remap the domain)
      set_domain(domain_type(vars));

      // create a faux domain with the size of the dimensions ordered correctly. this
      // is essentially a permute operation.
      domain_type dom;
      for(size_t i=0; i<vars.size(); ++i) {
        domain_type d1(variable_type(i, vars[i].size()));
        dom += d1;
      }

      for(size_t i=0; i < data.size(); ++i) { 
        assignment_type asg(dom, i);

        // permute the assignment
        std::vector<size_t> asgs(asg.begin(), asg.end());
        assignment_type fast_asg(vars, asgs);
        set_logP(fast_asg, data.at(i));
      }
    }

    /** Construct a unary table factor over the given var */
    dense_table_impl(const variable_type &var) { 
      // Construct the arguments (which will remap the domain)
      set_domain(domain_type(var));
    }

    /** Construct a unary dense table over the given var */
    dense_table_impl(const variable_type& var, std::vector<double>& logd) { 
      // Construct the arguments (which will remap the domain)
      set_domain(domain_type(var));

      set_data(logd);
    }

    /** Construct a unary dense table over the given var */
    dense_table_impl(const variable_type& var, double const* const begin, 
        double const* const end) 
    { 
      // Construct the arguments (which will remap the domain)
      set_domain(domain_type(var));

      set_data(begin, end - begin);
    }

  // NOTE currently, implementing the (big) three isnt strictly necessary
    /** Construct a copy */
    dense_table_impl(const dense_table_impl& other) :
      _args(other._args), _data(other._data) { }

    /** Destructor */
    virtual ~dense_table_impl() { }

    /** Standard assignment operator */
    dense_table_impl& operator=(const dense_table_impl& other) {
      if(this == &other) 
        return *this;

      _args = other._args;
      _data = other._data;
      return *this;
    }

  private:
    void set_data(const std::vector<double> &data) {
      DCHECK_EQ(_data.size(), data.size());
      DCHECK_EQ(_data.size(), _args.size());
      // i need this for copy_onto()
      //ASSERT_EQ(_args.num_vars(), 1);

      std::replace_copy_if(data.begin(), data.end(), _data.begin(), 
          isless(APPROX_LOG_ZERO()), APPROX_LOG_ZERO());
      //ASSERT_TRUE(is_finite());
    }
    void set_data(const double data[], const size_t n) {
      DCHECK_EQ(_data.size(), n);
      DCHECK_EQ(_data.size(), _args.size());
      // i need this for copy_onto()
      //ASSERT_EQ(_args.num_vars(), 1);

      std::replace_copy_if(data, data+n, _data.begin(), 
          isless(APPROX_LOG_ZERO()), APPROX_LOG_ZERO());
      //ASSERT_TRUE(is_finite());
    }

    struct isless {
      double _val;
      isless(const double& val) : _val(val) { }
      bool operator() (double number) { return number < _val; }
    };

  public:
    using table_base_t::APPROX_LOG_ZERO;

    // TODO implement operator== and operator!=

    using table_base_t::copy_onto;

    dense_table_impl& copy_onto(const dense_table_impl& other) {
      if(this == &other)
        return *this;

      DCHECK_EQ(args(), other.args());

      set_data(other._data);
      return *this;
    }

    void set_args(const domain_type& args) {
      _args = args;
      _data.resize(args.size());
    }
    void set_domain(const domain_type& args) {
      set_args(args);
    }

    inline const domain_type& args() const {
      return _args;
    }
    inline const domain_type& domain() const {
      return args();
    } 
    
    // NOTE index is serialized according to the linear indexing of the domain
    // TODO can i make this private
    inline const double& logP(const size_t index) const {
      DCHECK_LT(index, size());
      return _data[index];
    }
    const double& logP(const assignment_type& asg) const {
      if(asg.args() == args()) {
        // if the assignment index matches
        const size_t index(asg.linear_index());
        DCHECK_LT(index, size());
        return _data[index];
      } else {
        // Restrict the assignment to this domain
        const assignment_type sub_asg = asg.restrict(_args);
        DCHECK_LT(sub_asg.linear_index(), size());
        return _data[sub_asg.linear_index()];
      }
    }

  private:
    // clip values to be greater than or equal to APPROX_LOG_ZERO
    // NOTE not as efficient as exposing the value by reference, but safer.
    inline void set_logP(const size_t index, const double& val) {
      DCHECK_LT(index, size());
      _data[index] = std::max(val, APPROX_LOG_ZERO());
    }

  public:
    // NOTE be careful not to replace instances of logP(assignment_type) with 
    // logP(assignment_type.linear_index()) as these are not the same thing 
    // when the domains are different. 
    void set_logP(const assignment_type& asg, const double& val) {
      if(asg.args() == args()) {
        // if the assignment index matches
        const size_t index(asg.linear_index());
        DCHECK_LT(index, size());
        _data[index] = std::max(val, APPROX_LOG_ZERO());
      } else {
        // Restrict the assignment to this domain
        const assignment_type sub_asg = asg.restrict(_args);
        DCHECK_LT(sub_asg.linear_index(), size());
        _data[sub_asg.linear_index()] = std::max(val, APPROX_LOG_ZERO());
      }
    } // end of logP

    size_t size() const { 
      DCHECK_EQ(_args.size(), _data.size()); 
      return _args.size(); 
    }
    virtual size_t numel() const { return size(); } 

    size_t num_vars() const { return _args.num_vars(); }
    virtual size_t ndims() const { return num_vars(); } 

    virtual void zero() { std::fill(_data.begin(), _data.end(), 0); }

    void uniform() {
      std::fill(_data.begin(), _data.end(), log(1.0/size()));
    }
   
    void uniform(double value) {
      std::fill(_data.begin(), _data.end(), value);
    } 
  
    //! ensure that sum_x this(x) = 1 
    void normalize() {
      //ASSERT_TRUE(is_finite());
      // Compute the max value
      double max_value = logP(0);
      for(size_t i = 0; i < size(); ++i) {
        max_value = std::max(max_value, logP(i) );
      }
      // scale and compute normalizing constant
      double Z = 0.0;
      for(size_t i = 0; i < size(); ++i) {
        double val = logP(i) - max_value;
        set_logP(i, val);
        Z += exp(val);
      }
      // assert( !std::isinf(Z) );
      // assert( !std::isnan(Z) );
      // assert( Z > 0.0);
      const double logZ(log(Z));
      DASSERT_FALSE( std::isinf(logZ) );
      DASSERT_FALSE( std::isnan(logZ) );
      // Normalize
      for(size_t i = 0; i < size(); ++i) { set_logP( i, logP(i) - logZ ); }
      //ASSERT_TRUE(is_finite());
    } // End of normalize
    
    /** 
     * Ensure that the largest value in log form is zero.  This
     * prevents overflows on normalization. 
     */
    void shift_normalize() {
      //ASSERT_TRUE(is_finite());
      // Compute the max value
      double max_value = logP(0);
      for(size_t i = 0; i < size(); ++i) {
        max_value = std::max(max_value, logP(i));
      }
      for(size_t i = 0; i < size(); ++i) { set_logP( i, logP(i) - max_value ); }
      //ASSERT_TRUE(is_finite());
    }

    /**
     * Return false if any of the entries are not finite 
     */
    bool is_finite() const { 
      for(size_t i = 0; i < size(); ++i) {
        const bool is_inf( std::isinf( logP(i) ) );
        const bool is_nan( std::isnan( logP(i) ) );
        if( __builtin_expect( is_inf || is_nan, 0) ) return false;
      }
      return true;
    }


  public: 
    //! this(x) *= other(x);
    dense_table_impl& operator*=(const dense_table_impl& other) {
      return for_each_assignment(other, multiplies());
    }

//  //! Create a dense table on the fly
//  dense_table_impl operator*(const dense_table_impl& other) const {
//    dense_table_impl tbl = *this;
//    return tbl *= other;
//  }

    //! this(x) += other(x);
    // supports broadcasting of a sub-domain across the full domain 
    dense_table_impl& operator+=(const dense_table_impl& other) {
      return for_each_assignment(other, plus());
    }

    //! this(x) /= other(x);
    // supports broadcasting of a sub-domain across the full domain 
    dense_table_impl& operator/=(const dense_table_impl& other) {
      return for_each_assignment(other, divides());
    }

//  //! Create a dense table on the fly
//  dense_table_impl operator/(const dense_table_impl& other) const {
//    dense_table_impl tbl = *this;
//    return tbl /= other;
//  }

  private:
    // NOTE we assume we are in log space
    struct divides {
      double operator()(const double& a, const double& b) const {
        return a - b; 
      }
    };
    struct multiplies {
      double operator()(const double& a, const double& b) const {
        return a + b; 
      }
    };
    struct plus {
      double operator()(const double& a, const double& b) const {
        double out = log( exp(a) + exp(b) );
        DASSERT_FALSE(std::isinf( out ));
        DASSERT_FALSE(std::isnan( out ));
        return out;
      }
    };

    template<class Func>
    inline dense_table_impl& for_each_assignment(const dense_table_impl& other, 
        const Func& f) {
      //ASSERT_TRUE(is_finite());
      if(args() == other.args()) {
        DCHECK_EQ(size(), other.size());
        // More verctorizable version
        for(size_t i = 0; i < size(); ++i) {
          double val = f(logP(i), other.logP(i));
          set_logP( i, val );
          //logP(i) -= other.logP(i);
        }
      } else { 
        // other domain must be a subset of this domain
        DCHECK_EQ((args() + other.args()).num_vars(), num_vars());

        typename domain_type::const_iterator asg = args().begin();
        typename domain_type::const_iterator end = args().end();
        for( ; asg != end; ++asg) {
          double val = f(logP(asg->linear_index()), other.logP(*asg));
          set_logP( asg->linear_index(), val );
        }
      }
      //ASSERT_TRUE(is_finite());
      return *this;
    }

  public:
    // Currently unused
    //! this(x) = sum_y joint(x,y) * other(y) 
    void convolve(const dense_table_impl& joint,
                         const dense_table_impl& other) {
      // ensure that both tables have the same domain
      DCHECK_EQ(args() + other.args(), joint.args());
      // Initialize the table to zero so we can use it as an accumulator
      uniform(0);
      typename domain_type::const_iterator asg = joint.args().begin();
      typename domain_type::const_iterator end = joint.args().end();
      for( ; asg != end; ++asg) {
        const double value =
            exp(joint.logP(asg->linear_index()) + other.logP(*asg));
        DASSERT_FALSE(std::isinf( value ));
        DASSERT_FALSE(std::isnan( value ));
        // NOTE durring this accumulation, the table is not in log space
        //logP(*asg) += value;
        set_logP(*asg, logP(*asg) + value);
      }

      for(size_t i = 0; i < size(); ++i) {
        double sum = logP(i);
        DCHECK_GE(sum, 0.0);
        if(sum == 0) { set_logP( i, APPROX_LOG_ZERO() ); }
        else { set_logP( i, log(sum) ); }
      }
    }
    

    //! this(x) = other(x, y = asg) 
    void condition(const dense_table_impl& other,
                   const assignment_type& asg) {
      DCHECK_EQ(args(), other.args() - asg.args());
        
      // create a fast assignment starting from the '0' assignment
      // of args() and the conditioning assignment of asg
      fast_discrete_assignment<MAX_DIM> fastyasg(assignment_type(args()) & asg);
      // transpose the remaining assignments to the start
      fastyasg.transpose_to_start(args());
      
      typename domain_type::const_iterator xasg = args().begin(); 
      typename domain_type::const_iterator end = args().end();
      for( ; xasg != end; ++xasg) {
        // REVIEW should this be other.logP(fastasg)? since other and fastyasg 
        // dont have the same domain? (would still need to convert fastasg to 
        // a discrete_assignment)
        set_logP( xasg->linear_index(), other.logP(fastyasg.linear_index()) );
        ++fastyasg;
      }
    }


    //! this(x) = this(x) other(x, y = asg) 
    void times_condition(const dense_table_impl& other,
                         const assignment_type& asg) {
      //assert(args() == other.args() - asg.args());
        
      // create a fast assignment starting from the '0' assignment
      // of args() and the conditioning assignment of asg
      fast_discrete_assignment<MAX_DIM> fastyasg(assignment_type(args()) & asg);
      // transpose the remaining assignments to the start
      fastyasg.transpose_to_start(args());
      if(asg.num_vars() == 0) {
        *this *= other;
      } else {
        typename domain_type::const_iterator xasg = args().begin(); 
        typename domain_type::const_iterator end = args().end();
        for( ; xasg != end; ++xasg) {
          // REVIEW should this be other.logP(fastasg)? since other and fastyasg 
          // dont have the same domain? (would still need to convert fastasg to 
          // a discrete_assignment)
          double val = logP(xasg->linear_index()) + other.logP(fastyasg.linear_index());
          set_logP( xasg->linear_index(), val );
          ++fastyasg;
        }
      }
    }
    
    using table_base_t::marginalize;
    
    //! msg(x) = sum_y this(x,y) 
    void marginalize(dense_table_impl& msg) const {
      // No need to marginalize
      if(args() == msg.args()) {
        // Just copy and return
        msg = *this;
        return;
      }
      // Compute the domain to remove
      domain_type ydom = args() - msg.args();
      DCHECK_GT(ydom.num_vars(), 0);
          
      fast_discrete_assignment<MAX_DIM> fastyasg(args());
      fastyasg.transpose_to_start(ydom);
      // count the number of elements in ydom
      size_t numel = ydom.size();
      // Loop over x
      typename domain_type::const_iterator xasg = msg.args().begin(); 
      typename domain_type::const_iterator end = msg.args().end();
      for( ; xasg != end; ++xasg) {
        double sum = 0;
        for(size_t i = 0;i < numel; ++i) {
          sum += exp(logP(fastyasg.linear_index()));
          ++fastyasg;
        }
        DASSERT_FALSE( std::isinf(sum) );
        DASSERT_FALSE( std::isnan(sum) );
        DCHECK_GE(sum, 0.0);
        if(sum == 0) 
          msg.set_logP( xasg->linear_index(), APPROX_LOG_ZERO() );
        else 
          msg.set_logP( xasg->linear_index(), log(sum) );
      }
    }
      
    using table_base_t::MAP;

    //! msg(x) = max_y this(x,y)
    void MAP(dense_table_impl& msg) const {
      //ASSERT_TRUE(is_finite());
      // No need to marginalize
      if(args() == msg.args()) {
        // Just copy and return
        msg = *this;
        return;
      }
      // Compute the domain to remove
      domain_type ydom = args() - msg.args();
      DCHECK_GT(ydom.num_vars(), 0);
        
      fast_discrete_assignment<MAX_DIM> fastyasg(args());
      fastyasg.transpose_to_start(ydom);
      // count the number of elements in ydom
      size_t numel = ydom.size();
      // Loop over x
      typename domain_type::const_iterator xasg = msg.args().begin();
      typename domain_type::const_iterator end = msg.args().end();
      for( ; xasg != end; ++xasg) {
        double maxval = APPROX_LOG_ZERO();
        for(size_t i = 0;i < numel; ++i) {
          maxval = std::max(maxval, logP(fastyasg.linear_index()));
          ++fastyasg;
        }
        msg.set_logP( xasg->linear_index(), maxval );
      }
      //ASSERT_TRUE(is_finite());
    }

    //! This = other * damping + this * (1-damping) 
    void damp(const dense_table_impl& other, const double& damping) {
      //ASSERT_TRUE(is_finite());
      // This table must be over the same domain as the other
      if(damping == 0) return;
      DCHECK_EQ(args(), other.args());
      DCHECK_GT(damping, 0.0);
      DCHECK_LT(damping, 1.0);
      for(size_t i = 0; i < size(); ++i) {
        double val = damping * exp(other.logP(i)) + (1-damping) * exp(logP(i));
        DCHECK_GE(val, 0);
        if(val == 0) { set_logP( i, APPROX_LOG_ZERO() ); }
        else { set_logP( i, log(val) ); }
        DASSERT_FALSE( std::isinf(logP(i)) );
        DASSERT_FALSE( std::isnan(logP(i)) );
      }
      //ASSERT_TRUE(is_finite());
    }


    //! compute the l_inf norm (i.e. the max diff) between two tables
    double linf_diff(const dense_table_impl& other) const {
      //ASSERT_TRUE(is_finite());
      // This table must be over the same domain as the other
      DCHECK_EQ(args(), other.args());  
      double max_diff = 0;
      for(size_t i = 0; i < size(); ++i) {
        double diff = fabs(exp(other.logP(i)) - exp(logP(i)));
        if (diff > max_diff) {
          max_diff = diff;
        }
      }
      //ASSERT_TRUE(is_finite());
      return max_diff;
    }

    //! compute the average l1 norm between two tables
    double l1_diff(const dense_table_impl& other) const {
      //ASSERT_TRUE(is_finite());
      // This table must be over the same domain as the other
      DCHECK_EQ(args(), other.args());  
      double sum = 0;
      for(size_t i = 0; i < size(); ++i) {
        sum += fabs(exp(other.logP(i)) - exp(logP(i)));
      }
      //ASSERT_TRUE(is_finite());
      return sum / size(); // TODO the l1 norm is not normalized
    }

    //! compute the l1 norm in log space
    double l1_logdiff(const dense_table_impl& other) const {
      DCHECK_EQ(args(), other.args());
      double sum = 0; 
      for(size_t i = 0; i < size(); ++i) {
        sum += fabs(other.logP(i) - logP(i));
      }
      return sum / size(); // TODO the l1 norm is not normalized
    }

    //! argmax(): return the assignment of the largest value
    assignment_type max_asg() const {
      typename domain_type::iterator max_asg = args().begin();
      double max_value = logP(max_asg->linear_index());

      typename domain_type::const_iterator asg = args().begin(); 
      typename domain_type::const_iterator end = args().end();
      for( ; asg != end; ++asg) {
        if(logP(asg->linear_index()) > max_value) {
          max_value = logP(asg->linear_index());
          *max_asg = *asg;
        }
      }
      return *max_asg;
    }

    //! return the linear index of the largest value
    size_t max_index() const {
      return max_asg().linear_index();
    }

    /**
     * Compute the expectation of the dense table
     */
    inline void expectation(std::vector<double>& values) const {
      values.clear();
      values.resize(num_vars(), 0);
      double sum = 0;
      typename domain_type::const_iterator asg = args().begin(); 
      typename domain_type::const_iterator end = args().end();
      for( ; asg != end; ++asg) {
      const double scale = exp(logP(asg->linear_index()));
        sum += scale;
        typename assignment_type::const_iterator asg_it = asg.begin();
        for(size_t i = 0; i < num_vars(); ++i) {
          values[i] += asg_it[i] * scale;
        }
      }
      // Rescale for normalization
      for(size_t i = 0; i < num_vars(); ++i)  values[i] /= sum;
    } // end of expectation

    /**
     * Draw a sample from the dense table
     */
    inline assignment_type sample() const {
      DCHECK_GT(size(), 0);
      // This table must be normalized
      const double t = graphlab::random::rand01();
      DCHECK_GE( t, 0 );
      DCHECK_LT( t, 1 );
      double sum = 0;
      for(size_t i = 0; i < size(); ++i) {
        sum += exp( logP(i) );
        if(t <=sum) { return assignment_type(args(), i); }
        DCHECK_LT(sum, 1);
      }
      // Unreachable
      throw("Invalid state reached in sample()");
      return assignment_type();
    } // end of sample
    

    /**
     * Construct a binary agreement factor
     */
    void set_as_agreement(const double& lambda) {
      DCHECK_EQ(num_vars(), 2);

      typename domain_type::const_iterator asg = args().begin();
      typename domain_type::const_iterator end = args().end();
      for( ; asg != end; ++asg) {
        typename assignment_type::const_iterator asg_it = asg->begin();
        const int diff = abs( int(asg_it[0]) - int(asg_it[1]) );
        if( diff > 0) { set_logP( asg->linear_index(), -lambda ); }
        else { set_logP( asg->linear_index(), 0 ); }
      }
    } // end of set_as_agreement
    
    void set_as_laplace(const double& lambda) {
      DCHECK_EQ(num_vars(), 2);
      typename domain_type::const_iterator asg = args().begin();
      typename domain_type::const_iterator end = args().end();
      for( ; asg != end; ++asg) {
        typename assignment_type::const_iterator asg_it = asg->begin();
        const int diff = abs( int(asg_it[0]) - int(asg_it[1]) );
        set_logP( asg->linear_index(), -diff * lambda );
      }
    } // end of set_as_laplace

  public:
    void load(graphlab::iarchive& arc) {
      arc >> _args;
      arc >> _data;
    }

    void save(graphlab::oarchive& arc) const {
      arc << _args;
      arc << _data;
    }
    
  private:
  //! The domain of the table (the arity of the table, along with its cardinality) 
    domain_type _args;
    // NOTE _data is ordered according to the linear indexing of our domain, 
    // which is ordered such that the variable with the lowest id iterates 
    // fastest. weird! i know!
    std::vector<double> _data;
  }; // End of dense table



// REVIEW could move the these methods into dense_table_impl (which 
// would become dense_table again) like in sparse_table 
  template<size_t MAX_DIM> 
  class dense_table : public dense_table_impl<MAX_DIM> {
    typedef dense_table_impl<MAX_DIM>  dense_table_impl_t;

    // "using typename" parses but doesn't work in GCC < 4.7. Use "typedef typename" instead.
    //using typename dense_table_impl_t::table_base_t;
    //using typename dense_table_impl_t::variable_type;
    //using typename dense_table_impl_t::domain_type;
    //using typename dense_table_impl_t::assignment_type;
    typedef typename dense_table_impl_t::table_base_t     table_base_t;
    typedef typename dense_table_impl_t::variable_type    variable_type;
    typedef typename dense_table_impl_t::domain_type      domain_type;
    typedef typename dense_table_impl_t::assignment_type  assignment_type;

  public:
    /** Construct an empty dense table */
    dense_table() : dense_table_impl_t() { }
    
    /** Construct a dense table over the given domain */
    dense_table(const domain_type& dom) :
        dense_table_impl_t(dom) { }

    /** Construct a dense table over the given domain and distribution */
    dense_table(const domain_type& dom, const std::vector<double> &data) :
        dense_table_impl_t(dom, data) { }

    /** Construct a dense table over the given domain */
    dense_table(const std::vector<variable_type> &args) : 
        dense_table_impl_t(args) { }

    /** Construct a dense table over the given domain */
    dense_table(const std::vector<variable_type> &args, 
        const std::vector<double> &data) : dense_table_impl_t(args, data) { }

    // REVIEW make these static factory methods: 
    // static dense_table<MAX_DIM>& unary_table_node(...) or use the virtual
    // constructor idiom or something
    /** Construct a unary dense table over the given var */
    dense_table(const variable_type &var) : dense_table_impl_t(var) { }

    /** Construct a unary dense table over the given var */
    dense_table(const variable_type& var, std::vector<double>& logd) : 
        dense_table_impl_t(var, logd) { }

    /** Construct a unary dense table over the given var */
    dense_table(const variable_type& var, 
        double const* const begin, double const* const end) :
        dense_table_impl_t(var, begin, end) { }

    /** Construct a copy */
    dense_table(const dense_table& other) :
        dense_table_impl_t(other) { }

    virtual ~dense_table() { }

    // REVIEW currently, this isnt necessary
    /** Standard assignment operator */
    dense_table& operator=(const dense_table& other) {
      if(this == &other) 
        return *this;

      dense_table_impl_t::operator=(other);

      return *this;
    }

    friend std::ostream& operator<<(std::ostream& out,
                              const dense_table<MAX_DIM>& tbl) {
      out << "Dense Table: " << tbl.args() << "{" << std::endl;
      typename domain_type::const_iterator asg = tbl.args().begin();
      typename domain_type::const_iterator end = tbl.args().end();
      for( ; asg != end; ++asg) {
        out << "\tLogP(" << *asg << ")=" << tbl.logP(*asg) << std::endl;
      }
      out << "}";

//      dense_table<MAX_DIM>::const_iterator asg = tbl.args().begin();
//      dense_table<MAX_DIM>::const_iterator end = tbl.args().end();
//      for( ; asg != end; ++asg) {
//        out << tbl.logP(*asg) << " ";
//      }

      return out;
    }

  // virtual methods
  public:
    typedef dense_table const *const const_ptr;

    virtual dense_table& deep_copy(const table_base_t& base) {
      if(this == &base) return *this;

      // ensure we are dealing with a dense_table
      const_ptr other = dynamic_cast<const_ptr>(&base);
      if(other == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      *this = *other;
      return *this;
    }

    using dense_table_impl_t::copy_onto;

    virtual dense_table& copy_onto(const table_base_t& base) {
      if(this == &base) return *this;

      // ensure we are dealing with a dense_table
      const_ptr other = dynamic_cast<const_ptr>(&base);
      if(other == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      dense_table_impl_t::copy_onto(*other);
      return *this;
    }

    virtual const variable_type& var(const size_t index) const {
      return dense_table_impl_t::domain().var(index);
    }
/*
    //! this(x) += other(x);
    virtual dense_table& plus_equals(const table_base_t& base) {
      // ensure we are dealing with a dense_table
      const_ptr other = dynamic_cast<const_ptr>(&base);
      if(other == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      *this += *other;
      return *this;
    }
*/
    //! this(x) *= other(x);
    virtual dense_table& times_equals(const table_base_t& base) {
      // ensure we are dealing with a dense_table
      const_ptr other = dynamic_cast<const_ptr>(&base);
      if(other == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      *this *= *other;

      return *this;
    }

    //! this(x) /= other(x);
    virtual dense_table& divide_equals(const table_base_t& base) {
      // ensure we are dealing with a dense_table
      const_ptr other = dynamic_cast<const_ptr>(&base);
      if(other == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      *this /= *other;

      return *this;
    }

    //! (out(x) = this(x)) * other(x);
    virtual void times(const table_base_t& base, 
        table_base_t& out_base) const {

      // ensure we are dealing with a dense_table
      dense_table *const out = 
          dynamic_cast<dense_table *const>(&out_base);
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

      // ensure we are dealing with a dense_table
      dense_table *const out = 
          dynamic_cast<dense_table *const>(&out_base);
      if(out == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      *out = *this; // deep copy
      out->divide_equals(base);
    }

    using dense_table_impl_t::marginalize;

    virtual void marginalize(table_base_t& base) const {
      // ensure we are dealing with a dense_table
      dense_table* msg = dynamic_cast<dense_table*>(&base);
      if(msg == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      dense_table_impl_t::marginalize(*msg);
    }

    using dense_table_impl_t::MAP;

    virtual void MAP(table_base_t& base) const {
      // ensure we are dealing with a dense_table
      dense_table* msg = dynamic_cast<dense_table*>(&base);
      if(msg == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      dense_table_impl_t::MAP(*msg);
    }
    
  public:
    virtual void load(graphlab::iarchive& arc) { 
      dense_table_impl_t::load(arc); 
    }

    virtual void save(graphlab::oarchive& arc) const {
      dense_table_impl_t::save(arc); 
    }

    virtual std::ostream& print(std::ostream& out = std::cout) const {
      // ensure we are dealing with a dense_table
      const_ptr tbl = dynamic_cast<const_ptr>(this);
      if(tbl == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      out << *tbl;
      return out;
    }
  };


}; // end of namespace graphlab



#include <graphlab/macros_undef.hpp>
#endif // DENSE_TABLE_HPP

