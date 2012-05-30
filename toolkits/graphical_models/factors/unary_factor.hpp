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


#ifndef UNARY_FACTOR_HPP
#define UNARY_FACTOR_HPP


/**
 * This file contains the definitions of some of the basic factor
 * types needed for loopy belief propagation.  This is demo code and
 * is intentionally kept as simple as possible.
 *
 *  \author Joseph Gonzalez
 */


// Including Standard Libraries
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

// Random number generation
#include <graphlab/util/random.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

// Basic unary factor
#include "binary_factor.hpp"


// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>



  
/**
 * A unary factor is a table over a single variable and is associated
 * with edge variable in the pairwise markov random field.  Unary
 * factors are also used to represent messages. All data is
 * represented in log form.
 */
class unary_factor {  
public:

  /**
   * Construct a unary factor with a given variable and arrity.
   */
  unary_factor(uint32_t var = 0, uint16_t arity = 0) :
    _var(var), _data(arity) {}

  // /**
  //  * Copy constructor for a unary factor
  //  */
  // unary_factor(const unary_factor& other) :
  //   _var(other._var), _data(other._data) { }

  // unary_factor& operator=(const unary_factor& other) {
  //   _var = other._var;
  //   _data = other._data;
  //   return *this;
  // }

  /**
   * Increase the size of the factor
   */
  inline void resize(uint16_t arity) {
    _data.resize(arity);
  }
        
  uint32_t& var() { return _var;  }
  const uint32_t& var() const { return _var; }
  uint16_t arity() const { return (uint16_t)_data.size(); }
  size_t size() const { return _data.size(); }
  bool empty() const { return _data.empty(); }


  inline double& logP(size_t asg) {
    assert(asg < arity()); return _data[asg];
  } // end of logP for a unary factor

  inline const double& logP(size_t asg) const {
    assert(asg < arity()); return _data[asg];
  } // end of logP for a unary factor 

    /** zero out the factor */
    // inline void zero() {
    //   for(size_t asg = 0; asg < arity(); ++asg) logP(asg) = 0;
    // }

  inline void uniform(double value = 0) {
    for(uint16_t asg = 0; asg < arity(); ++asg) logP(asg) = value;
  }


  /** ensure that sum_x this(x) = 1 */
  inline void normalize() {
    assert(arity() > 0);
    // Compute the max value
    double max_value = logP(0);
    for(uint16_t asg = 0; asg < arity(); ++asg) 
      max_value = std::max(max_value, logP(asg));
    assert( !std::isinf(max_value) );
    assert( !std::isnan(max_value) );
    // scale and compute normalizing constant
    double Z = 0.0;
    for(uint16_t asg = 0; asg < arity(); ++asg) 
      Z += std::exp(logP(asg) -= max_value);
    assert( !std::isinf(Z) );
    assert( !std::isnan(Z) );
    assert( Z > 0.0);
    double logZ = std::log(Z);
    // Normalize
    for(uint16_t asg = 0; asg < arity(); ++asg)
      logP(asg) -= logZ;
  } // End of normalize

  inline unary_factor& operator*=(const unary_factor& other) {
    ASSERT_EQ(arity(), other.arity());
    for(uint16_t asg = 0; asg < arity(); ++asg)
      logP(asg) += other.logP(asg);
    return *this;
  }

    /** this(x) += other(x); */
  inline unary_factor& operator+=(const unary_factor& other) {
    ASSERT_EQ(arity(), other.arity());
    for(uint16_t asg = 0; asg < arity(); ++asg)
      logP(asg) = log(exp(logP(asg)) + exp(other.logP(asg)));
    return *this;
  } // end plus


    /** this(x) /= other(x); */
  inline unary_factor& operator/=(const unary_factor& other) {
    ASSERT_EQ(arity(), other.arity());
    for(uint16_t asg = 0; asg < arity(); ++asg)
      logP(asg) -= other.logP(asg);
    return *this;
  } // end of divide
  
    /** this(x) = sum_y fact(x,y) * other(y) */
  inline void convolve(const binary_factor& bin_fact,
                       const unary_factor& other) {
    // Compute C(x) = Sum_y A(x,y) B(y)
    for(uint16_t x = 0; x < arity(); ++x) {
      double sum = 0.0;
      for(uint16_t y = 0; y < other.arity(); ++y) {          
        sum += std::exp(bin_fact.logP(var(), x, other.var(), y) +
                        other.logP(y));
      }
      assert( !(sum < 0.0) );
      // Gaurd against zeros
      if(sum == 0) sum = std::numeric_limits<double>::min();
      logP(x) = std::log(sum);
    }
  }
  
  /** this(x) = this(x) * fact(x, asg) */
  inline void condition(const binary_factor& bin_fact,
                        uint16_t asg) {
    uint32_t other_var =
      var() != bin_fact.var1()? bin_fact.var1() : bin_fact.var2();
    for(uint16_t x = 0; x < arity(); ++x) 
      logP(x) += bin_fact.logP(var(), x, other_var, asg);    
  } // end of condition
  
  
    /** This = other * damping + this * (1-damping) */
  inline void damp(const unary_factor& other, double damping) {
    assert(arity() == other.arity());
    if(damping == 0) return;
    assert(damping >= 0.0);
    assert(damping < 1.0);
    for(uint16_t asg = 0; asg < arity(); ++asg) 
      logP(asg) = std::log(damping * std::exp(other.logP(asg)) + 
                           (1.0 - damping) * std::exp(logP(asg)));  
  } // end of damp
  
  
    /** Compute the residual between two unary factors */
  inline double residual(const unary_factor& other) const {  
    assert(arity() == other.arity());
    double residual = 0;
    for(uint16_t asg = 0; asg < arity(); ++asg) 
      residual += std::abs(std::exp(logP(asg)) - 
                           std::exp(other.logP(asg)));
    return residual / arity();
  } // end of residual
  
  
    /** get the max assignment*/
  inline size_t max_asg() const {  
    size_t max_asg = 0;
    double max_value = logP(0);
    for(uint16_t asg = 0; asg < arity(); ++asg) { 
      if(logP(asg) > max_value) {
        max_value = logP(asg);
        max_asg = asg;
      }
    }
    return max_asg;
  } // end of max asg
  
    /** Get the expected assignment */
  inline double expectation() const {  
    double sum = 0;
    double s2 = 0;
    for(uint16_t asg = 0; asg < arity(); ++asg)  {
      sum += asg * std::exp(logP(asg));       
      s2 += std::exp(logP(asg));       
    }
    return sum / s2;;
  } // end of expectation
  
    /** Draw a random sample from the factor */
  inline size_t sample() const {  
    // Using the cdf method to generate a random sample
    assert(arity() > 0);
    // double t = static_cast<double>(rand()) / RAND_MAX;
    double t = graphlab::random::rand01();
    assert( t >= 0);
    assert(t < 1);
    double sum = 0.0;
    for(uint16_t asg = 0; asg < arity(); ++asg) {
      sum += exp(logP(asg));
      if(t <= sum) return asg;
      assert(sum < 1);
    }
    // We were unable to draw a sample;
    assert(false);
  } // end of sample

  void save(graphlab::oarchive &oarc) const {
    oarc << _var << _data;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> _var >> _data;
  }
private:
  uint32_t _var;
  std::vector<double> _data;
};  // End of unary factor



inline std::ostream& operator<<(std::ostream& out, 
                                const unary_factor& fact) {
  out << "Unary Factor(" << fact.arity()  << ")"
      << std::endl;
  for(size_t i = 0; i < fact.arity(); ++i) {
    out << fact.logP(i) << " ";
  }
  out << std::endl;
  return out;
} // end of operator<<


#include <graphlab/macros_undef.hpp>
#endif

