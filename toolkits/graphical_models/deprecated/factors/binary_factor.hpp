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


#ifndef BINARY_FACTOR_HPP
#define BINARY_FACTOR_HPP


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
#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/serialization/serialization_includes.hpp>

// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>



/**
 * A binary factor is a table over a pair of variables and is
 * assocaited with each edge in a pairwise markov random field.  All
 * data is represented in log form.
 */
class binary_factor {
  
public:
    
  binary_factor(uint32_t var1 = 0,
                uint16_t arity1 = 0,
                uint32_t var2 = 0,
                uint16_t arity2 = 0) :
    _var1(var1), _arity1(arity1), _var2(var2), _arity2(arity2),
    _data(arity1 * arity2) { }

  binary_factor(const binary_factor& other) :
    _var1(other._var1), _arity1(other._arity1),
    _var2(other._var2), _arity2(other._arity2),
    _data(other._data) { }

  binary_factor& operator=(const binary_factor& other) {
    _var1 = other._var1;
    _arity1 = other._arity1;
    _var2 = other._var2;
    _arity2 = other._arity2;
    _data = other._data;
    return *this;
  }

  void resize(uint16_t arity1, uint16_t arity2) {
    _arity1 = arity1;
    _arity2 = arity2;
    _data.resize(_arity1 * _arity2);
  }
  
  uint32_t& var1() { return _var1;  } 
  const uint32_t& var1() const { return _var1; }
  uint32_t& var2() { return _var2;  } 
  const uint32_t& var2() const { return _var2; } 

  const uint16_t& arity1() const { return _arity1; }
  const uint16_t& arity2() const { return _arity2; } 

  /** Get the value of the factor.  In var1 == var2 the variables
      are ignored. */
  double& logP(uint32_t x1, uint16_t asg1,
               uint32_t x2, uint16_t asg2) {
    // If the factor is not symmetric then we may have to match the
    // arguments
    if( _var1 != _var2 ) {
      assert((x1 == var1() && x2 == var2()) ||
             (x2 == var1() && x1 == var2()));
      if(x1 == var2() && x2 == var1()) std::swap(asg1, asg2);
    }
    assert( asg1 < arity1() );
    assert( asg2 < arity2() );
    // return value
    return _data[asg1 + asg2 * arity1()];
  } // end of logP for a binary factor


  const double& logP(uint32_t x1, uint16_t asg1,
                     uint32_t x2, uint16_t asg2) const {
    // If the factor is not symmetric then we may have to match the
    // arguments
    if( _var1 != _var2 ) {
      assert((x1 == var1() && x2 == var2()) ||
             (x2 == var1() && x1 == var2()));
      if(x1 == var2() && x2 == var1()) std::swap(asg1, asg2);
    }
    ASSERT_LT( asg1 , arity1() );
    ASSERT_LT( asg2 , arity2() );
    // return value
    return _data[asg1 + asg2 * arity1()];
  } // end of logP for a binary factor

  
  double& logP(uint16_t asg1, uint16_t asg2) {
    ASSERT_LT( asg1 , arity1() );
    ASSERT_LT( asg2 , arity2() );
    return _data[asg1 + asg2 * arity1()];
  } // end of logP for a binary factor


  const double& logP(uint16_t asg1, uint16_t asg2) const {
    assert( asg1 < arity1() );
    assert( asg2 < arity2() );
    return _data[asg1 + asg2 * arity1()];
  } // end of logP for a binary factor

    /** ensure that sum_x this(x) = 1 */
  inline void normalize() {
    assert(arity1() > 0);
    assert(arity2() > 0);
    // Compute the max value
    double max_value = logP(0,0);
    for(uint16_t asg1 = 0; asg1 < arity1(); ++asg1) 
      for(uint16_t asg2 = 0; asg2 < arity2(); ++asg2)  
        max_value = std::max(max_value, 
                             logP(asg1, asg2));
    assert( !std::isinf(max_value) );
    assert( !std::isnan(max_value) );
    // scale and compute normalizing constant
    double Z = 0.0;
    for(uint16_t asg1 = 0; asg1 < arity1(); ++asg1) 
      for(uint16_t asg2 = 0; asg2 < arity2(); ++asg2)  
        Z += std::exp(logP(asg1, asg2) -= max_value);
    assert( !std::isinf(Z) );
    assert( !std::isnan(Z) );
    assert( Z > 0.0);
    double logZ = std::log(Z);
    // Normalize
    for(uint16_t asg1 = 0; asg1 < arity1(); ++asg1) 
      for(uint16_t asg2 = 0; asg2 < arity2(); ++asg2)  
        logP(asg1, asg2) -= logZ;
  } // End of normalize

  
  
  void set_as_agreement(double lambda) {
    for(uint16_t i = 0; i < arity1(); ++i) { 
      for(uint16_t j = 0; j < arity2(); ++j) { 
        if( i != j) logP(i,j) = -lambda;
        else logP(i,j) = 0;
      }
    }
  } // end of set_as_agreement
  
  void set_as_laplace(double lambda) {
    for(uint16_t i = 0; i < arity1(); ++i) { 
      for(uint16_t j = 0; j < arity2(); ++j) { 
        logP(i,j) = -std::abs(double(i) - double(j)) * lambda;
      }
    }
  } // end of set_as_laplace



    /**
     * Compute the Mooji Kappen Message derivative. 
     */
  double mk_derivative() const {
    double max_value = -std::numeric_limits<double>::max();
    for(uint16_t a = 0; a < arity1(); ++a) {
      for(uint16_t b = 0; b < arity2(); ++b) {
        for(uint16_t x = 0; x < arity1(); ++x) {
          for(uint16_t y = 0; y < arity2(); ++y) {
            if(a != x && b != y) {
              double value =
                (logP(a,b) + logP(x,y) - (logP(x,b) + logP(a,y)))/4.0;
              value = std::tanh(value);
              max_value = std::max(max_value, value);
            }
          }
        }
      }
    }
    return max_value;
  }

  //! Compute the Ihler dynamic range
  double ihler_dynamic_range() const {
    double min_value = *std::min_element(_data.begin(), _data.end());
    double max_value = *std::max_element(_data.begin(), _data.end());
    return std::exp((max_value - min_value)/2);
  }

    

  //! Print the factor description
  void printP(std::ostream& out) const {
    out << "Binary Factor(v_" << var1() << " in {1..."
        << arity1() << "}, " 
        << ", v_ " << var2() << " in {1..." 
        << arity2() << "})" << std::endl;
    for(uint16_t i = 0; i < arity1(); ++i) {
      for(uint16_t j = 0; j < arity2(); ++j) {
        out << std::exp(logP(i,j)) << " ";
      }
      out << std::endl;
    }
  }

  //! Save the factor to a file
  void save(graphlab::oarchive &oarc) const {
    oarc << _var1 << _arity1 
         << _var2 << _arity2
         << _data;
  }

  //! Load the factor from a file
  void load(graphlab::iarchive &iarc) {
    iarc >> _var1 >> _arity1 
         >> _var2 >> _arity2
         >> _data;
  }


    
private:
  uint32_t _var1;
  uint16_t _arity1;
  uint32_t _var2;
  uint16_t _arity2;
  std::vector<double> _data;
    
}; // end of class binary_factor



inline std::ostream& operator<<(std::ostream& out, 
                                const binary_factor& fact) {
  out << "Binary Factor(v_" << fact.var1() << " in {1..."
      << fact.arity1() << "}, " 
      << ", v_ " << fact.var2() << " in {1..." 
      << fact.arity2() << "})" << std::endl;
  for(uint16_t i = 0; i < fact.arity1(); ++i) {
    for(uint16_t j = 0; j < fact.arity2(); ++j) {
      out << fact.logP(i,j) << " ";
    }
    out << std::endl;
  }
  return out;
} // end of operator<<


#include <graphlab/macros_undef.hpp>
#endif

