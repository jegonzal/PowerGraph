/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_FACTORS_HPP
#define GRAPHLAB_FACTORS_HPP


/**
 * This file contains the definitions of some of the basic factor
 * types needed for loopy belief propagation.  This is demo code and
 * is intentionally kept as simple as possible.
 *
 *  \author Joseph Gonzalez
 */



// INCLUDES ===================================================================>

// Including Standard Libraries
#include <cassert>
#include <cmath>

#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <set>


// Random number generation
#include <graphlab/util/random.hpp>

#include <graphlab/serialization/serialization_includes.hpp>

#include <graphlab/factors/discrete_variable.hpp>
#include <graphlab/factors/discrete_domain.hpp>
#include <graphlab/factors/discrete_assignment.hpp>
#include <graphlab/factors/fast_discrete_assignment.hpp>



// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

namespace graphlab {





  // static const double EPSILON = std::numeric_limits<double>::min();
  // static const double LOG_EPSILON = log(EPSILON);
  // static const double MAX_DOUBLE =  std::numeric_limits<double>::max();
  //  static const double LOG_MAX_DOUBLE = std::log(MAX_DOUBLE);




  
  /**
   * A factor over up to max_dim dimensions
   */
  template<size_t MAX_DIM>
  class table_factor {
  public:
    typedef discrete_variable            variable_type;
    typedef discrete_domain<MAX_DIM>     domain_type;
    typedef discrete_assignment<MAX_DIM> assignment_type;
    

    /** Construct an empty table factor */
    table_factor() { }
    
    /** Construct a table factor over the given domain */
    table_factor(const domain_type& dom) :
      _args(dom), _data(dom.size()) {  }

    /** Construct a copy */
    table_factor(const table_factor& other) :
      _args(other._args), _data(other._data) { }

    /** Standard assignment operator */
    table_factor& operator=(const table_factor& other) {
      _args = other._args;
      _data = other._data;
      return *this;
    }

    void set_args(const domain_type& args) {
      _args = args;
      _data.resize(args.size());
    }
    
    const domain_type& args() const {
      return _args;
    }

    const double& logP(size_t index) const {
      ASSERT_LT(index, _data.size());
      return _data[index];
    }
    
    double& logP(size_t index) {
      ASSERT_LT(index, _data.size());
      return _data[index];
    }

    const double& logP(const assignment_type& asg) const {
      if(asg.args() == args()) {
        // if the assignment index matches
        const size_t index(asg.linear_index());
        ASSERT_LT(index, _data.size());
        return _data[index];
      } else {
        // Restrict the assignment to this domain
        const assignment_type sub_asg = asg.restrict(_args);
        const size_t index(sub_asg.linear_index());
        ASSERT_LT(index, _data.size());
        return _data[index];
      }
    }
    

    double& logP(const assignment_type& asg) {
      if(asg.args() == args()) {
        // if the assignment index matches
        const size_t index(asg.linear_index());
        ASSERT_LT(index, _data.size());
        return _data[index];
      } else {
        // Restrict the assignment to this domain
        const assignment_type sub_asg = asg.restrict(_args);
        const size_t index(sub_asg.linear_index());
        ASSERT_LT(index, _data.size());
        return _data[sub_asg.linear_index()];
      }
    } // end of logP

    size_t size() const { return _args.size(); }

    size_t num_vars() const { return _args.num_vars(); }

    double* begin() { return _data.begin(); }
    const double* begin() const { return _data.begin(); }
    

    
    double* end() { return _data.end(); }
    const double* end() const { return _data.end(); }
    
    void zero() { std::fill(_data.begin(), _data.end(), 0); }
        
    void uniform() {
      std::fill(_data.begin(), _data.end(),
                log(1.0/_data.size()));
    }
 
    void uniform(double value) {
      std::fill(_data.begin(), _data.end(), value);
    } 

    
    //! ensure that sum_x this(x) = 1 
    void normalize() {
      // Compute the max value
      double max_value = logP(0);
      for(size_t i = 0; i < _data.size(); ++i) 
        max_value = std::max(max_value, _data[i] );
      // scale and compute normalizing constant
      double Z = 0.0;
      for(size_t i = 0; i < _data.size(); ++i)
        Z += exp( _data[i] -= max_value );
      // assert( !std::isinf(Z) );
      // assert( !std::isnan(Z) );
      // assert( Z > 0.0);
      const double logZ(log(Z));
      ASSERT_FALSE( std::isinf(logZ) );
      ASSERT_FALSE( std::isnan(logZ) );
      // Normalize
      for(size_t i = 0; i < _data.size(); ++i) _data[i] -= logZ;
    } // End of normalize
    

    /** 
     * Ensure that the largest value in log form is zero.  This
     * prevents overflows on normalization. 
     */
    void shift_normalize() {
      // Compute the max value
      double max_value = logP(0);
      for(size_t i = 0; i < _data.size(); ++i) 
        max_value = std::max(max_value, _data[i]);
      for(size_t i = 0; i < _data.size(); ++i) _data[i] -= max_value;
    }


    /**
     * Return false if any of the entries are not finite 
     */
    bool is_finite() { 
      for(size_t i = 0; i < _data.size(); ++i) {
        const bool is_inf( std::isinf( _data[i] ) );
        const bool is_nan( std::isnan( _data[i] ) );
        if( __builtin_expect( is_inf || is_nan, 0) ) return false;
      }
      return true;
    }


    //! this(x) += other(x);
    inline table_factor& operator+=(const table_factor& other) {
      if(args() == other.args()) {
        ASSERT_EQ(_data.size(), other._data.size());
        // More verctorizable version
        for(size_t i = 0; i < _data.size(); ++i) 
          _data[i] = log( exp(_data[i]) + exp(other._data[i]) );
      } else { 
        for(assignment_type asg = args().begin(); asg < args().end(); ++asg) {
          logP(asg.linear_index()) =
            log( exp(logP(asg.linear_index())) + exp(other.logP(asg)) );
          // ASSERT_FALSE(std::isinf( logP(asg.linear_index()) ));
          // ASSERT_FALSE(std::isnan( logP(asg.linear_index()) ));
        }
      }
      return *this;
    }

    

    //! this(x) *= other(x);
    inline table_factor& operator*=(const table_factor& other) {
      if(args() == other.args()) {
        ASSERT_EQ(_data.size(), other._data.size());
        // More verctorizable version
        for(size_t i = 0; i < _data.size(); ++i) _data[i] += other._data[i];
      } else {
        for(assignment_type asg = args().begin(); asg < args().end(); ++asg) {
          logP(asg.linear_index()) += other.logP(asg); 
          // ASSERT_FALSE(std::isinf( logP(asg.linear_index()) ));
          // ASSERT_FALSE(std::isnan( logP(asg.linear_index()) ));
        }
      }
      return *this;
    }

    //! Create a table factor on the fly
    table_factor operator*(const table_factor& other) const {
      table_factor factor = *this;
      return factor *= other;
    }

    //! this(x) /= other(x);
    inline table_factor& operator/=(const table_factor& other) {
      if(args() == other.args()) {
        ASSERT_EQ(_data.size(), other._data.size());
         // More verctorizable version
        for(size_t i = 0; i < _data.size(); ++i) _data[i] -= other._data[i];
      } else { 
        for(assignment_type asg = args().begin(); asg < args().end(); ++asg) {
          logP(asg.linear_index()) -= other.logP(asg);
          // ASSERT_FALSE(std::isinf( logP(asg.linear_index()) ));
          // ASSERT_FALSE(std::isnan( logP(asg.linear_index()) ));
        }
      }
      return *this;
    }

    //! Create a table factor on the fly
    table_factor operator/(const table_factor& other) const {
      table_factor factor = *this;
      return factor /= other;
    }

    
    // Currently unused
    //! this(x) = sum_y joint(x,y) * other(y) 
    inline void convolve(const table_factor& joint,
                         const table_factor& other) {
      // ensure that both factors have the same domain
      ASSERT_EQ(args() + other.args(), joint.args());
      // Initialize the factor to zero so we can use it as an
      // accumulator
      uniform(0);
      for(assignment_type asg = joint.args().begin();
          asg < joint.args().end(); ++asg) {
        const double value =
          exp(joint.logP(asg.linear_index()) + other.logP(asg));
        ASSERT_FALSE(std::isinf( value ));
        ASSERT_FALSE(std::isnan( value ));
        logP(asg) += value;
      }

      for(size_t i = 0; i < _data.size(); ++i) {
        double& sum = _data[i];
        ASSERT_GE(sum, 0.0);
        if(sum == 0) sum = -std::numeric_limits<double>::max();
        else sum = log(sum);
      }
      
      // // ensure that both factors have the same domain
      // assert(args() + other.args() == joint.args());
      // for(assignment_type xasg = args().begin(); 
      //     xasg < args().end(); ++xasg) {
      //   double sum = 0;
      //   for(assignment_type yasg = other.args().begin(); 
      //       yasg < other.args().end(); ++yasg) {
      //     assignment_type joint_asg = xasg & yasg;
      //     sum += std::exp(joint.logP(joint_asg.linear_index()) + 
      //                     other.logP(yasg.linear_index()));
      //   }
      //   assert( !std::isinf(sum) );
      //   assert( !std::isnan(sum) );
      //   assert(sum >= 0.0);
      //   if(sum == 0) logP(xasg.linear_index()) = -std::numeric_limits<double>::max();
      //   else logP(xasg.linear_index()) = std::log(sum);
      // }
    }
    

    //! this(x) = other(x, y = asg) 
//     inline void condition(const table_factor& other,
//                           const assignment_type& asg) {
//       assert(args() == other.args() - asg.args());
//       for(assignment_type xasg = args().begin(); 
//           xasg < args().end(); ++xasg) {
//         assignment_type joint = xasg & asg;
//         assert(joint.args() == other.args());
//         logP(xasg.linear_index()) = other.logP(joint.linear_index());        
//       }
//     }




    //! this(x) = other(x, y = asg) 
    inline void condition(const table_factor& other,
                          const assignment_type& asg) {
      ASSERT_EQ(args(), other.args() - asg.args());
      
      // create a fast assignment starting from the '0' assignment
      // of args() and the conditioning assignment of asg
      fast_discrete_assignment<MAX_DIM> fastyasg(args().begin() & asg);
      // transpose the remaining assignments to the start
      fastyasg.transpose_to_start(args());
      
      for(assignment_type xasg = args().begin(); 
          xasg < args().end(); ++xasg) {
        logP(xasg.linear_index()) = other.logP(fastyasg.linear_index());        
        ++fastyasg;
      }
    }



//     inline void times_condition(const table_factor& other,
//                                 const assignment_type& asg) {
//       //      assert(args() == other.args() - asg.args());
//       if(asg.num_vars() == 0) {
//         *this *= other;
//       } else {
//         for(assignment_type xasg = args().begin(); 
//             xasg < args().end(); ++xasg) {
//           assignment_type joint = xasg & asg;
//           assert(joint.args() == other.args());
//           logP(xasg.linear_index()) += other.logP(joint.linear_index());        
//         }
//       }
//     }


    //! this(x) = this(x) other(x, y = asg) 
    inline void times_condition(const table_factor& other,
                                const assignment_type& asg) {
      //      assert(args() == other.args() - asg.args());
      
      // create a fast assignment starting from the '0' assignment
      // of args() and the conditioning assignment of asg
      fast_discrete_assignment<MAX_DIM> fastyasg(args().begin() & asg);
      // transpose the remaining assignments to the start
      fastyasg.transpose_to_start(args());
      if(asg.num_vars() == 0) {
        *this *= other;
      } else {
        for(assignment_type xasg = args().begin(); 
            xasg < args().end(); ++xasg) {
          logP(xasg.linear_index()) += other.logP(fastyasg.linear_index());        
          ++fastyasg;
        }
      }
    }
   
    //! this(x) = sum_y joint(x,y) 
//     inline void marginalize(const table_factor& joint) {
//       // No need to marginalize
//       if(args() == joint.args()) {
//         // Just copy and return
//         *this = joint;
//         return;
//       }
//       // Compute the domain to remove
//       domain_type ydom = joint.args() - args();
//       assert(ydom.num_vars() > 0);
//       
//       // Loop over x
//       for(assignment_type xasg = args().begin(); 
//           xasg < args().end(); ++xasg) {
//         double sum = 0;
//         for(assignment_type yasg = ydom.begin(); 
//             yasg < ydom.end(); ++yasg) {
//           assignment_type joint_asg = xasg & yasg;
//           sum += exp(joint.logP(joint_asg.linear_index()));
//         }
//         assert( !std::isinf(sum) );
//         assert( !std::isnan(sum) );
//         assert(sum >= 0.0);
//         if(sum == 0) 
//           logP(xasg.linear_index()) = -std::numeric_limits<double>::max();
//         else logP(xasg.linear_index()) = log(sum);
//       }
//     }

   
    //! this(x) = sum_y joint(x,y) 
    inline void marginalize(const table_factor& joint) {
      // No need to marginalize
      if(args() == joint.args()) {
        // Just copy and return
        *this = joint;
        return;
      }
      // Compute the domain to remove
      domain_type ydom = joint.args() - args();
      ASSERT_GT(ydom.num_vars(), 0);
      
      fast_discrete_assignment<MAX_DIM> fastyasg(joint.args());
      fastyasg.transpose_to_start(ydom);
      // count the number of elements in ydom
      size_t numel = 1;
      for (size_t i = 0;i < ydom.num_vars(); ++i) {
        numel *= ydom.var(i).size();
      }
      // Loop over x
      for(assignment_type xasg = args().begin(); 
          xasg < args().end(); ++xasg) {
        double sum = 0;
        for(size_t i = 0;i < numel; ++i) {
          sum += exp(joint.logP(fastyasg.linear_index()));
          ++fastyasg;
        }
        ASSERT_FALSE( std::isinf(sum) );
        ASSERT_FALSE( std::isnan(sum) );
        ASSERT_GE(sum, 0.0);
        if(sum == 0) 
          logP(xasg.linear_index()) = -std::numeric_limits<double>::max();
        else logP(xasg.linear_index()) = log(sum);
      }
    }
    

    //! This = other * damping + this * (1-damping) 
    inline void damp(const table_factor& other, double damping) {
      // This factor must be over the same dimensions as the other
      // factor
      if(damping == 0) return;
      ASSERT_EQ(args(), other.args());  
      ASSERT_GT(damping, 0.0);
      ASSERT_LT(damping, 1.0);
      for(size_t i = 0; i < args().size(); ++i) {
        double val = damping * exp(other.logP(i)) + 
          (1-damping) * exp(logP(i));
        ASSERT_GE(val, 0);
        if(val == 0) logP(i) = -std::numeric_limits<double>::max();
        else logP(i) = log(val);
        ASSERT_FALSE( std::isinf(logP(i)) );
        ASSERT_FALSE( std::isnan(logP(i)) );
      }
    }


    //! compute the average l1 norm between to factors
    inline double l1_diff(const table_factor& other) const {
      // This factor must be over the same dimensions as the other
      // factor
      ASSERT_EQ(args(), other.args());  
      double sum = 0;
      for(size_t i = 0; i < args().size(); ++i) {
        sum += fabs(exp(other.logP(i)) - exp(logP(i)));
      }
      return sum / args().size();
    }


    //! compute the l1 norm in log space
    inline double l1_logdiff(const table_factor& other) const {
      ASSERT_EQ(args(), other.args());
      double sum = 0; 
      for(size_t i = 0; i < args().size(); ++i) {
        sum += fabs(other.logP(i) - logP(i));
      }
      return sum / args().size();
    }


    inline assignment_type max_asg() const {
      assignment_type max_asg = args().begin();
      double max_value = logP(max_asg.linear_index());
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        if(logP(asg.linear_index()) > max_value) {
          max_value = logP(asg.linear_index());
          max_asg = asg;
        }
      }
      return max_asg;
    }

    /**
     * Compute the expectation of the table factor
     */
    inline void expectation(std::vector<double>& values) const {
      values.clear();
      values.resize(num_vars(), 0);
      double sum = 0;
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        const double scale = exp(logP(asg.linear_index()));
        sum += scale;
        for(size_t i = 0; i < num_vars(); ++i) {
          values[i] += asg.asg_at(i) * scale;
        }
      }
      // Rescale for normalization
      for(size_t i = 0; i < num_vars(); ++i)  values[i] /= sum;
    } // end of expectation


    /**
     * Draw a sample from the table factor
     */
    inline assignment_type sample() const {
      ASSERT_GT(size(), 0);
      // This factor must be normalized
      const double t = random::rand01();
      ASSERT_GE( t, 0 );
      ASSERT_LT( t, 1 );
      double sum = 0;
      for(size_t i = 0; i < _data.size(); ++i) {
        sum += exp(  logP(i)  );
        if(t <= sum) return assignment_type(args(), i) ;
        ASSERT_LT(sum, 1);
      }
      // Unreachable
      throw("Invalid state reached in sample()");
      return assignment_type();
    } // end of sample
    
    /**
     * Construct a binary agreement factor
     */
    void set_as_agreement(double lambda) {
      ASSERT_EQ(num_vars(), 2);
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        const int diff = abs( int(asg.asg(0)) - int(asg.asg(1)) );
        if( diff > 0) logP(asg.linear_index()) = -lambda;
        else logP(asg.linear_index()) = 0;
      }
    } // end of set_as_agreement
    


    void set_as_laplace(double lambda) {
      ASSERT_EQ(num_vars(), 2);
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        const int diff = abs( int(asg.asg_at(0)) - int(asg.asg_at(1)) );
        logP(asg.linear_index()) = -diff * lambda;
      }
    } // end of set_as_laplace

    void load(iarchive& arc) {
      arc >> _args;
      arc >> _data;      
    }

    void save(oarchive& arc) const {
      arc << _args;
      arc << _data;      
    }
    
  private:
    domain_type _args;
    std::vector<double> _data;
  };  // End of unary factor

  
}; // end of graphlab namespace







template<size_t MAX_DIM>
std::ostream& operator<<(std::ostream& out,
                         const graphlab::table_factor<MAX_DIM>& factor) {
  typedef graphlab::table_factor<MAX_DIM> factor_type;
  typedef typename factor_type::assignment_type assignment_type;
  out << "Table Factor: " << factor.args() << "{" << std::endl;
  for(assignment_type asg = factor.args().begin();
      asg < factor.args().end(); ++asg) {
    out << "\tLogP(" << asg << ")=" << factor.logP(asg) << std::endl;
  }
  return out << "}";
}





#include <graphlab/macros_undef.hpp>
#endif







