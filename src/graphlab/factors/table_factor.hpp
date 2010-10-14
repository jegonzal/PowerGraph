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

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>



// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

namespace graphlab {

  /** represents a discrete variable */
  struct variable {
    uint32_t id;
    uint32_t arity;
    variable() : id(0), arity(0) { }
    variable(uint32_t id, uint32_t arity) : id(id), arity(arity) { }
    bool operator<(const variable& other) const { return id < other.id; }
    bool operator==(const variable& other) const {
      bool equal_id = id == other.id;
      // if the id are equal than the arity better be equal
      assert((equal_id && (arity == other.arity)) || !equal_id);
      return equal_id;
    }
    bool operator!=(const variable& other) const { return id != other.id; }
    void load(iarchive& arc) { arc >> id; arc >> arity; }
    void save(oarchive& arc) const { arc << id; arc << arity; }
  };












  template<size_t MAX_DIM> class assignment;
















  
  /** Represents a domain over a collection of variables */
  template<size_t MAX_DIM>
  class domain {
  public:
    typedef assignment<MAX_DIM> assignment_type;
    
    //! Make an empy domain
    domain() : _num_vars(0), _size(0) { }
    //! Make a single variable domain
    domain(const variable& v1) :
      _num_vars(1), _size(0) {
      assert(_num_vars <= MAX_DIM);
      _vars[0] = v1;
      recompute_size();
    }

    //! Make a two variable domain
    domain(const variable& v1, const variable& v2) :
      _num_vars(2), _size(0) {
      assert(_num_vars <= MAX_DIM);
      assert(v1 != v2);
      if(v1 < v2) {
        _vars[0] = v1;
        _vars[1] = v2;
      } else {
        _vars[0] = v2;
        _vars[1] = v1;
      }
      recompute_size();
    }

    //! Make a three variable domain
    domain(const variable& v1,
           const variable& v2,
           const variable& v3) :
      _num_vars(3), _size(0) {
      assert(_num_vars <= MAX_DIM);
      assert(v1 != v2);
      assert(v2 != v3);
      assert(v1 != v3);
      
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
      } else if( v3 < v1 && v1 < v2) {
        _vars[0] = v3;
        _vars[1] = v1;
        _vars[2] = v2;
      } else { assert(false); }
      recompute_size();
    }

    //! Make a domain from a vector of variables
    domain(const std::vector<variable>& variables) :
      _num_vars(variables.size()), _size(0) {
      assert(_num_vars <= MAX_DIM);     
      for(size_t i = 0; i < _num_vars; ++i)       
        _vars[i] = variables[i];
      std::sort(_vars, _vars + _num_vars);
      recompute_size();
    }

    //! Make a domain from a set of variables
    domain(const std::set<variable>& variables) :
      _num_vars(variables.size()), _size(0) {
      assert(_num_vars <= MAX_DIM); 
      size_t i = 0; 
      foreach(const variable& var, variables) _vars[i++] = var;
      recompute_size();
    }


    //! add the domain to this domain
    domain& operator+=(const domain& other) {
      if(other.num_vars() == 0) return *this;
      domain backup = *this;
      _num_vars = 0;
      _size = 0;
      for(size_t i = 0, j = 0; 
          i < backup.num_vars() || j < other.num_vars(); ) {
        assert(_num_vars <= MAX_DIM);
        // Both 
        if(i < backup.num_vars() && j < other.num_vars()) {
          if(backup.var(i) < other.var(j))  
            _vars[_num_vars++] = backup.var(i++);
          else if(other.var(j) < backup.var(i))  
            _vars[_num_vars++] = other.var(j++);
          else { _vars[_num_vars++] = backup.var(i++); j++; }
        } else if(i < backup.num_vars()) {
          _vars[_num_vars++] = backup.var(i++);
        } else if(j < other.num_vars()) {
          _vars[_num_vars++] = other.var(j++);
        } else {
          // Unreachable
          assert(false);
        }
      }
      recompute_size();
      return *this;
    }
    
    //! add the other domain to this domain
    domain operator+(const domain& other) const {
      domain dom = *this;
      return dom += other;
    }

    
    //! subtract the other domain from this domain
    domain& operator-=(const domain& other) {
      if(other.num_vars() == 0) return *this;
      
      size_t tmp_num_vars = 0;
      for(size_t i = 0, j = 0; i < _num_vars; ++i ) {
        // advance the other index
        for( ; j < other._num_vars && _vars[i].id > other._vars[j].id; ++j);
        if(!(j < other._num_vars && _vars[i].id == other._vars[j].id)) {
          _vars[tmp_num_vars++] = _vars[i];
        }
      }
      _num_vars = tmp_num_vars;
      recompute_size();
      return *this;
    }

    //! subtract the other domain from this domain
    domain operator-(const domain& other) const {
      domain dom = *this;
      return dom -= other;
    }


    domain intersect(const domain& other) const {
      domain new_dom;
      new_dom._num_vars = 0;
      for(size_t i = 0, j = 0;
          i < num_vars() && j < other.num_vars(); ) {
        if(_vars[i] == other._vars[j]) {
          // new domain gets the variable
          new_dom._vars[new_dom._num_vars] = _vars[i];
          // Everyone advances
          new_dom._num_vars++;  i++; j++;
        } else {
          // otherwise increment one of the variables          
          if(_vars[i] < other._vars[j]) i++; else j++;
        }
      }
      new_dom.recompute_size();      
      return new_dom;
    }
    

    //! Get the number of variables
    size_t num_vars() const { return _num_vars; }

    //! Get the ith variable
    const variable& var(size_t index) const {
      assert(index < _num_vars);
      return _vars[index];
    }
    /** get the index of the variable or returns number of variables
        if the index is not found */
    size_t var_location(size_t var_id) const {
      size_t location = _num_vars;
      for(size_t i = 0; i < _num_vars && !(location < _num_vars); ++i) {
        if(_vars[i].id == var_id) location = i;
      }
      return location;
    }
    
    //! determine the number of assignments
    size_t size() const { return _size; }


    //! test whether two domains are equal
    bool operator==(const domain& other) const {
      bool equal = num_vars() == other.num_vars();
      for(size_t i = 0; equal && i < num_vars(); ++i) {
        equal = equal && var(i) == other.var(i);
        if(!equal) return false;
      }
      if(equal) assert(size() == other.size());
      return equal;
    }
    
    //!  test whether two domains are not equal
    bool operator!=(const domain& other) const {
      return !(*this == other);
    }


    //! Get the first assignment in the domain
    assignment_type begin() const;
    //! Get the second assignment in the domain
    assignment_type end() const;

    void load(iarchive& arc) {
      arc >> _num_vars;
      assert(_num_vars <= MAX_DIM);
      for(size_t i = 0; i < _num_vars; ++i) arc >> _vars[i];
      recompute_size();
    }
    
    void save(oarchive& arc) const {
      arc << _num_vars;
      for(size_t i = 0; i < _num_vars; ++i) arc << _vars[i];
    }

  private:
    //! Recompute the size of the linear table over this domain
    void recompute_size() {
      if(_num_vars > 0) {
        _size = 1;
        for(size_t i = 0; i < _num_vars; ++i) {
          // Require variables to be sorted order
          if(i > 0) assert( _vars[ i-1] < _vars[i]  );
          // and have positive arity
          assert(_vars[i].arity > 0);
          _size *= _vars[i].arity;
        }
      } else {
        _size = 0;
      }
    }
  
    size_t _num_vars;
    variable _vars[MAX_DIM];
    size_t _size;
  };














  

  
  template<size_t MAX_DIM>
  class assignment {
  public:
    typedef domain<MAX_DIM> domain_type;
    
    //! Construct an empty assignment
    assignment() : _index(0) { }  

    //! Construct a zero assignment over the domain
    assignment(const domain_type& args) :
      _args(args), _index(0) {
      for(size_t i = 0; i < args.num_vars(); ++i) 
        _asgs[i] = 0;
    }

    //! Construct a zero assignment over the domain
    assignment(const domain_type& args, size_t index) :
      _args(args), _index(index) {
      assert(index < _args.size());
      recompute_asgs();
    }
    
    //! construct an assignment from two variables
    assignment(const variable& v1, size_t asg1) :
      _args(v1), _index(0) {
      set_asg(v1.id, asg1);
    }

    
    //! construct an assignment from two variables
    assignment(const variable& v1, size_t asg1, 
               const variable& v2, size_t asg2) :  
      _args(v1, v2), _index(0) {
      set_asg(v1.id, asg1);
      set_asg(v2.id, asg2);      
    }

    //! Construct an assignment from a vector of variables and a
    //! vector of values
    assignment(const domain_type& args,
               const std::vector<size_t>& values) :
      _args(args), _index(0) {
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        assert(values[i] < args.var(i).arity);
        _asgs[i] = values[i];        
      }
      recompute_linear_index();
    }
    
    
    //! Construct the union of two assignments
    inline assignment& operator&=(const assignment& asg2) {
      assignment asg1 = *this;
      const domain_type& dom1 = asg1.args();
      const domain_type& dom2 = asg2.args();
      _args = dom1 + dom2;
      _index = 0;
      size_t i = 0, j1 = 0, j2 = 0;
      for( ; i < _args.num_vars() && 
             (j1 < dom1.num_vars() || j2 < dom2.num_vars()); 
           ++i) {
        // If the the two assignments share a same variable
        if(j1 < dom1.num_vars() && 
           _args.var(i) == dom1.var(j1) && 
           j2 < dom2.num_vars() &&
           _args.var(i) == dom2.var(j2)) {
          // Then they must have the same assignment
          assert(asg1._asgs[j1] == asg2._asgs[j2]);
          _asgs[i] = asg1._asgs[j1];
          ++j1; ++j2;
        } else if(j1 < dom1.num_vars() &&
                  _args.var(i) == dom1.var(j1) ) {
          _asgs[i] = asg1._asgs[j1];
          ++j1;
        } else if(j2 < dom2.num_vars() &&
                  _args.var(i) == dom2.var(j2) ) {
          _asgs[i] = asg2._asgs[j2];
          ++j2;
        } else {
          // Unreachable state
          assert(false);
        }
      }
      assert(i == _args.num_vars());
      assert(j1 == dom1.num_vars());
      assert(j2 == dom2.num_vars());
      recompute_linear_index();
      return *this;
    }

    
    // Construct the union of two assignments
    assignment operator&(const assignment& other) const {
      assignment new_asg = *this;
      return new_asg &= other;
    }
    
    //! Get the variable in the assignment
    const domain_type& args() const { return _args; }

    //! get the number of variables
    size_t num_vars() const { return _args.num_vars(); }

    //! get the size of the assignment
    size_t size() const { return _args.size(); }
    
    //! Get the next assignment
    assignment& operator++() {
      assert(_index < _args.size());
      // Increment the index
      ++_index;
      // Update the assignments
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        _asgs[i]= ((_asgs[i] + 1) % _args.var(i).arity);
        if(_asgs[i] > 0) { return *this; }
      }
      // Reached end
      make_end();
      return *this;
    }

    //! Uniformly sample a new index value
    void uniform_sample() {
      set_index( std::rand() % size() );
    }
    
    //! Get the index of this assignment
    size_t linear_index() const { return _index; }
    
    size_t asg(size_t var_id) const {
      size_t index = _args.var_location(var_id);
      assert(index < _args.num_vars());
      return _asgs[index];
    }

    size_t asg_at(size_t index) const {
      assert(index < _args.num_vars());
      return _asgs[index];
    }
    
    void set_asg(size_t var_id, size_t value) {
      size_t index = _args.var_location(var_id);
      assert(index < _args.num_vars());
      assert(value < _args.var(index).arity);
      _asgs[index] = value;
      recompute_linear_index();
    }

    void set_asg_at(size_t index, size_t value) {
      assert(index < _args.num_vars());
      assert(value < _args.var(index).arity);
      _asgs[index] = value;
      recompute_linear_index();
    }

    void set_index(size_t index) {
      assert(index < _args.size());
      _index = index;
      recompute_asgs();
    }



    //! Tests whether two assignments are equal
    bool operator==(const assignment& other) const {
      return _index == other._index;
    }
    //! Tests whether two assignments are not equal
    bool operator!=(const assignment& other) const {
      return _index != other._index;
    }
    //! Tests whether this assignment is < other
    bool operator<(const assignment& other) const {
      return _index < other._index;
    }
    //! Make this an ending assignment
    void make_end() {
      _index = -1;
      for(size_t i = 0; i < _args.num_vars(); ++i)
        _asgs[i] = _args.var(i).arity;
    }

    //! Restrict the assignment to an assignment over the subdomain
    assignment restrict(const domain_type& sub_domain) const {
      assignment other_asg(sub_domain);
      size_t index = 0;
      // Map the variables 
      for(size_t i = 0; i < _args.num_vars() && 
            index < sub_domain.num_vars(); ++i) {
        if(sub_domain.var(index) == _args.var(i)) {
          other_asg._asgs[index] = _asgs[i];
          index++;
        }
      }
      assert(index == sub_domain.num_vars());
      // Recompute the index
      other_asg.recompute_linear_index();
      return other_asg;
    } // end of restrict

    //! Update the variables in this assignment with the values in the
    //! other assignment
    void update(const assignment& other) {
      for(size_t i = 0, j = 0;
          i < num_vars() && j < other.num_vars(); ) {
        if(_args.var(i) == other._args.var(j)) {
          _asgs[i] = other._asgs[j]; i++; j++;
        } else {
          if(_args.var(i) < other._args.var(j)) i++; else j++;
        }
      }
      recompute_linear_index();
    }
    
    void load(iarchive& arc) {
      arc >> _args;
      arc >> _index;
      recompute_asgs();
    }
    
    void save(oarchive& arc) const {
      arc << _args;
      arc << _index;
    }

  private:

    //! Recompute the index from the assignment
    void recompute_linear_index() {
      size_t multiple = 1;
      // Clear the index
      _index = 0;
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        _index += multiple * _asgs[i];
        assert(_args.var(i).arity > 0);
        multiple *= _args.var(i).arity;
      }
    }
    
    //! Recompute the assignments from the index
    void recompute_asgs() {
      assert(_index < _args.size());
      size_t quotient = _index;
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        _asgs[i] = quotient % _args.var(i).arity;
        quotient /= _args.var(i).arity;
        assert(_asgs[i] < _args.var(i).arity);

      }
    }


    domain_type _args;
    size_t _asgs[MAX_DIM];
    size_t _index;
  };










  
  template<size_t MAX_DIM>
  assignment<MAX_DIM> domain<MAX_DIM>::begin() const {
    return assignment<MAX_DIM>(*this);
  }

  template<size_t MAX_DIM>
  assignment<MAX_DIM> domain<MAX_DIM>::end() const {
    assignment<MAX_DIM> ret(*this);
    ret.make_end();
    return ret;
  }



  static const double EPSILON = std::numeric_limits<double>::min();
  static const double LOG_EPSILON = std::log(EPSILON);
  static const double MAX_DOUBLE =  std::numeric_limits<double>::max();
  static const double LOG_MAX_DOUBLE = std::log(MAX_DOUBLE);















  
  /**
   * A factor over up to max_dim dimensions
   */
  template<size_t MAX_DIM>
  class table_factor {
  public:
   
    typedef domain<MAX_DIM>     domain_type;
    typedef assignment<MAX_DIM> assignment_type;
    
    table_factor() { }
    
    table_factor(const domain_type& dom) :
      _args(dom), _data(dom.size()) { 
      //   uniform();
    }

    table_factor(const table_factor& other) :
      _args(other._args), _data(other._data) { }

    table_factor& operator=(const table_factor& other) {
      _args = other._args;
      _data = other._data;
      return *this;
    }

    void set_args(const domain_type& args) {
      _args = args;
      _data.resize(args.size());
      //      uniform();
    }
    
    const domain_type& args() const {
      return _args;
    }

    const double& logP(size_t index) const {
      assert(index < _args.size());
      return _data[index];
    }
    
    double& logP(size_t index) {
      assert(index < _args.size());
      return _data[index];
    }

    const double& logP(const assignment_type& asg) const {
      if(asg.args() == args()) {
        // if the assignment index matches
        size_t index = asg.linear_index();
        assert(index < _data.size());
        return _data[index];
      } else {
        // Restrict the assignment to this domain
        assignment_type sub_asg = asg.restrict(_args);
        return _data[sub_asg.linear_index()];
      }
    }
    

    double& logP(const assignment_type& asg) {
      if(asg.args() == args()) {
        // if the assignment index matches
        size_t index = asg.linear_index();
        assert(index < _data.size());
        return _data[index];
      } else {
        // Restrict the assignment to this domain
        assignment_type sub_asg = asg.restrict(_args);
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
                std::log(1.0/_data.size()));
    }
 
    void uniform(double value) {
      std::fill(_data.begin(), _data.end(), value);
    } 


    
    //! ensure that sum_x this(x) = 1 
    void normalize() {
      // Compute the max value
      double max_value = logP(0);
      for(size_t asg = 0; asg < size(); ++asg) 
        max_value = std::max(max_value, logP(asg));
      assert( !std::isinf(max_value) );
      assert( !std::isnan(max_value) );
      // scale and compute normalizing constant
      double Z = 0.0;
      for(size_t asg = 0; asg < size(); ++asg) {
        logP(asg) -= max_value;
        Z += std::exp(logP(asg));
      }
      assert( !std::isinf(Z) );
      assert( !std::isnan(Z) );
      assert( Z > 0.0);
      double logZ = std::log(Z);
      assert( !std::isinf(logZ) );
      assert( !std::isnan(logZ) );
      // Normalize
      for(size_t asg = 0; asg < size(); ++asg) {
        logP(asg) -= logZ;
        if(logP(asg) < LOG_EPSILON) logP(asg) = -MAX_DOUBLE;
      }
    } // End of normalize
    

    //! this(x) += other(x);
    inline table_factor& operator+=(const table_factor& other) {
      for(assignment_type asg = args().begin(); asg < args().end(); ++asg) {
        logP(asg.linear_index()) =
          log( exp(logP(asg.linear_index())) + exp(other.logP(asg)) );
        assert(logP(asg.linear_index()) <= LOG_MAX_DOUBLE);
        //  assert( !std::isinf( logP(asg.linear_index()) ) );
        //  assert( !std::isnan( logP(asg.linear_index()) ) );
      }
      return *this;
    }

    

    //! this(x) *= other(x);
    inline table_factor& operator*=(const table_factor& other) {
      for(assignment_type asg = args().begin(); asg < args().end(); ++asg) {
        logP(asg.linear_index()) += other.logP(asg);
        assert(logP(asg.linear_index()) <= LOG_MAX_DOUBLE);
        assert( !std::isinf( logP(asg.linear_index()) ) );
        assert( !std::isnan( logP(asg.linear_index()) ) );
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
      for(assignment_type asg = args().begin(); asg < args().end(); ++asg) {
        logP(asg.linear_index()) -= other.logP(asg);
        assert( !std::isinf( logP(asg.linear_index()) ) );
        assert( !std::isnan( logP(asg.linear_index()) ) );
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
      assert(args() + other.args() == joint.args());
      // Initialize the factor to zero so we can use it as an
      // accumulator
      uniform(0);
      for(assignment_type asg = joint.args().begin();
          asg < joint.args().end(); ++asg) {
        const double value =
          std::exp(joint.logP(asg.linear_index()) + other.logP(asg));
        assert( !std::isinf(value) );
        assert( !std::isnan(value) );
        logP(asg) += value;
      }

      for(size_t i = 0; i < _data.size(); ++i) {
        double& sum = _data[i];
        assert(sum >= 0.0);
        if(sum == 0) sum = -MAX_DOUBLE;
        else sum = std::log(sum);
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
      //   if(sum == 0) logP(xasg.linear_index()) = -MAX_DOUBLE;
      //   else logP(xasg.linear_index()) = std::log(sum);
      // }
    }
    

    //! this(x) = other(x, y = asg) 
    inline void condition(const table_factor& other,
                          const assignment_type& asg) {
      assert(args() == other.args() - asg.args());
      for(assignment_type xasg = args().begin(); 
          xasg < args().end(); ++xasg) {
        assignment_type joint = xasg & asg;
        assert(joint.args() == other.args());
        logP(xasg.linear_index()) = other.logP(joint.linear_index());        
      }
    }


    //! this(x) = this(x) other(x, y = asg) 
    inline void times_condition(const table_factor& other,
                                const assignment_type& asg) {
      //      assert(args() == other.args() - asg.args());
      if(asg.num_vars() == 0) {
        *this *= other;
      } else {
        for(assignment_type xasg = args().begin(); 
            xasg < args().end(); ++xasg) {
          assignment_type joint = xasg & asg;
          //          assert(joint.args() == other.args());
          logP(xasg.linear_index()) += other.logP(joint);        
        }
      }
    }

   
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
      assert(ydom.num_vars() > 0);
      
      // Loop over x
      for(assignment_type xasg = args().begin(); 
          xasg < args().end(); ++xasg) {
        double sum = 0;
        for(assignment_type yasg = ydom.begin(); 
            yasg < ydom.end(); ++yasg) {
          assignment_type joint_asg = xasg & yasg;
          sum += std::exp(joint.logP(joint_asg.linear_index()));
        }
        assert( !std::isinf(sum) );
        assert( !std::isnan(sum) );
        assert(sum >= 0.0);
        if(sum == 0) logP(xasg.linear_index()) = -MAX_DOUBLE;
        else logP(xasg.linear_index()) = std::log(sum);
      }
    }

    //! This = other * damping + this * (1-damping) 
    inline void damp(const table_factor& other, double damping) {
      // This factor must be over the same dimensions as the other
      // factor
      assert(args() == other.args());  
      if(damping == 0) return;
      assert(damping >= 0.0);
      assert(damping < 1.0);
      for(size_t i = 0; i < args().size(); ++i) {
        double val = damping * std::exp(other.logP(i)) + 
          (1-damping) * std::exp(logP(i));
        assert(val >= 0);
        if(val == 0) logP(i) = -MAX_DOUBLE;
        else logP(i) = std::log(val);
        assert( !std::isinf(logP(i)) );
        assert( !std::isnan(logP(i)) );
      }
    }


    //! compute the average l1 norm between to factors
    inline double residual(const table_factor& other) const {
      // This factor must be over the same dimensions as the other
      // factor
      assert(args() == other.args());  
      double sum = 0;
      for(size_t i = 0; i < args().size(); ++i) {
        sum += std::abs(std::exp(other.logP(i)) - std::exp(logP(i)));
      }
      return sum / args().size();
    }

    inline double log_residual(const table_factor& other) const {
      assert(args() == other.args());
      double sum = 0; 
      for(size_t i = 0; i < args().size(); ++i) {
        sum += std::abs(other.logP(i) - logP(i));
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

    inline void expectation(std::vector<double>& values) const {
      values.clear();
      values.resize(num_vars(), 0);
      double sum = 0;
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        double scale = std::exp(logP(asg.linear_index()));
        sum += scale;
        for(size_t i = 0; i < num_vars(); ++i) {
          values[i] += asg.asg_at(i) * scale;
        }
      }
      // Rescale for normalization
      for(size_t i = 0; i < num_vars(); ++i)  values[i] /= sum;
    }


    inline assignment_type sample() const {
      assert(size() > 0);
      // This factor must be normalized
      double t = random::rand01();
      assert( t >= 0 && t < 1);
      double sum = 0;
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        sum += std::exp(logP(asg.linear_index()));
        if(t <= sum) return asg;
        assert(sum < 1);
      }
      std::cout << "{";
      for(size_t i = 0; i < num_vars(); ++i) {
        std::cout << args().var(i).id << " ";
      }
      std::cout << "}"  << std::endl;
      // Unreachable
      assert(false);
    }
    
    void set_as_agreement(double lambda) {
      assert(num_vars() == 2);
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        int diff = std::abs( int(asg.asg(0)) - int(asg.asg(1)) );
        if( diff > 0) logP(asg.linear_index()) = -lambda;
        else logP(asg.linear_index()) = 0;
      }
    } // end of set_as_agreement
    


    void set_as_laplace(double lambda) {
      assert(num_vars() == 2);
      for(assignment_type asg = args().begin(); 
          asg < args().end(); ++asg) {
        int diff = std::abs( int(asg.asg_at(0)) - int(asg.asg_at(1)) );
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




std::ostream& operator<<(std::ostream& out, const graphlab::variable& var);


template<size_t MAX_DIM>
std::ostream& operator<<(std::ostream& out,
                         const graphlab::domain<MAX_DIM>& dom) {
  out << "{";
  for(size_t i = 0; i < dom.num_vars(); ++i) {
    out << dom.var(i);
    if( i < dom.num_vars()-1 ) out << ", ";
  }
  return out << "}";  
}

template<size_t MAX_DIM>
std::ostream& operator<<(std::ostream& out,
                         const graphlab::assignment<MAX_DIM>& asg) {
  out << "{";
  for(size_t i = 0; i < asg.args().num_vars(); ++i) {
    out << "v_" << asg.args().var(i).id
        << "=" << asg.asg_at(i);
    if(i < asg.args().num_vars() - 1) out << ", ";
  }
  out << "}=" << asg.linear_index();
  return out;
}


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







