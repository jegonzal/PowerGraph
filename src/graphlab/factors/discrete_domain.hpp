#ifndef GRAPHLAB_DISCRETE_DOMAIN
#define GRAPHLAB_DISCRETE_DOMAIN


#include <graphlab/factors/discrete_variable.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  //! Predeclearation of assignment used for discrete domain
  template<size_t MAX_DIM> class discrete_assignment;

  /**
   * This class respresents a discrete discrete_domain over a set of variables.
   */
  template<size_t MAX_DIM>
  class discrete_domain {
  public:

    typedef discrete_assignment<MAX_DIM> assignment_type;    

    //! Make an empy domain
    discrete_domain() : _num_vars(0), _size(0) { }
    //! Make a single variable discrete_domain
    discrete_domain(const discrete_variable& v1) :
      _num_vars(1), _size(0) {
      assert(_num_vars <= MAX_DIM);
      _vars[0] = v1;
      recompute_size();
    }

    //! Make a two variable discrete_domain
    discrete_domain(const discrete_variable& v1, const discrete_variable& v2) :
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

    //! Make a three variable discrete_domain
    discrete_domain(const discrete_variable& v1,
           const discrete_variable& v2,
           const discrete_variable& v3) :
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

    //! Make a discrete_domain from a vector of variables
    discrete_domain(const std::vector<discrete_variable>& variables) :
      _num_vars(variables.size()), _size(0) {
      assert(_num_vars <= MAX_DIM);     
      for(size_t i = 0; i < _num_vars; ++i)       
        _vars[i] = variables[i];
      std::sort(_vars, _vars + std::min(MAX_DIM, _num_vars) );
      recompute_size();
    }

    //! Make a discrete_domain from a set of variables
    discrete_domain(const std::set<discrete_variable>& variables) :
      _num_vars(variables.size()), _size(0) {
      assert(_num_vars <= MAX_DIM); 
      size_t i = 0; 
      foreach(const discrete_variable& var, variables) _vars[i++] = var;
      recompute_size();
    }


    
    discrete_domain& operator+=(const discrete_variable& var) {
      if(_vars[_num_vars - 1] < var) {
        _vars[_num_vars] = var;
        _num_vars++;
        recompute_size();
        return *this;
      }
      return operator+=(discrete_domain(var));
    }

    //! add the other discrete_domain to this discrete_domain
    discrete_domain operator+(const discrete_variable& var) const {
      discrete_domain dom = *this;
      return dom += var;
    }



    //! add the discrete_domain to this discrete_domain
    discrete_domain& operator+=(const discrete_domain& other) {
      if(other.num_vars() == 0) return *this;
      discrete_domain backup = *this;
      _num_vars = 0;
      _size = 0;
      for(size_t i = 0, j = 0; 
          i < backup.num_vars() || j < other.num_vars(); ) {
        assert(_num_vars <= MAX_DIM);
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
          // Unreachable
          assert(false);
        }
      }
      recompute_size();
      return *this;
    }
    
    //! add the other discrete_domain to this discrete_domain
    discrete_domain operator+(const discrete_domain& other) const {
      discrete_domain dom = *this;
      return dom += other;
    }

    
    //! subtract the other discrete_domain from this discrete_domain
    discrete_domain& operator-=(const discrete_domain& other) {
      if(other.num_vars() == 0) return *this;
      
      size_t tmp_num_vars = 0;
      for(size_t i = 0, j = 0; i < _num_vars; ++i ) {
        // advance the other index
        for( ; j < other._num_vars && _vars[i].id() > other._vars[j].id(); ++j);
        if(!(j < other._num_vars && _vars[i].id() == other._vars[j].id())) {
          _vars[tmp_num_vars++] = _vars[i];
        }
      }
      _num_vars = tmp_num_vars;
      recompute_size();
      return *this;
    }

    //! subtract the other discrete_domain from this discrete_domain
    discrete_domain operator-(const discrete_domain& other) const {
      discrete_domain dom = *this;
      return dom -= other;
    }


    discrete_domain intersect(const discrete_domain& other) const {
      discrete_domain new_dom;
      new_dom._num_vars = 0;
      for(size_t i = 0, j = 0;
          i < num_vars() && j < other.num_vars(); ) {
        if(_vars[i] == other._vars[j]) {
          // new discrete_domain gets the variable
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
    const discrete_variable& var(size_t index) const {
      assert(index < _num_vars);
      return _vars[index];
    }
    /** get the index of the variable or returns number of variables
        if the index is not found */
    size_t var_location(size_t var_id) const {
      size_t location = _num_vars;
      for(size_t i = 0; i < _num_vars && !(location < _num_vars); ++i) {
        if(_vars[i].id() == var_id) location = i;
      }
      return location;
    }
    
    //! determine the number of assignments
    size_t size() const { return _size; }


    //! test whether two discrete_domains are equal
    bool operator==(const discrete_domain& other) const {
      bool equal = num_vars() == other.num_vars();
      for(size_t i = 0; equal && i < num_vars(); ++i) {
        equal = equal && var(i) == other.var(i);
        if(!equal) return false;
      }
      if(equal) assert(size() == other.size());
      return equal;
    }
    
    //!  test whether two discrete_domains are not equal
    bool operator!=(const discrete_domain& other) const {
      return !(*this == other);
    }


    //! Get the first assignment in the discrete_domain
    assignment_type begin() const;
    //! Get the second assignment in the discrete_domain
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
    //! Recompute the size of the linear table over this discrete_domainy
    void recompute_size() {
      if(_num_vars > 0) {
        _size = 1;
        for(size_t i = 0; i < _num_vars; ++i) {
          // Require variables to be sorted order
          if(i > 0) assert( _vars[ i-1] < _vars[i]  );
          // and have positive arity
          assert(_vars[i].size() > 0);
          _size *= _vars[i].size();
        }
      } else {
        _size = 0;
      }
    }
  
    size_t _num_vars;
    discrete_variable _vars[MAX_DIM];
    size_t _size;
  };


}; // END OF NAMESPACE


template<size_t MAX_DIM>
std::ostream& operator<<(std::ostream& out,
                         const graphlab::discrete_domain<MAX_DIM>& dom) {
  out << "{";
  for(size_t i = 0; i < dom.num_vars(); ++i) {
    out << dom.var(i);
    if( i < dom.num_vars()-1 ) out << ", ";
  }
  return out << "}";  
}





#include <graphlab/macros_undef.hpp>

#endif
