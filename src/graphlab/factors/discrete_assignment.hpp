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

#ifndef GRAPHLAB_DISCRETE_ASSIGNMENT
#define GRAPHLAB_DISCRETE_ASSIGNMENT


#include <graphlab/factors/discrete_variable.hpp>
#include <graphlab/factors/discrete_domain.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {

  
  template<size_t MAX_DIM>
  class discrete_assignment {
  public:
    typedef discrete_domain<MAX_DIM> domain_type;
    typedef discrete_variable        variable_type;
    //! Construct an empty discrete_assignment
    discrete_assignment() : _index(0) { }  

    //! Construct a zero discrete_assignment over the domain
    discrete_assignment(const domain_type& args) :
      _args(args), _index(0) {
      for(size_t i = 0; i < args.num_vars(); ++i) 
        _asgs[i] = 0;
    }

    //! Construct a zero discrete_assignment over the domain
    discrete_assignment(const domain_type& args, size_t index) :
      _args(args), _index(index) {
      assert(index < _args.size());
      recompute_asgs();
    }
    
    //! construct an discrete_assignment from two variables
    discrete_assignment(const variable_type& v1, size_t asg1) :
      _args(v1), _index(asg1) {
      assert(asg1 < v1.size());
      _asgs[0] = asg1;
    }

    
    //! construct an discrete_assignment from two variables
    discrete_assignment(const variable_type& v1, size_t asg1, 
                        const variable_type& v2, size_t asg2) :  
      _args(v1, v2), _index(0) {
      set_asg(v1.id(), asg1);
      set_asg(v2.id(), asg2);      
    }

    //! Construct an discrete_assignment from a vector of variables and a
    //! vector of values
    discrete_assignment(const domain_type& args,
                        const std::vector<size_t>& values) :
      _args(args), _index(0) {
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        assert(values[i] < args.var(i).size());
        _asgs[i] = values[i];        
      }
      recompute_linear_index();
    }
    
    
    // //! Construct the union of two discrete_assignments
    // inline discrete_assignment& operator&=(const discrete_assignment& asg2) {
    //   discrete_assignment asg1 = *this;
    //   const domain_type& dom1 = asg1.args();
    //   const domain_type& dom2 = asg2.args();
    //   _args = dom1 + dom2;
    //   _index = 0;
    //   size_t i = 0, j1 = 0, j2 = 0;
    //   for( ; i < _args.num_vars() && 
    //          (j1 < dom1.num_vars() || j2 < dom2.num_vars()); 
    //        ++i) {
    //     // If the the two discrete_assignments share a same variable
    //     if(j1 < dom1.num_vars() && 
    //        _args.var(i) == dom1.var(j1) && 
    //        j2 < dom2.num_vars() &&
    //        _args.var(i) == dom2.var(j2)) {
    //       // Then they must have the same discrete_assignment
    //       //          assert(asg1._asgs[j1] == asg2._asgs[j2]);
    //       _asgs[i] = asg1._asgs[j1];
    //       ++j1; ++j2;
    //     } else if(j1 < dom1.num_vars() &&
    //               _args.var(i) == dom1.var(j1) ) {
    //       _asgs[i] = asg1._asgs[j1];
    //       ++j1;
    //     } else if(j2 < dom2.num_vars() &&
    //               _args.var(i) == dom2.var(j2) ) {
    //       _asgs[i] = asg2._asgs[j2];
    //       ++j2;
    //     } else {
    //       // Unreachable state
    //       assert(false);
    //     }
    //   }
    //   assert(i == _args.num_vars());
    //   assert(j1 == dom1.num_vars());
    //   assert(j2 == dom2.num_vars());
    //   recompute_linear_index();
    //   return *this;
    // }
    // // Construct the union of two discrete_assignments
    // discrete_assignment operator&(const discrete_assignment& other) const {
    //   discrete_assignment new_asg = *this;
    //   return new_asg &= other;
    // }




    //! Construct the union of two discrete_assignments
    inline discrete_assignment operator&(const discrete_assignment& other) const {
      discrete_assignment result(args() + other.args());
      // Require disjoint discrete_assignments
      //      assert(args().size() + other.args().size() == result.size());
      size_t i = 0, j = 0, k = 0;
      while(i < num_vars() && j < other.num_vars()) {        
        // extra increment if necessary
        assert(k < result.num_vars());
        result._asgs[k] = 
          (result.args().var(k) == args().var(i))?
          asg_at(i) : other.asg_at(j);
        // if the variables are the same then the discrete_assignments must
        // also be the same
        assert(!(args().var(i) == other.args().var(j)) ||
               (asg_at(i) == other.asg_at(j)));               
        // move indexs
        i += (args().var(i) == result.args().var(k));
        j += (other.args().var(j) == result.args().var(k));
        k++;
      }
      while(i < num_vars()) 
        result._asgs[k++] = asg_at(i++);
      while(j < other.num_vars()) 
        result._asgs[k++] = other.asg_at(j++); 
      // recompute the linear index of the result
      result.recompute_linear_index();
      return result;
    }
    
    // Construct the union of two discrete_assignments
    inline discrete_assignment& operator&=(const discrete_assignment& other) {
      discrete_assignment tmp = *this & other;
      *this = tmp;
      return *this;
    }
    
    //! Get the variable in the discrete_assignment
    const domain_type& args() const { return _args; }

    //! get the number of variables
    size_t num_vars() const { return _args.num_vars(); }

    //! get the size of the discrete_assignment
    size_t size() const { return _args.size(); }
    
    //! Get the next discrete_assignment
    discrete_assignment& operator++() {
      assert(_index < _args.size());
      // Increment the index
      ++_index;
      // Update the discrete_assignments
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        _asgs[i]= ((_asgs[i] + 1) % _args.var(i).size());
        if(_asgs[i] > 0) { return *this; }
      }
      // Reached end
      make_end();
      return *this;
    }

    //! Uniformly sample a new index value
    void uniform_sample() {
      set_index( random::fast_uniform(size_t(0), size() - 1)  );
    }
    
    //! Get the index of this discrete_assignment
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
      assert(value < _args.var(index).size());
      _asgs[index] = value;
      recompute_linear_index();
    }

    void set_asg_at(size_t index, size_t value) {
      assert(index < _args.num_vars());
      assert(value < _args.var(index).size());
      _asgs[index] = value;
      recompute_linear_index();
    }

    void set_index(size_t index) {
      assert(index < _args.size());
      _index = index;
      recompute_asgs();
    }



    //! Tests whether two discrete_assignments are equal
    bool operator==(const discrete_assignment& other) const {
      return _index == other._index;
    }
    //! Tests whether two discrete_assignments are not equal
    bool operator!=(const discrete_assignment& other) const {
      return _index != other._index;
    }
    //! Tests whether this discrete_assignment is < other
    bool operator<(const discrete_assignment& other) const {
      return _index < other._index;
    }
    //! Make this an ending discrete_assignment
    void make_end() {
      _index = -1;
      // for(size_t i = 0; i < _args.num_vars(); ++i)
      //   _asgs[i] = _args.var(i).size();
    }

    //! Restrict the discrete_assignment to an discrete_assignment over the subdomain
    discrete_assignment restrict(const domain_type& sub_domain) const {
      discrete_assignment other_asg(sub_domain);
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

    //! Update the variables in this discrete_assignment with the values in the
    //! other discrete_assignment
    void update(const discrete_assignment& other) {
      for(size_t i = 0, j = 0;
          i < num_vars() && j < other.num_vars(); ) {
        if(_args.var(i) == other._args.var(j)) {
          _asgs[i] = other._asgs[j]; i++; j++;
        }
        while(i < num_vars() &&
              _args.var(i) < other.args().var(j)) i++;
        while(j < other.num_vars() && 
              other.args().var(j) < _args.var(i)) j++;
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

    //! Recompute the index from the discrete_assignment
    void recompute_linear_index() {
      size_t multiple = 1;
      // Clear the index
      _index = 0;
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        _index += multiple * _asgs[i];
	//        assert(_args.var(i).nasgs > 0);
        multiple *= _args.var(i).size();
      }
    }
    
    //! Recompute the discrete_assignments from the index
    void recompute_asgs() {
      assert(_index < _args.size());
      size_t quotient = _index;
      for(size_t i = 0; i < _args.num_vars(); ++i) {
        _asgs[i] = quotient % _args.var(i).size();
        quotient /= _args.var(i).size();
        // assert(_asgs[i] < _args.var(i).size());
      }
    }


    domain_type _args;
    uint16_t _asgs[MAX_DIM];
    uint32_t _index;
  };


  
  template<size_t MAX_DIM>
  discrete_assignment<MAX_DIM> discrete_domain<MAX_DIM>::begin() const {
    return discrete_assignment<MAX_DIM>(*this);
  }

  template<size_t MAX_DIM>
  discrete_assignment<MAX_DIM> discrete_domain<MAX_DIM>::end() const {
    discrete_assignment<MAX_DIM> ret(*this);
    ret.make_end();
    return ret;
  }


}; //end of namespace graphlab




template<size_t MAX_DIM>
std::ostream& operator<<(std::ostream& out,
                         const graphlab::discrete_assignment<MAX_DIM>& asg) {
  out << "{";
  for(size_t i = 0; i < asg.args().num_vars(); ++i) {
    out << "v_" << asg.args().var(i).id()
        << "=" << asg.asg_at(i);
    if(i < asg.args().num_vars() - 1) out << ", ";
  }
  out << "}=" << asg.linear_index();
  return out;
}




#include <graphlab/macros_undef.hpp>
#endif

