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


#ifndef SPARSE_INDEX_HPP
#define SPARSE_INDEX_HPP


#include <stdint.h>
#include <iostream>


namespace graphlab {


// some version of this could be shared with discrete_assignment.
template<size_t MAX_DIM>
class sparse_index {
public:
  typedef uint16_t*        iterator;
  typedef const uint16_t*  const_iterator;

public: 
  sparse_index() : _num_vars(0) { } 
  explicit sparse_index(size_t num_vars) : 
      _num_vars(num_vars) 
  { 
    DCHECK_LE(_num_vars, MAX_DIM);
    for(size_t i=0; i<_num_vars; ++i)
      _asg[i] = 0;
  }
  sparse_index(size_t const* const begin, size_t const* const end) {
    _num_vars = end - begin;
    DCHECK_LE(_num_vars, MAX_DIM);
    for(size_t i=0; i<_num_vars; ++i)
      _asg[i] = begin[i];
  }
  explicit sparse_index(std::vector<size_t> asg) : 
      _num_vars(asg.size()) 
  {
    DCHECK_LE(_num_vars, MAX_DIM);     
    for(size_t i=0; i<_num_vars; ++i)
      _asg[i] = asg[i];
  }

public: 
  // SEE http://stackoverflow.com/questions/4421706/operator-overloading
  //! Tests whether two sparse_index are equal 
  bool operator==(const sparse_index& other) const {
    DCHECK_EQ(num_vars(), other.num_vars());
    
    for(size_t i=0; i<_num_vars; ++i) 
      if(_asg[i] != other._asg[i]) return false; 
    return true;
  }
  //! Tests whether this sparse_index is < other
  bool operator<(const sparse_index& other) const {
    DCHECK_EQ(num_vars(), other.num_vars());
  
    for(size_t i=0; i<_num_vars; ++i) {
      if(_asg[i] > other._asg[i]) return false;
      else if(_asg[i] < other._asg[i]) return true;
    }
    return false;
  }
  //! Tests whether two sparse_indexs are not equal
  bool operator!=(const sparse_index& other) const {
    return !this->operator==(other);
  }
  //! Tests whether this sparse_index is > other
  bool operator>(const sparse_index& other) const {
    return other.operator<(*this);
  }
  //! Tests whether this sparse_index is <= other
  bool operator<=(const sparse_index& other) const {
    return !this->operator>(other);
  }
  //! Tests whether this sparse_index is >= other
  bool operator>=(const sparse_index& other) const {
    return !this->operator<(other);
  }

  iterator begin() { return &_asg[0]; }
  iterator end() { return &_asg[_num_vars]; }
  const_iterator begin() const { return &_asg[0]; }
  const_iterator end() const { return &_asg[_num_vars]; }

// TODO i dont like these methods being public, but there is little i 
// can do about it since sparse_index doesnt know about the domain
  //! Return the assignment at the specified index 
  // NOTE this is not the assignment for a given variable id; 
  // be mindful of variable reordering 
  inline size_t asg_at(const size_t index) const {
    DCHECK_LT(index, _num_vars);
    return _asg[index];
  }

  //! Set the assignment at the specified index 
  // NOTE this is the assignment for a given index, not for a given variable id; 
  // be mindful of variable reordering 
  inline void set_asg_at(const size_t index, const size_t value) {
    DCHECK_LT(index, _num_vars);
    _asg[index] = value;
  }
  
  inline size_t num_vars() const { return _num_vars; }

  void load(graphlab::iarchive& arc) {
    arc >> _num_vars;
    ASSERT_LE(_num_vars, MAX_DIM);
    for(size_t i = 0; i < _num_vars; ++i) arc >> _asg[i];
  }
  
  void save(graphlab::oarchive& arc) const {
    arc << _num_vars;
    for(size_t i = 0; i < _num_vars; ++i) arc << _asg[i];
  }

  friend std::ostream& operator<<(std::ostream& out, 
      const sparse_index& sa) {
    for(size_t i=0; i < sa._num_vars; ++i) {
      out << sa._asg[i]; 
      if(i < sa._num_vars - 1) out << ", ";
    }
    return out;
  }

private:
  size_t   _num_vars;
  uint16_t _asg[MAX_DIM];
};

} // end of namespace graphlab

#endif // SPARSE_INDEX_HPP
