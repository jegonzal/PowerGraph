#ifndef GRAPHLAB_FAST_SET_HPP
#define GRAPHLAB_FAST_SET_HPP


#include <iostream>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>


namespace graphlab {


  /**
   * This class implements a dense set of fixed maximum size which
   * support quick operations with stack allocation.
   */
  template<size_t MAX_DIM, typename T>
  class fast_set {
  public: // typedefs
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator; 

  public:
    // empty set
    fast_set() : nelems(0) { }
    fast_set(const T& elem) : nelems(1) { values[0] = elem; }

    T* begin() { return values; }
    T* end() { return values + nelems; }

    const T* begin() const { return values; }
    const T* end() const { return values + nelems; }

    size_t size() const { return nelems; }


    bool contains(const T& elem) {
      for(size_t i = 0; i < nelems; ++i) 
        if(values[i] == elem) return true;
      return false;
    }


    const T& operator[](size_t index) const {
      assert(index < nelems);
      return values[index];
    }

    fast_set operator+(const fast_set& other) const {
      fast_set result;
      size_t i = 0, j = 0;
      while(i < size() && j < other.size()) {
        if(values[i] < other.values[j])  // This comes first
          result.values[result.nelems++] = values[i++];
        else if (values[i] > other.values[j])  // other comes first
          result.values[result.nelems++] = other.values[j++];
        else {  // both are same
          result.values[result.nelems++] = values[i++]; j++;
        }
      }
      // finish writing this
      while(i < size()) result.values[result.nelems++] = values[i++];
      // finish writing other
      while(j < other.size()) result.values[result.nelems++] = other.values[j++];
      // We could have segfaulted but check now anyways
      assert(result.nelems <= MAX_DIM);    
      return result;
    }

    fast_set& operator+=(const fast_set& other) {
      *this = *this + other;
      return *this;
    }


    fast_set& operator-=(const fast_set& other) {
      if(other.size() == 0) return *this;    
      // Backup the old nelems and reset nelems
      size_t old_nelems = size(); nelems = 0;
      for(size_t i = 0, j = 0; i < old_nelems; ++i ) {
        // advance the other index
        for( ; j < other.size() && values[i] > other.values[j]; ++j);
        // otherwise check equality or move forward
        if(j >= other.size() || (values[i] != other.values[j])) 
          values[nelems++] = values[i];
      }
      assert(nelems <= MAX_DIM);
      return *this;
    }

    fast_set operator-(const fast_set& other) const {
      fast_set result = *this;
      result -= other;
      return result;
    }


    fast_set operator*(const fast_set& other) const {
      fast_set result; 
      for(size_t i = 0, j = 0;
          i < size() && j < other.size(); ) {
        if(values[i] == other.values[j]) {
          result.values[result.nelems++] = values[i++]; j++;
        } else if(values[i] < other.values[j]) i++; else j++;    
      }
      assert(result.nelems <= MAX_DIM);
      return result;
    }


    fast_set operator*=(const fast_set& other)  {
      *this = *this * other;
      return *this;
    }

    void load(iarchive& arc) {
      arc >> nelems;
      assert(nelems <= MAX_DIM);
      for(size_t i = 0; i < nelems; ++i) {
        arc >> values[i];
        if( i > 0 ) assert(values[i] > values[i-1]);
      }
    }
    
    void save(oarchive& arc) const {
      arc << nelems;
      for(size_t i = 0; i < nelems; ++i) arc << values[i];
    }


 
  private:
    size_t nelems;
    T values[MAX_DIM];

  }; // end of fast_set

}; // end of graphlab namespace

template<size_t MAX_DIM, typename T>
std::ostream&
operator<<(std::ostream& out, const graphlab::fast_set<MAX_DIM, T>& set) {
  out << "{";
  for(size_t i = 0; i < set.size(); ++i) {
    out << set[i];
    if(i + 1 < set.size()) out << ", "; 
  }
  out << "}";
  return out;
}




#endif
