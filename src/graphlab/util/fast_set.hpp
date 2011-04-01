#ifndef GRAPHLAB_FAST_SET_HPP
#define GRAPHLAB_FAST_SET_HPP


#include <iostream>
#include <set>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>



#include <graphlab/macros_def.hpp>
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
    //! Construct an empty set
    fast_set() : nelems(0) { }

    //! Create a stack set with just one element
    fast_set(const T& elem) : nelems(1) { values[0] = elem; }
    
    /**
     * Create a stack from an std set by adding each element one at a
     * time
     */
    template<typename OtherT>
    fast_set(const std::set<OtherT>& other) : nelems(other.size()) { 
      assert(nelems <= MAX_DIM);
      size_t index = 0;
      foreach(const OtherT& elem, other) values[index++] = elem;
    }

    //! Get the begin iterator
    inline T* begin() { return values; }

    //! get the end iterator
    inline T* end() { return values + nelems; }

    //! Get the begin iterator
    const T* begin() const { return values; }

    //! Get the end iterator
    const T* end() const { return values + nelems; }

    //! get the size of the set
    inline size_t size() const { return nelems; }

    //! determine if there are any elements in the set
    bool empty() const { return size() == 0; }

    
    //! test whether the set contains the given element
    bool contains(const T& elem) const {
      for(size_t i = 0; i < nelems; ++i) 
        if(values[i] == elem) return true;
      return false;
    }

    //! test whether the set contains the given set of element
    bool contains(const fast_set& other) const {
      return (*this - other).empty();
    }


    /**
     * Test if this set is contained in other.  If so this returns
     * true. 
     */
    bool operator<=(const fast_set& other) const {
      return (*this - other).empty();
    }


    //! insert an element into this set
    void insert(const T& elem) {
      *this += elem;
    }


    //! insert a range of elements into this set
    void insert(const T* begin, const T* end) {
      size_t other_size = end - begin;  
      assert(other_size <= MAX_DIM);
      fast_set other;
      other.nelems = other_size;
      for(size_t i = 0; i < other_size; ++i) {
        other.values[i] = begin[i];
        if(i > 0) assert(other.values[i] > other.values[i-1]);
      }
      *this += other;
    }


    //! remove an element from the set
    void erase(const T& elem) { *this -= elem; }
    

    //! get the element at a particular location
    const T& operator[](size_t index) const {
      assert(index < nelems);
      return values[index];
    }


    

    //! Take the union of two sets
    inline fast_set operator+(const fast_set& other) const {
      fast_set result;
      size_t i = 0, j = 0;
      while(i < size() && j < other.size()) {
        assert(result.nelems < MAX_DIM);
        if(values[i] < other.values[j])  // This comes first
          result.values[result.nelems++] = values[i++];
        else if (values[i] > other.values[j])  // other comes first
          result.values[result.nelems++] = other.values[j++];
        else {  // both are same
          result.values[result.nelems++] = values[i++]; j++;
        }
      }
      // finish writing this
      while(i < size()) {
        assert(result.nelems < MAX_DIM);
        result.values[result.nelems++] = values[i++];
      }
      // finish writing other
      while(j < other.size()) {
        assert(result.nelems < MAX_DIM);
        result.values[result.nelems++] = other.values[j++];
      }
      return result;
    }


    //! Add the other set to this set
    inline fast_set& operator+=(const fast_set& other) {
      *this = *this + other;
      return *this;
    }


    //! Add an element to this set
    inline fast_set& operator+=(const T& elem) {
      // Find where elem should be inserted
      size_t index = 0;
      for(; index < nelems && values[index] < elem; ++index);      
      assert(index < MAX_DIM);

      // if the element already exists return
      if(index < nelems && values[index] == elem) return *this;

      // otherwise the element does not exist so add it at the current
      // location and increment the number of elements
      T tmp = elem;
      nelems++;
      assert(nelems <= MAX_DIM);
      // Insert the element at index swapping out the rest of the
      // array
      for(size_t i = index; i < nelems; ++i) 
        std::swap(values[i], tmp);      
      // Finished return
      return *this;
    }



    //! Remove the other set from this set
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

    //! subtract the right set form the left set
    fast_set operator-(const fast_set& other) const {
      fast_set result = *this;
      result -= other;
      return result;
    }

    //! Take the intersection of two sets
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

    //! Take the intersection of two sets
    fast_set operator*=(const fast_set& other)  {
      *this = *this * other;
      return *this;
    }

    //! Load the set form file
    void load(iarchive& arc) {
      arc >> nelems;
      assert(nelems <= MAX_DIM);
      for(size_t i = 0; i < nelems; ++i) {
        arc >> values[i];
        if( i > 0 ) assert(values[i] > values[i-1]);
      }
    }
    
    //! Save the set to file
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



#include <graphlab/macros_undef.hpp>
#endif
