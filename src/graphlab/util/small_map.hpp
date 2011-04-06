#ifndef GRAPHLAB_SMALL_MAP_HPP
#define GRAPHLAB_SMALL_MAP_HPP


#include <iostream>
#include <set>
#include <iterator>
#include <algorithm>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/small_set.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  template<typename KeyT, typename ValueT>
  struct  small_map_less {
    bool operator()(const std::pair<KeyT, ValueT>& a, 
                    const std::pair<KeyT, ValueT>& b) const {
      return a.first < b.first;
    } 
  };
  
  template<size_t MAX_DIM, typename KeyT, typename ValueT>
  class small_map : 
    public small_set< MAX_DIM, std::pair<KeyT, ValueT>, 
                      small_map_less<KeyT, ValueT> > {
  public: 
    typedef small_map_less<KeyT, ValueT> less;
    typedef small_set< MAX_DIM, std::pair<KeyT, ValueT>, less> base;
    typedef typename base::value_type value_type;
    typedef typename base::iterator   iterator;
    typedef typename base::const_iterator const_iterator;
    typedef KeyT   key_type;
    typedef ValueT mapped_type;

    using base::begin;
    using base::end;
    using base::operator+;

  public:
    //! construct an empty map
    small_map() { }
    
    //! Construct a map with a single element
    small_map(const key_type& key, const mapped_type& value) :
      base(std::make_pair(key, value)) {  }
    
    //! Lookup an element in the map
    const mapped_type& operator[](const key_type& key) const {
      value_type* const ptr = 
        std::lower_bound(begin(), end(), 
                         std::make_pair(key, mapped_type()),
                         less());
      ASSERT_NE(ptr, end());
      ASSERT_TURE(ptr->first == key);
      return ptr->second;
    }

    //! Lookup an element in the map
    mapped_type& operator[](const key_type& key) {
      value_type* ptr = 
        std::lower_bound(begin(), end(), 
                         std::make_pair(key, mapped_type()),
                         less());
      if(ptr != end() && (key == ptr->first) ) { return ptr->second; }
      // Add the entry to the map
      (*this) += std::make_pair(key, ValueT());
      ptr = 
        std::lower_bound(begin(), end(),
                         std::make_pair(key, mapped_type()),
                         less());
      ASSERT_NE(ptr, end());
      ASSERT_TRUE(ptr->first == key);
      return ptr->second;
    }

    mapped_type& safe_find(const key_type& key) {
      value_type* const ptr = 
        std::lower_bound(begin(), end(), 
                         std::make_pair(key, mapped_type()),
                         less());
      ASSERT_NE(ptr, end());
      ASSERT_TRUE(ptr->first == key);
      return ptr->second;
    }
    


  }; // end of small map


}; // end graphlab

template<size_t MAX_DIM, typename KeyT, typename ValueT>
std::ostream& 
operator<<(std::ostream& out, 
           const graphlab::small_map<MAX_DIM, KeyT, ValueT>& map) {
  typedef std::pair<KeyT, ValueT> pair_type;
  size_t index = 0;
  out << '{';
  foreach(const pair_type& pair, map) {
    out << pair.first << "->" << pair.second;
    if(index+1 < map.size()) out << ", ";
  }
  return out << '}';
}
#include <graphlab/macros_undef.hpp>
#endif
