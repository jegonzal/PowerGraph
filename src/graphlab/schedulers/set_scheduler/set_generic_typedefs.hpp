#ifndef SET_GENERIC_TYPEDEFS_HPP
#define SET_GENERIC_TYPEDEFS_HPP
#include <set>
#include <vector>

#include <graphlab/graph/graph.hpp>

#include <graphlab/extern/bitmagic/bm.h>

#include <graphlab/util/stl_util.hpp>

namespace graphlab {
  /* typedef std::set<vertex_id_t> ss_set_type;
     typedef std::set<vertex_id_t>::iterator ss_set_type_iterator ;
  
     inline ss_set_type_iterator begin(const ss_set_type &s) {
     return s.begin();
     }
  
     inline ss_set_type_iterator end(const ss_set_type &s) {
     return s.end();
     }
     inline void ss_insert(ss_set_type &s, vertex_id_t v) {
     s.insert(v);
     }
  
     inline void ss_union(ss_set_type &lhs, const ss_set_type &rhs) {
     lhs = set_union(lhs, rhs);
     }

     inline void ss_intersect(ss_set_type &lhs, const ss_set_type &rhs) {
     lhs = set_intersect(lhs, rhs);
     }
  
     inline void ss_subtract(ss_set_type &lhs, const ss_set_type &rhs) {
     lhs = set_difference(lhs, rhs);
     }
  
     inline size_t ss_size(const ss_set_type &s) {
     return s.size();
     }
  
     inline bool ss_contains(const ss_set_type &s, vertex_id_t v) {
     return s.find(v) != s.end();
     }*/
  


  struct ss_set {
    typedef bm::bvector<> type;
    typedef bm::bvector<>::enumerator iterator;
  };
  typedef ss_set::type ss_set_type;
  typedef ss_set::iterator ss_set_type_iterator;
  
  inline ss_set_type_iterator begin(const ss_set_type &s) {
    return s.first();
  }
  
  inline ss_set_type_iterator end(const ss_set_type &s) {
    return s.end();
  }
  
  inline void ss_insert(ss_set_type &s, vertex_id_t v) {
    s.set_bit(v);
  }

  inline void ss_remove(ss_set_type &s, vertex_id_t v) {
    s.clear_bit(v);
  }
  
  inline void ss_union(ss_set_type &lhs, const ss_set_type &rhs) {
    lhs |= rhs;
  }

  inline void ss_intersect(ss_set_type &lhs, const ss_set_type &rhs) {
    lhs &= rhs;
  }
  
  inline void ss_subtract(ss_set_type &lhs, const ss_set_type &rhs) {
    lhs -= rhs;
  }
  
  inline size_t ss_size(const ss_set_type &s) {
    return s.count();
  }
  
  inline bool ss_contains(const ss_set_type &s, vertex_id_t v) {
    return s.test(v);
  }
  
  
  

  struct ss_small_set {
    typedef std::vector<vertex_id_t> type;
    typedef std::vector<vertex_id_t>::const_iterator iterator;
  };

  typedef ss_small_set::type ss_small_set_type;
  typedef ss_small_set::iterator ss_small_set_type_iterator;

  
  inline ss_small_set_type_iterator begin(const ss_small_set_type &s) {
    return s.begin();
  }
  
  inline ss_small_set_type_iterator end(const ss_small_set_type &s) {
    return s.end();
  }
  
  inline void ss_insert(ss_small_set_type &s, vertex_id_t v) {
    s.push_back(v);
  }
  
  inline void ss_union(ss_small_set_type &lhs, const ss_small_set_type &rhs) {
    std::copy(rhs.begin(), rhs.end(), std::inserter(lhs, lhs.end()));
  }
  
  inline size_t ss_size(const ss_small_set_type &s) {
    return s.size();
  }
  
  inline bool ss_contains(const ss_small_set_type &s, vertex_id_t v) {
    return find(s.begin(), s.end(), v) != s.end();
  }
  
  inline void ss_remove(ss_small_set_type &s, vertex_id_t v) {
    std::vector<vertex_id_t>::iterator i = find(s.begin(), s.end(), v);
    if (i == s.end()) return;
    else {
      *i = s[s.size() - 1];
      s.resize(s.size() - 1);
    }
    
  }
  
};

#endif
