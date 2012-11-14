#ifndef GRAPHLAB_EXTENSION_DATA_HPP
#define GRAPHLAB_EXTENSION_DATA_HPP

#include <graphlab/graph/distributed_graph.hpp>
#include <boost/variant.hpp>
#include <Eigen/Dense>
#include <string>
#include <map>
#include "../collaborative_filtering/eigen_serialization.hpp"
#include "MurmurHash3.h"

// Here I define the basic "var" variant type.
// which is basically a boost::variant around a double, string, vector and matrix

namespace graphlab {
namespace extension {

// the key here is that we are going to lock down the 
// type system for GraphLab to a wrapped boost::variant
typedef boost::variant<double, 
                       std::string,
                       Eigen::VectorXd,
                       Eigen::MatrixXd> var;

extern var& operator+=(var& value, const var& other);
   

/**
 * Gets a typechecked value from the variant var
 */
template <typename T>
inline const T& get(const var& v) {
  const T* val = boost::get<const T>(&v);
  if (val == NULL) {
    logger_once(LOG_ERROR, "Reading invalid type from var");
    static T t;
    return t;
  }
  return *val;
}

/**
 * Gets a typechecked value from the variant var
 */
template <typename T>
inline T& get(var& v) {
  T* val = boost::get<T>(&v);
  if (val == NULL) {
    logger_once(LOG_ERROR, "Reading invalid type from var");
    static T t;
    return t;
  }
  return *val;
}



//////////////////////////////////////////////////////////////
// Here we briefly escape out tp the global namespace       //
// to define serializers and deserializers for the var      //
//////////////////////////////////////////////////////////////
}}

BEGIN_OUT_OF_PLACE_SAVE(oarc, graphlab::extension::var, value) {
  if ( const double* val = boost::get<double>( &value) ) {
    oarc << char(1) << (*val);
  } else if ( const std::string* val = boost::get<std::string>( &value) ) {
    oarc << char(2) << (*val);
  } else if ( const Eigen::VectorXd* val = boost::get<Eigen::VectorXd>( &value) ) {
    oarc << char(3) << (*val);
  } else if ( const Eigen::MatrixXd* val = boost::get<Eigen::MatrixXd>( &value) ) {
    oarc << char(4) << (*val);
  }
} END_OUT_OF_PLACE_SAVE()


BEGIN_OUT_OF_PLACE_LOAD(iarc, graphlab::extension::var, value) {
  char content_type;
  iarc >> content_type;
  if (content_type == 1) {
    double val; iarc >> val; value = val;
  } else if (content_type == 2) {
    std::string val; iarc >> val; value = val;
  } else if (content_type == 3) {
    Eigen::VectorXd val; iarc >> val; value = val;
  } else if (content_type == 4) {
    Eigen::MatrixXd val; iarc >> val; value = val;
  }
} END_OUT_OF_PLACE_LOAD()
namespace graphlab { namespace extension { 

//////////////////////////////////////////////////////////////
// Returning to your regular progamming                     //
//////////////////////////////////////////////////////////////



typedef uint32_t key_id_type;

/** vars is a dynamic struct with a mapping from string->var.
  * Internally, it is stored as key_id_type->var where the key is
  * a hash value of the string.
  * we assume that the murmurhash will never collide for the small
  * namespaces considered.
 */
inline key_id_type get_id_from_name(const char* key) {
  uint32_t ret = 0;
  MurmurHash3_x86_32((void*)key, strlen(key), 12345, (void*)(&ret));
  return ret;
}
// overload for string
inline key_id_type get_id_from_name(const std::string& key) {
  uint32_t ret = 0;
  MurmurHash3_x86_32((void*)key.c_str(), key.length(), 12345, (void*)(&ret));
  return ret;
}
// overload for int
inline key_id_type get_id_from_name(key_id_type key) {
  return key;
}



/**
 * A dynamic struct storing mappings from string->var where var
 * is a variant.
 * fields can be accessed with operator() or ".field()"
 */
struct vars {
  std::vector<std::pair<key_id_type, var*> > table;
  static var empty_var;
  simple_spinlock lock;
  vars() { }
  ~vars() {
    }

  void clear() {
    lock.lock();
    for(const std::pair<key_id_type, var*>& p: table) {
      delete p.second;
    }
    table.clear();
    lock.unlock();
  }

  vars& operator=(const vars& v) {
    clear();
    for (size_t i = 0;i < v.table.size(); ++i) {
      field(v.table[i].first) = *(v.table[i].second);
    }
    return *this;
  }

  void save(oarchive& oarc) const {
    lock.lock();
    oarc << (size_t)table.size();
    for(const std::pair<key_id_type, var*>& p: table) {
      oarc << p.first << (*p.second);
    }
    lock.unlock();
  }
  
  void load(iarchive& iarc) {
    size_t tsize;
    iarc >> tsize;
    for (size_t i = 0;i < tsize; ++i) {
      key_id_type key; iarc >> key;
      iarc >> field(key);
    }
  }

  var& operator()(const std::string& key) {
    return field(key);
  }
  const var& operator()(const std::string& key) const {
    return field(key);
  }
  var& operator()(const char* key) {
    return field(key);
  }
  const var& operator()(const char* key) const {
    return field(key);
  }
  var& operator()(key_id_type key) {
    return field(key);
  }
  const var& operator()(key_id_type key) const {
    return field(key);
  }

  var& field(const std::string& _key) {
    key_id_type key = get_id_from_name(_key);
    return field(key);
  }

  const var& field(const std::string& _key) const {
    key_id_type key = get_id_from_name(_key);
    return field(key);
  }
  
  var& field(const char* _key) {
    key_id_type key = get_id_from_name(_key);
    return field(key);
  }
  const var& field(const char* _key) const {
    key_id_type key = get_id_from_name(_key);
    return field(key);
  }
   
  var& field(key_id_type key) {
    lock.lock();

    for(std::pair<key_id_type , var*>& p: table) {
      if (p.first == key) {
        lock.unlock();
        return *(p.second);
      }
    }
    var* ret = new var;
    // force slow resize to limit memory usage
    // assume that field creation is not a common operation.
    table.reserve(table.size() + 1);
    table.push_back(std::make_pair(key, ret));
    lock.unlock();
    return *ret;
  }

  const var& field(key_id_type key) const {
    for(const std::pair<key_id_type , var*>& p: table) {
      if (p.first == key) {
        return *(p.second);
      }
    }
    return empty_var;
  }

}; 



typedef distributed_graph<vars, vars> internal_graph_type;

} // namespace extension
} // namespace graphlab

#endif
