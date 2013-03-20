#ifndef GRAPHLAB_SERIALIZATION_CONDITIONAL_SERIALIZE_HPP
#define GRAPHLAB_SERIALIZATION_CONDITIONAL_SERIALIZE_HPP
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>
namespace graphlab {

template <typename T>
struct conditional_serialize {
  bool hasval;
  T val;

  conditional_serialize(): hasval(false) { }
  conditional_serialize(T& val): hasval(true), val(val) { }

  conditional_serialize(const conditional_serialize& cs): hasval(cs.hasval), val(cs.val) { }
  conditional_serialize& operator=(const conditional_serialize& cs) {
    hasval = cs.hasval;
    val = cs.val;
    return (*this);
  }
  void save(oarchive& oarc) const {
    oarc << hasval;
    if (hasval) oarc << val;
  }

  void load(iarchive& iarc) {
    iarc >> hasval;
    if (hasval) iarc >> val;
  }
};

};

#endif
