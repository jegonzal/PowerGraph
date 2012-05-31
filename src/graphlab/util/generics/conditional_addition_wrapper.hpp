#ifndef GRAPHLAB_CONDITIONAL_ADDITION_WRAPPER_HPP
#define GRAPHLAB_CONDITIONAL_ADDITION_WRAPPER_HPP
#include <algorithm>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>
namespace graphlab {

template <typename T>
struct conditional_addition_wrapper {
 public:
  bool has_value;
  T value;
  conditional_addition_wrapper():has_value(false), value(T()) {};
  explicit conditional_addition_wrapper(const T& t,
                                        bool has_value = true)
          :has_value(has_value), value(t) {};
  
  void set(const T& t) {
    value = t;
    has_value = true;
  }
  void swap(T& t) {
    std::swap(value, t);
    has_value = true;
  }
  void clear() {
    has_value = false;
    value = T();
  }
  conditional_addition_wrapper& operator+=(
                  const conditional_addition_wrapper<T> &c) {
    if (has_value && c.has_value) {
      // if we both have value, do the regular +=
      value += c.value;
    }
    else if (!has_value && c.has_value) {
      // I have no value, but other has value. Use the other
      has_value = true;
      value = c.value;
    }
    return *this;
  }

  conditional_addition_wrapper& operator+=(const T &c) {
    if (has_value) {
      value += c;
    }
    else if (!has_value) {
      // I have no value, but other has value. Use the other
      has_value = true;
      value = c;
    }
    return *this;
  }


  void save(oarchive& oarc) const {
    oarc << has_value;
    if (has_value) oarc << value;
  }


  void load(iarchive& iarc) {
    iarc >> has_value;
    if (has_value) iarc >> value;
    else value = T();
  }
  
};
} // namespace graphlab
#endif