#ifndef BLOOM_FILTER_HPP
#define BLOOM_FILTER_HPP
#include <graphlab/util/dense_bitset.hpp>

template <size_t len, size_t probes>
class fixed_bloom_filter {
 private:
  fixed_dense_bitset<len> bits;
 public:
  inline fixed_bloom_filter() { }
  
  inline void clear() {
    bits.clear();
  }
  
  inline void insert(uint64_t i) {
    for (size_t i = 0;i < probes; ++i) {
      bits.set_bit_unsync(i % len);
      i = i * 0x9e3779b97f4a7c13LL;
    }
  }
  
  inline bool may_contain(size_t i) {
    for (size_t i = 0;i < probes; ++i) {
      if (bits.get_bit_unsync(i % len) == false) return false;
      i = i * 0x9e3779b97f4a7c13LL;
    }
    return true;
  }

};

#endif
