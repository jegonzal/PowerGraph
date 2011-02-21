#ifndef GRAPHLAB_DENSE_BITSET_HPP
#define GRAPHLAB_DENSE_BITSET_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
  
  /**  \ingroup util_internal
   *  Implements an atomic dense bitset
   */
  class dense_bitset {
  public:
    dense_bitset() : array(NULL), len(0), arrlen(0) {
    }

    dense_bitset(size_t size) : array(NULL), len(size) {
      resize(size);
      clear();
    }

    dense_bitset(const dense_bitset &db) {
      array = NULL;
      len = 0;
      arrlen = 0;
      *this = db;
    }
    
    ~dense_bitset() {free(array);}
  
    inline dense_bitset& operator=(const dense_bitset& db) {
      resize(db.size());
      len = db.len;
      arrlen = db.arrlen;
      memcpy(array, db.array, sizeof(size_t) * arrlen);
      return *this;
    }
  
    inline void resize(size_t n) {
      len = n;
      //need len bits
      arrlen = next_powerof2(n) / sizeof(size_t) + 1;
      array = (size_t*)realloc(array, sizeof(size_t) * arrlen);
    }
  
    inline void clear() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = 0;
    }

    inline void fill() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = -1;
    }

    inline void prefetch(uint32_t b) const{
      __builtin_prefetch(&(array[b / (8 * sizeof(size_t))]));
    }
    inline bool get(uint32_t b) const{
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      return array[arrpos] & (size_t(1) << size_t(bitpos));
    }

    //! Set the bit returning the old value
    inline bool set_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      return __sync_fetch_and_or(array + arrpos, mask) & mask;
    }
    
    //! Set the bit returning the old value
    inline bool set_bit_unsync(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      bool ret = array[arrpos] & mask;
      array[arrpos] |= mask;
      return ret;
    }

    //! Set the state of the bit returning the old value
    inline bool set(uint32_t b, bool value) {
      if (value) return set_bit(b);
      else return clear_bit(b);
    }

    //! Set the state of the bit returning the old value
    inline bool set_unsync(uint32_t b, bool value) {
      if (value) return set_bit_unsync(b);
      else return clear_bit_unsync(b);
    }


    //! Clear the bit returning the old value
    inline bool clear_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t test_mask(size_t(1) << size_t(bitpos)); 
      const size_t clear_mask(~test_mask); 
      return __sync_fetch_and_and(array + arrpos, clear_mask) & test_mask;
    }

    inline bool clear_bit_unsync(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t test_mask(size_t(1) << size_t(bitpos)); 
      const size_t clear_mask(~test_mask); 
      bool ret = array[arrpos] & test_mask;
      array[arrpos] &= clear_mask;
      return ret;
    }

    inline bool first_bit(uint32_t &b) {
      for (size_t i = 0; i < arrlen; ++i) {
        if (array[i]) {
          b = i * (sizeof(size_t) * 8) + first_bit_in_block(array[i]);
          return true;
        }
      }
      return false;
    }

    inline bool next_bit(uint32_t &b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      //try to find the next bit in this block
      bitpos = next_bit_in_block(bitpos, array[arrpos]);
      if (bitpos != 0) {
        b = arrpos * (sizeof(size_t) * 8) + bitpos;
        return true;
      }
      else {
        // we have to loop through the rest of the array
        for (size_t i = arrpos + 1; i < arrlen; ++i) {
          if (array[i]) {
            b = i * (sizeof(size_t) * 8) + first_bit_in_block(array[i]);
            return true;
          }
        }
      }
      return false;
    }

    inline size_t size() const {
      return len;
    }
    
    inline void save(oarchive& oarc) const {
      oarc <<len << arrlen;
      if (arrlen > 0) serialize(oarc, array, arrlen*sizeof(size_t));
    }

    inline void load(iarchive& iarc) {
      if (array != NULL) free(array);
      array = NULL;
      iarc >> len >> arrlen;
      if (arrlen > 0) {
        array = (size_t*)malloc(arrlen*sizeof(size_t));
        deserialize(iarc, array, arrlen*sizeof(size_t));
      }
    }

  private:
   
    inline size_t next_powerof2(size_t val) {
      --val;
      val = val | (val >> 1);
      val = val | (val >> 2);
      val = val | (val >> 4);
      val = val | (val >> 8);
      val = val | (val >> 16);
#ifdef _LP64
      val = val | (val >> 32);
#endif
      return val + 1; 
    }
  
 
    inline static void bit_to_pos(uint32_t b, uint32_t& arrpos, uint32_t& bitpos) {
      // the compiler better optimize this...
      arrpos = b / (8 * sizeof(size_t));
      bitpos = b & (8 * sizeof(size_t) - 1);
    }
  
    // returns 0 on failure
    inline size_t next_bit_in_block(const uint32_t& b, const size_t& block) {
      // use CAS to set the bit
      size_t belowselectedbit = size_t(-1) - (((size_t(1) << b) - 1)|(size_t(1)<<b));
      size_t x = block & belowselectedbit ;
      if (x == 0) return 0;
      else return __builtin_ctzl(x);
    }

    // returns 0 on failure
    inline size_t first_bit_in_block(const size_t& block) {
      // use CAS to set the bit
      if (block == 0) return 0;
      else return __builtin_ctzl(block);
    }

    size_t* array;
    size_t len;
    size_t arrlen;

    template <int len>
    friend class fixed_dense_bitset;
  };

  
  
  
  template <int len>
  class fixed_dense_bitset {
  public:
    fixed_dense_bitset() {
      clear();
    }
    
    fixed_dense_bitset(const fixed_dense_bitset &db) {
      *this = db;
    }
    
    ~fixed_dense_bitset() {}
  
    inline fixed_dense_bitset& operator=(const fixed_dense_bitset& db) {
      memcpy(array, db.array, sizeof(size_t) * arrlen);
      return *this;
    }
  
  
    inline void clear() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = 0;
    }

    inline void fill() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = -1;
    }

    inline void prefetch(uint32_t b) const{
      __builtin_prefetch(&(array[b / (8 * sizeof(size_t))]));
    }
    
    inline bool get(uint32_t b) const{
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      return array[arrpos] & (size_t(1) << size_t(bitpos));
    }

    //! Set the bit returning the old value
    inline bool set_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      return __sync_fetch_and_or(array + arrpos, mask) & mask;
    }
    
    //! Set the bit returning the old value
    inline bool set_bit_unsync(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      bool ret = array[arrpos] & mask;
      array[arrpos] |= mask;
      return ret;
    }

    //! Set the state of the bit returning the old value
    inline bool set(uint32_t b, bool value) {
      if (value) return set_bit(b);
      else return clear_bit(b);
    }

    //! Set the state of the bit returning the old value
    inline bool set_unsync(uint32_t b, bool value) {
      if (value) return set_bit_unsync(b);
      else return clear_bit_unsync(b);
    }


    //! Clear the bit returning the old value
    inline bool clear_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t test_mask(size_t(1) << size_t(bitpos)); 
      const size_t clear_mask(~test_mask); 
      return __sync_fetch_and_and(array + arrpos, clear_mask) & test_mask;
    }

    inline bool clear_bit_unsync(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t test_mask(size_t(1) << size_t(bitpos)); 
      const size_t clear_mask(~test_mask); 
      bool ret = array[arrpos] & test_mask;
      array[arrpos] &= clear_mask;
      return ret;
    }

    inline bool first_bit(uint32_t &b) {
      for (size_t i = 0; i < arrlen; ++i) {
        if (array[i]) {
          b = i * (sizeof(size_t) * 8) + first_bit_in_block(array[i]);
          return true;
        }
      }
      return false;
    }

    inline bool next_bit(uint32_t &b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      //try to find the next bit in this block
      bitpos = next_bit_in_block(bitpos, array[arrpos]);
      if (bitpos != 0) {
        b = arrpos * (sizeof(size_t) * 8) + bitpos;
        return true;
      }
      else {
        // we have to loop through the rest of the array
        for (size_t i = arrpos + 1; i < arrlen; ++i) {
          if (array[i]) {
            b = i * (sizeof(size_t) * 8) + first_bit_in_block(array[i]);
            return true;
          }
        }
      }
      return false;
    }

    inline size_t size() const {
      return len;
    }
    
    inline void save(oarchive& oarc) const {
      oarc <<len << arrlen;
      if (arrlen > 0) serialize(oarc, array, arrlen*sizeof(size_t));
    }

    inline void load(iarchive& iarc) {
      if (array != NULL) free(array);
      array = NULL;
      size_t l;
      size_t arl;
      iarc >> l >> arl;
      ASSERT_EQ(l, len);
      ASSERT_EQ(arl, arrlen);
      if (arrlen > 0) {
        deserialize(iarc, array, arrlen*sizeof(size_t));
      }
    }

  private:
   
    inline size_t next_powerof2(size_t val) {
      --val;
      val = val | (val >> 1);
      val = val | (val >> 2);
      val = val | (val >> 4);
      val = val | (val >> 8);
      val = val | (val >> 16);
#ifdef _LP64
      val = val | (val >> 32);
#endif
      return val + 1; 
    }
  
 
    inline static void bit_to_pos(uint32_t b, uint32_t &arrpos, uint32_t &bitpos) {
      // the compiler better optimize this...
      arrpos = b / (8 * sizeof(size_t));
      bitpos = b & (8 * sizeof(size_t) - 1);
    }
  

    // returns 0 on failure
    inline size_t next_bit_in_block(const uint32_t &b, const size_t &block) {
      // use CAS to set the bit
      size_t belowselectedbit = size_t(-1) - (((size_t(1) << b) - 1)|(size_t(1)<<b));
      size_t x = block & belowselectedbit ;
      if (x == 0) return 0;
      else return __builtin_ctzl(x);
    }

    // returns 0 on failure
    inline size_t first_bit_in_block(const size_t &block) {
      // use CAS to set the bit
      if (block == 0) return 0;
      else return __builtin_ctzl(block);
    }

    static const size_t arrlen;
    size_t array[len / sizeof(size_t) + 1];
  };

  template<int len>
  const size_t fixed_dense_bitset<len>::arrlen = len / sizeof(size_t) + 1;
}
#endif
