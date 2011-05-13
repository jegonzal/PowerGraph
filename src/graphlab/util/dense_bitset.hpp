/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_DENSE_BITSET_HPP
#define GRAPHLAB_DENSE_BITSET_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
  
  /**  \ingroup util
   *  Implements an atomic dense bitset
   */
  class dense_bitset {
  public:

    /// Constructs a bitset of 0 length
    dense_bitset() : array(NULL), len(0), arrlen(0) {
    }

    /// Constructs a bitset with 'size' bits. All bits will be cleared.
    dense_bitset(size_t size) : array(NULL), len(size) {
      resize(size);
      clear();
    }

    /// Make a copy of the bitset db
    dense_bitset(const dense_bitset &db) {
      array = NULL;
      len = 0;
      arrlen = 0;
      *this = db;
    }
    
    /// destructor
    ~dense_bitset() {free(array);}
  
    /// Make a copy of the bitset db
    inline dense_bitset& operator=(const dense_bitset& db) {
      resize(db.size());
      len = db.len;
      arrlen = db.arrlen;
      memcpy(array, db.array, sizeof(size_t) * arrlen);
      return *this;
    }
  
    /** Resizes the current bitset to hold n bits.
    Existing bits will not be changed. If the array size is increased,
    the value of the new bits are undefined
    */
    inline void resize(size_t n) {
      len = n;
      //need len bits
      arrlen = next_powerof2(n) / sizeof(size_t) + 1;
      array = (size_t*)realloc(array, sizeof(size_t) * arrlen);
    }
  
    /// Sets all bits to 0
    inline void clear() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = 0;
    }
    
    /// Sets all bits to 1
    inline void fill() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = (size_t)-1;
    }

    /// Prefetches the word containing the bit b
    inline void prefetch(uint32_t b) const{
      __builtin_prefetch(&(array[b / (8 * sizeof(size_t))]));
    }
    
    /// Returns the value of the bit b
    inline bool get(uint32_t b) const{
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      return array[arrpos] & (size_t(1) << size_t(bitpos));
    }

    //! Atomically sets the bit at position b to true returning the old value
    inline bool set_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      return __sync_fetch_and_or(array + arrpos, mask) & mask;
    }
    
    /** Set the bit at position b to true returning the old value.
        Unlike set_bit(), this uses a non-atomic set which is faster,
        but is unsafe if accessed by multiple threads.
    */
    inline bool set_bit_unsync(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      bool ret = array[arrpos] & mask;
      array[arrpos] |= mask;
      return ret;
    }

    //! Atomically sets the state of the bit to the new value returning the old value
    inline bool set(uint32_t b, bool value) {
      if (value) return set_bit(b);
      else return clear_bit(b);
    }

    /** Set the state of the bit returning the old value.
      This version uses a non-atomic set which is faster, but
      is unsafe if accessed by multiple threads.
    */
    inline bool set_unsync(uint32_t b, bool value) {
      if (value) return set_bit_unsync(b);
      else return clear_bit_unsync(b);
    }


    //! Atomically set the bit at b to false returning the old value
    inline bool clear_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t test_mask(size_t(1) << size_t(bitpos)); 
      const size_t clear_mask(~test_mask); 
      return __sync_fetch_and_and(array + arrpos, clear_mask) & test_mask;
    }

    /** Clears the state of the bit returning the old value.
      This version uses a non-atomic set which is faster, but
      is unsafe if accessed by multiple threads.
    */
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

    /** Returns true with b containing the position of the 
        first bit set to true.
        If such a bit does not exist, this function returns false.
    */
    inline bool first_bit(uint32_t &b) {
      for (size_t i = 0; i < arrlen; ++i) {
        if (array[i]) {
          b = (uint32_t)(i * (sizeof(size_t) * 8)) + first_bit_in_block(array[i]);
          return true;
        }
      }
      return false;
    }

    /** Where b is a bit index, this function will return in b,
        the position of the next bit set to true, and return true.
        If all bits after b are false, this function returns false.
    */
    inline bool next_bit(uint32_t &b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      //try to find the next bit in this block
      bitpos = next_bit_in_block(bitpos, array[arrpos]);
      if (bitpos != 0) {
        b = (uint32_t)(arrpos * (sizeof(size_t) * 8)) + bitpos;
        return true;
      }
      else {
        // we have to loop through the rest of the array
        for (size_t i = arrpos + 1; i < arrlen; ++i) {
          if (array[i]) {
            b = (uint32_t)(i * (sizeof(size_t) * 8)) + first_bit_in_block(array[i]);
            return true;
          }
        }
      }
      return false;
    }

    ///  Returns the number of bits in this bitset
    inline size_t size() const {
      return len;
    }
    
    /// Serializes this bitset to an archive
    inline void save(oarchive& oarc) const {
      oarc <<len << arrlen;
      if (arrlen > 0) serialize(oarc, array, arrlen*sizeof(size_t));
    }

    /// Deserializes this bitset from an archive
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
    inline uint32_t next_bit_in_block(const uint32_t& b, const size_t& block) {
      // use CAS to set the bit
      size_t belowselectedbit = size_t(-1) - (((size_t(1) << b) - 1)|(size_t(1)<<b));
      size_t x = block & belowselectedbit ;
      if (x == 0) return 0;
      else return (uint32_t)__builtin_ctzl(x);
    }

    // returns 0 on failure
    inline uint32_t first_bit_in_block(const size_t& block) {
      // use CAS to set the bit
      if (block == 0) return 0;
      else return (uint32_t)__builtin_ctzl(block);
    }

    size_t* array;
    size_t len;
    size_t arrlen;

    template <int len>
    friend class fixed_dense_bitset;
  };

  
  
  /**
  Like bitset, but of a fixed length as defined by the template parameter
  */
  template <int len>
  class fixed_dense_bitset {
  public:
    /// Constructs a bitset of 0 length
    fixed_dense_bitset() {
      clear();
    }
    
   /// Make a copy of the bitset db
    fixed_dense_bitset(const fixed_dense_bitset &db) {
      *this = db;
    }
    
    /// destructor
    ~fixed_dense_bitset() {}
  
    /// Make a copy of the bitset db
    inline fixed_dense_bitset& operator=(const fixed_dense_bitset& db) {
      memcpy(array, db.array, sizeof(size_t) * arrlen);
      return *this;
    }
  
    /// Sets all bits to 0
    inline void clear() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = 0;
    }
    
    /// Sets all bits to 1
    inline void fill() {
      for (size_t i = 0;i < arrlen; ++i) array[i] = -1;
    }

    /// Prefetches the word containing the bit b
    inline void prefetch(uint32_t b) const{
      __builtin_prefetch(&(array[b / (8 * sizeof(size_t))]));
    }
    
    /// Returns the value of the bit b
    inline bool get(uint32_t b) const{
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      return array[arrpos] & (size_t(1) << size_t(bitpos));
    }

    //! Atomically sets the bit at b to true returning the old value
    inline bool set_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      return __sync_fetch_and_or(array + arrpos, mask) & mask;
    }
    
    /** Set the bit at position b to true returning the old value.
        Unlike set_bit(), this uses a non-atomic set which is faster,
        but is unsafe if accessed by multiple threads.
    */
    inline bool set_bit_unsync(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      bool ret = array[arrpos] & mask;
      array[arrpos] |= mask;
      return ret;
    }

    /** Set the state of the bit returning the old value.
      This version uses a non-atomic set which is faster, but
      is unsafe if accessed by multiple threads.
    */
    inline bool set(uint32_t b, bool value) {
      if (value) return set_bit(b);
      else return clear_bit(b);
    }

    /** Set the state of the bit returning the old value.
      This version uses a non-atomic set which is faster, but
      is unsafe if accessed by multiple threads.
    */
    inline bool set_unsync(uint32_t b, bool value) {
      if (value) return set_bit_unsync(b);
      else return clear_bit_unsync(b);
    }


    //! Atomically set the bit at b to false returning the old value
    inline bool clear_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t test_mask(size_t(1) << size_t(bitpos)); 
      const size_t clear_mask(~test_mask); 
      return __sync_fetch_and_and(array + arrpos, clear_mask) & test_mask;
    }

    /** Clears the state of the bit returning the old value.
      This version uses a non-atomic set which is faster, but
      is unsafe if accessed by multiple threads.
    */
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

    /** Returns true with b containing the position of the 
        first bit set to true.
        If such a bit does not exist, this function returns false.
    */
    inline bool first_bit(uint32_t &b) {
      for (size_t i = 0; i < arrlen; ++i) {
        if (array[i]) {
          b = (uint32_t)(i * (sizeof(size_t) * 8)) + first_bit_in_block(array[i]);
          return true;
        }
      }
      return false;
    }

    /** Where b is a bit index, this function will return in b,
        the position of the next bit set to true, and return true.
        If all bits after b are false, this function returns false.
    */
    inline bool next_bit(uint32_t &b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      //try to find the next bit in this block
      bitpos = next_bit_in_block(bitpos, array[arrpos]);
      if (bitpos != 0) {
        b = (uint32_t)(arrpos * (sizeof(size_t) * 8)) + bitpos;
        return true;
      }
      else {
        // we have to loop through the rest of the array
        for (size_t i = arrpos + 1; i < arrlen; ++i) {
          if (array[i]) {
            b = (uint32_t)(i * (sizeof(size_t) * 8)) + first_bit_in_block(array[i]);
            return true;
          }
        }
      }
      return false;
    }
    
    ///  Returns the number of bits in this bitset
    inline size_t size() const {
      return len;
    }
    
    /// Serializes this bitset to an archive
    inline void save(oarchive& oarc) const {
      oarc <<len << arrlen;
      if (arrlen > 0) serialize(oarc, array, arrlen*sizeof(size_t));
    }

    /// Deserializes this bitset from an archive
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
    inline uint32_t next_bit_in_block(const uint32_t &b, const size_t &block) {
      // use CAS to set the bit
      size_t belowselectedbit = size_t(-1) - (((size_t(1) << b) - 1)|(size_t(1)<<b));
      size_t x = block & belowselectedbit ;
      if (x == 0) return 0;
      else return (uint32_t)__builtin_ctzl(x);
    }

    // returns 0 on failure
    inline uint32_t first_bit_in_block(const size_t &block) {
      // use CAS to set the bit
      if (block == 0) return 0;
      else return (uint32_t)__builtin_ctzl(block);
    }

    static const size_t arrlen;
    size_t array[len / sizeof(size_t) + 1];
  };

  template<int len>
  const size_t fixed_dense_bitset<len>::arrlen = len / sizeof(size_t) + 1;
}
#endif

