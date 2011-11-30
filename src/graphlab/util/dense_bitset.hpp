/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
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
    
    //! Atomically xors a bit with 1
    inline bool xor_bit(uint32_t b) {
      // use CAS to set the bit
      uint32_t arrpos, bitpos;
      bit_to_pos(b, arrpos, bitpos);
      const size_t mask(size_t(1) << size_t(bitpos)); 
      return __sync_fetch_and_xor(array + arrpos, mask) & mask;
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

    struct bit_pos_iterator {
      typedef std::forward_iterator_tag iterator_category;
      typedef uint32_t value_type;
      typedef uint32_t difference_type;
      typedef const uint32_t* pointer;
      typedef const uint32_t& reference;
      uint32_t pos;
      const dense_bitset* db;
      bit_pos_iterator():pos(-1),db(NULL) {}
      bit_pos_iterator(const dense_bitset* const db, uint32_t pos):pos(pos),db(db) {}
      
      uint32_t operator*() const {
        return pos;
      }
      uint32_t operator++(){
        if (db->next_bit(pos) == false) pos = (uint32_t)(-1);
        return pos;
      }
      uint32_t operator++(int){
        uint32_t prevpos = pos;
        if (db->next_bit(pos) == false) pos = (uint32_t)(-1);
        return prevpos;
      }
      bool operator==(const bit_pos_iterator& other) const {
        ASSERT_TRUE(db == other.db);
        return other.pos == pos;
      }
      bool operator!=(const bit_pos_iterator& other) const {
        ASSERT_TRUE(db == other.db);
        return other.pos != pos;
      }
    };
    
    typedef bit_pos_iterator iterator;
    typedef bit_pos_iterator const_iterator;

    
    bit_pos_iterator begin() const {
      uint32_t pos;
      if (first_bit(pos) == false) pos = uint32_t(-1);
      return bit_pos_iterator(this, pos);
    }
    
    bit_pos_iterator end() const {
      return bit_pos_iterator(this, (uint32_t)(-1));
    }

    /** Returns true with b containing the position of the 
        first bit set to true.
        If such a bit does not exist, this function returns false.
    */
    inline bool first_bit(uint32_t &b) const {
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
    inline bool next_bit(uint32_t &b) const {
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


    size_t popcount() const {
      const uint32_t* tmp = reinterpret_cast<const uint32_t*>(array);
      size_t ret = 0;
      for (size_t i = 0;i < arrlen * (sizeof(size_t) / sizeof(uint32_t)); ++i) {
        ret +=  __builtin_popcount(tmp[i]);
      }
      return ret;
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
    inline uint32_t next_bit_in_block(const uint32_t& b, const size_t& block) const {
      // use CAS to set the bit
      size_t belowselectedbit = size_t(-1) - (((size_t(1) << b) - 1)|(size_t(1)<<b));
      size_t x = block & belowselectedbit ;
      if (x == 0) return 0;
      else return (uint32_t)__builtin_ctzl(x);
    }

    // returns 0 on failure
    inline uint32_t first_bit_in_block(const size_t& block) const{
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
    fixed_dense_bitset(const fixed_dense_bitset<len> &db) {
      *this = db;
    }
    
    /// destructor
    ~fixed_dense_bitset() {}
  
    /// Make a copy of the bitset db
    inline fixed_dense_bitset<len>& operator=(const fixed_dense_bitset<len>& db) {
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


    struct bit_pos_iterator {
      typedef std::forward_iterator_tag iterator_category;
      typedef uint32_t value_type;
      typedef uint32_t difference_type;
      typedef const uint32_t* pointer;
      typedef const uint32_t& reference;
      uint32_t pos;
      const fixed_dense_bitset* db;
      bit_pos_iterator():pos(-1),db(NULL) {}
      bit_pos_iterator(const fixed_dense_bitset* const db, uint32_t pos):pos(pos),db(db) {}
      
      uint32_t operator*() const {
        return pos;
      }
      uint32_t operator++(){
        if (db->next_bit(pos) == false) pos = (uint32_t)(-1);
        return pos;
      }
      uint32_t operator++(int){
        uint32_t prevpos = pos;
        if (db->next_bit(pos) == false) pos = (uint32_t)(-1);
        return prevpos;
      }
      bool operator==(const bit_pos_iterator& other) const {
        ASSERT_TRUE(db == other.db);
        return other.pos == pos;
      }
      bool operator!=(const bit_pos_iterator& other) const {
        ASSERT_TRUE(db == other.db);
        return other.pos != pos;
      }
    };
    
    typedef bit_pos_iterator iterator;
    typedef bit_pos_iterator const_iterator;

    
    bit_pos_iterator begin() const {
      uint32_t pos;
      if (first_bit(pos) == false) pos = uint32_t(-1);
      return bit_pos_iterator(this, pos);
    }
    
    bit_pos_iterator end() const {
      return bit_pos_iterator(this, (uint32_t)(-1));
    }

    /** Returns true with b containing the position of the 
        first bit set to true.
        If such a bit does not exist, this function returns false.
    */
    inline bool first_bit(uint32_t &b) const {
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
    inline bool next_bit(uint32_t &b) const {
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
      size_t l;
      size_t arl;
      iarc >> l >> arl;
      ASSERT_EQ(l, len);
      ASSERT_EQ(arl, arrlen);
      if (arrlen > 0) {
        deserialize(iarc, array, arrlen*sizeof(size_t));
      }
    }

    size_t popcount() const {
      const uint32_t* tmp = reinterpret_cast<const uint32_t*>(array);
      size_t ret = 0;
      for (size_t i = 0;i < arrlen * (sizeof(size_t) / sizeof(uint32_t)); ++i) {
        ret +=  __builtin_popcount(tmp[i]);
      }
      return ret;
    }

  private:
   
    inline static size_t next_powerof2(size_t val) {
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
    inline uint32_t next_bit_in_block(const uint32_t &b, const size_t &block) const {
      // use CAS to set the bit
      size_t belowselectedbit = size_t(-1) - (((size_t(1) << b) - 1)|(size_t(1)<<b));
      size_t x = block & belowselectedbit ;
      if (x == 0) return 0;
      else return (uint32_t)__builtin_ctzl(x);
    }

    // returns 0 on failure
    inline uint32_t first_bit_in_block(const size_t &block) const {
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

