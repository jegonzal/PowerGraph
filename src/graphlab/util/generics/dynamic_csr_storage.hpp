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
#ifndef GRAPHLAB_DYNAMIC_CSR_STORAGE
#define GRAPHLAB_DYNAMIC_CSR_STORAGE

#include <iostream>
#include <vector>
#include <algorithm>

#include <graphlab/util/generics/counting_sort.hpp>
#include <graphlab/util/generics/block_linked_list.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <boost/iterator/permutation_iterator.hpp>

namespace graphlab {
  /**
   * A compact key-value(s) data structure using Compressed Sparse Row format.
   * The key has type size_t and can be assolicated with multiple values of valuetype.
   * The core operation of is querying the list of values associated with the query key *  and returns the begin and end iterators via <code>begin(id)</code>
   * and <code>end(id)</code>.
   */
  template <typename valuetype, typename sizetype=size_t, uint32_t blocksize=64>
  class dynamic_csr_storage {
   public:
     typedef block_linked_list<valuetype, blocksize> block_linked_list;
     typedef typename block_linked_list::iterator iterator;
     typedef typename block_linked_list::const_iterator const_iterator;
     typedef typename block_linked_list::blocktype blocktype;
     typedef valuetype value_type;

   public:
     dynamic_csr_storage() { }

     template<typename idtype>
     dynamic_csr_storage(const std::vector<idtype>& id_vec,
                         const std::vector<valuetype>& value_vec) {
        init(id_vec, value_vec);
     }

     template<typename idtype>
     void init(const std::vector<idtype>& id_vec,
               const std::vector<valuetype>& value_vec) {

      ASSERT_EQ(id_vec.size(), value_vec.size());

      std::vector<sizetype> permute_index;
      // Build index for id -> value 
      // Prefix of the counting array equals to the begin index for each id
      std::vector<sizetype> prefix;
      counting_sort(id_vec, permute_index, &prefix);

      // Fill in the value vector
      typedef boost::permutation_iterator<
               typename std::vector<valuetype>::const_iterator,
               typename std::vector<sizetype>::const_iterator> permute_iterator;

      permute_iterator _begin = boost::make_permutation_iterator(value_vec.begin(), permute_index.begin());
      permute_iterator _end = boost::make_permutation_iterator(value_vec.end(), permute_index.end());

      values.assign(_begin, _end);
      
      sizevec2ptrvec(prefix, value_ptrs);

      // Build the index pointers 
#ifdef DEBUG_CSR
      for (size_t i = 0; i < permute_index.size(); ++i)
        std::cout << permute_index[i] << " ";
      std::cout << std::endl;

      for (size_t i = 0; i < value_ptrs.size(); ++i)
        std::cout << prefix[i] << " ";
      std::cout << std::endl;

      for (permute_iterator it = _begin; it != _end; ++it) {
        std::cout << *it << " ";
      }
      std::cout << std::endl;

      for (size_t i = 0; i < num_keys(); ++i) {
        std::cout << i << ": ";
        iterator it = begin(i);
        while (it != end(i)) {
          std::cout << *it << " "; 
          ++it;
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
#endif
     }

     /**
      * Wrap the index vector and value vector into csr_storage.
      * Check the property of the input vector. 
      * Clean up the input on finish. 
      */
     void wrap(std::vector<sizetype>& valueptr_vec,
               std::vector<valuetype>& value_vec) {

       for (ssize_t i = 1; i < valueptr_vec.size(); ++i) {
         ASSERT_LE(valueptr_vec[i-1], valueptr_vec[i]);
         ASSERT_LT(valueptr_vec[i], value_vec.size());
       }

       values.assign(value_vec.begin(), value_vec.end());
       sizevec2ptrvec(valueptr_vec, value_ptrs);

       std::vector<valuetype>().swap(value_vec);
       std::vector<sizetype>().swap(valueptr_vec);
     }

     /// Number of keys in the storage.
     inline size_t num_keys() const { return value_ptrs.size(); }

     /// Number of values in the storage.
     inline size_t num_values() const { return values.size(); }

     /// Return iterator to the begining value with key == id 
     inline iterator begin(size_t id) {
       return id < num_keys() ? value_ptrs[id] : values.end();
     } 

     /// Return iterator to the ending+1 value with key == id 
     inline iterator end(size_t id) {
       return (id+1) < num_keys() ? value_ptrs[id+1] : values.end();
     }

     /// Return iterator to the begining value with key == id 
     inline const_iterator begin(size_t id) const {
       return id < num_keys() ? value_ptrs[id] : values.end();
     } 

     /// Return iterator to the ending+1 value with key == id 
     inline const_iterator end(size_t id) const {
       return (id+1) < num_keys() ? value_ptrs[id+1] : values.end();
     }


     /// Insert a new value to a given key
     template <typename idtype>
     void insert (const idtype& key, const valuetype& value) {
       // iterator to the insertion position
       iterator ins_iter = values.get_insert_iterator(end(key));

       // add blocks for new key
       while (key >= num_keys()) {
         value_ptrs.push_back(ins_iter);
       }

       // begin_ins_iter and end_ins_iterator point to 
       // defines the range of the new inserted element.
       iterator begin_ins_iter =  values.insert(ins_iter, value);
       iterator end_ins_iter = begin_ins_iter; ++end_ins_iter;

       const uint32_t mid = blocksize/2;

       // Update pointers. 
       // If the begin_ins_iter is moved, the inserted block is splitted.
       // Update pointers to the left of ins_iter:
       if (begin_ins_iter != ins_iter) {
         int scan = key;
         blocktype* oldptr = ins_iter.get_blockptr();
         while (scan >= 0
                && value_ptrs[scan].get_blockptr() == oldptr
                && value_ptrs[scan].get_offset() >= mid) {
           value_ptrs[scan].get_blockptr() = begin_ins_iter.get_blockptr();
           value_ptrs[scan].get_offset() -= mid;
           --scan;
         }
       }
       // Update pointers to the right of ins_iter. 
       // Base case: the pointer of ins_iter is mapped to end_ins_iter. 
       uint32_t oldoffset =  ins_iter.get_offset();
       iterator newiter =  end_ins_iter;
       for (size_t scan = key+1; scan < num_keys(); ++scan) {
         if (value_ptrs[scan] == values.end()) {
           value_ptrs[scan] = end_ins_iter;
         } else if (value_ptrs[scan].get_blockptr() == ins_iter.get_blockptr()) {
           while (oldoffset != value_ptrs[scan].get_offset()) {
             ++oldoffset;
             ++newiter;
           }
           value_ptrs[scan] = newiter;
         } else {
           break;
         }
       }
     }

     /// Insert a range of values to a given key
     template <typename idtype, typename InputIterator>
     void insert (const idtype& key, InputIterator first, InputIterator last) {
       if (last-first == 0) {
         return;
       }

       // iterator to the insertion position
       iterator ins_iter = values.get_insert_iterator(end(key));

       // add blocks for new key
       while (key >= num_keys()) {
         value_ptrs.push_back(ins_iter);
       }

       // begin_ins_iter and end_ins_iterator point to 
       // defines the range of the new inserted element.
       std::pair<iterator,iterator> iter_pair =  values.insert(ins_iter, first, last);
       iterator begin_ins_iter = iter_pair.first;
       iterator end_ins_iter =  iter_pair.second;

       value_ptrs[key] = begin_ins_iter;
       // Update pointers. 
      ASSERT_TRUE(begin_ins_iter == ins_iter);
       // Update pointers to the right of ins_iter. 
       // Base case: the pointer of ins_iter is mapped to end_ins_iter. 
       uint32_t oldoffset =  ins_iter.get_offset();
       iterator newiter =  end_ins_iter;
       for (size_t scan = key+1; scan < num_keys(); ++scan) {
         if (value_ptrs[scan] == values.end()) {
           value_ptrs[scan] = end_ins_iter;
         } else if (value_ptrs[scan].get_blockptr() == ins_iter.get_blockptr()) {
           while (oldoffset != value_ptrs[scan].get_offset()) {
             ++oldoffset;
             ++newiter;
           }
           value_ptrs[scan] = newiter;
         } else {
           break;
         }
       }
     }

     /// Debug print out the content of the storage;
     void print(std::ostream& out) const {
       for (size_t i = 0; i < num_keys(); ++i)  {
         const_iterator iter = begin(i);
          out << i << ": ";
          // out << "begin: " << iter.get_blockptr() << " " << iter.get_offset() << std::endl;
          // out << "end: " << end(i).get_blockptr() << " " << end(i).get_offset() << std::endl;
         while (iter != end(i)) {
           out << *iter <<  " ";
           ++iter;
         }
         out << std::endl;
       }
     }

     ////////////////////// Internal APIs /////////////////
   public:
     /**
      * \internal
      */
     const std::vector<iterator>& get_index() { return value_ptrs; }
     const block_linked_list& get_values() { return values; }

     void swap(dynamic_csr_storage<valuetype, sizetype>& other) {
       value_ptrs.swap(other.value_ptrs);
       values.swap(other.values);
     }

     void clear() {
       std::vector<iterator>().swap(value_ptrs);
       values.clear();
     }

     void load(iarchive& iarc) { }

     void save(oarchive& oarc) const { }

     size_t estimate_sizeof() const {
       return sizeof(value_ptrs) + sizeof(values) + sizeof(sizetype)*value_ptrs.size() + sizeof(valuetype) * values.size();
     }

     void meminfo(std::ostream& out) {
       out << "num values: " <<  (float)num_values()
                 << "\n num blocks: " << values.num_blocks()
                 << "\n block size: " << blocksize
                 << std::endl;
       out << "utilization: " <<  (float)num_values() / (values.num_blocks() * blocksize) << std::endl;
     }

     ///////////////////// Helper Functions /////////////
   private:
     // Convert integer pointers into block_linked_list::value_iterator
     // Assuming all blocks are fully packed.
     void sizevec2ptrvec (const std::vector<sizetype>& ptrs,
                          std::vector<iterator>& out) {
       ASSERT_EQ(out.size(), 0);
       out.reserve(ptrs.size());

       // for efficiency, we advance pointers based on the previous value
       // because block_linked_list is mostly forward_traversal.
       iterator it = values.begin();
       sizetype prev = 0;
       for (size_t i = 0; i < ptrs.size(); ++i) {
         sizetype cur = ptrs[i];
         it += (cur-prev);
         out.push_back(it);
         prev = cur; 
       }
     }

   private:
     std::vector<iterator> value_ptrs;
     block_linked_list values;
  }; // end of class
} // end of graphlab 
#endif
       // Update pointers to the left of ins_iter:
       // Base case: the pionter of ins_iter is mapped to begin_ins_iter
       // if (begin_ins_iter != ins_iter) {
       //   int scan = key-1;
       //   blocktype* oldptr = ins_iter.get_blockptr();
       //   while (scan >= 0 && value_ptrs[scan].get_blockptr() == oldptr) {
       //     // compute the relative distance of old pointers
       //     size_t dist = ins_iter.get_offset() - value_ptrs[scan].get_offset();

       //     // this distance should still hold for new pointers
       //     if (dist <= begin_ins_iter.get_offset()) {
       //       value_ptrs[scan].get_blockptr() = begin_ins_iter.get_blockptr();
       //       value_ptrs[scan].get_offset()  = begin_ins_iter.get_offset()-dist;
       //     } 
       //     else {
       //       // keep the old pointer, update offset
       //       value_ptrs[scan].get_offset() = 
       //           value_ptrs[scan].get_blockptr()->size() - (dist-begin_ins_iter.get_offset());
       //     }
       //     --scan;
       //   }
       // }


