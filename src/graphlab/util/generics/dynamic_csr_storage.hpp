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
       // add blocks for new key
       if (key >= num_keys()) {
         blocktype* newblock = values.append_block();
         while (key >= num_keys()) {
           value_ptrs.push_back(iterator(newblock, 0));
         }
       }

       iterator enditer = values.get_insert_iterator(end(key));
       blocktype* inptr = enditer.get_blockptr();

       bool split = inptr->is_full();
       values.insert(enditer, value);

       // scan right
       int scan = key;
       // Update affected pointers
       if (split) {
         // scan left
         while (scan >= 0
                && value_ptrs[scan].get_blockptr() == inptr
                && value_ptrs[scan].get_offset() >= (blocksize/2)) {
           value_ptrs[scan].get_blockptr() = inptr->next();
           value_ptrs[scan].get_offset() -= (blocksize/2);
           --scan;
         }
       }
       scan = key+1;
       while (scan < num_keys()
              && value_ptrs[scan].get_blockptr() == inptr) {
         ++value_ptrs[scan].get_offset();
         if (split && value_ptrs[scan].get_offset() >= (blocksize/2)) {
           value_ptrs[scan].get_blockptr() = inptr->next();
           value_ptrs[scan].get_offset() -= (blocksize/2);
         }
         ++scan;
       }
       // values.print(std::cerr);
       // print(std::cerr);
     }

     // /// Insert a range of values to a given key
     // template <typename idtype, typename InputIterator>
     // void insert (const idtype& key, InputIterator first, InputIterator last) {
     //   // add blocks for new key
     //   if (key >= num_keys()) {
     //     blocktype* newblock = values.append_block();
     //     while (key >= num_keys()) {
     //       value_ptrs.push_back(iterator(newblock, 0));
     //     }
     //   }

     // }

     /// Debug print out the content of the storage;
     void print(std::ostream& out) {
       for (size_t i = 0; i < num_keys(); ++i)  {
         iterator iter = begin(i);
          out << i << ": ";
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
     std::vector<iterator> get_index() { return value_ptrs; }
     block_linked_list get_values() { return values; }

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
