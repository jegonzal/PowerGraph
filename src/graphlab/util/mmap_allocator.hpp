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


#ifndef GRAPHLAB_MMAP_ALLOCATOR_HPP
#define GRAPHLAB_MMAP_ALLOCATOR_HPP

#include <vector>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/mmap_wrapper.hpp>

namespace graphlab {

typedef uint64_t mmap_allocator_offset_t;
  
namespace mmap_allocator_impl {
  
  struct mmap_file_header {
    uint64_t utilized_bytes;
  };
  
  struct mmap_vector_header {
    mmap_allocator_offset_t nextblock;
    uint64_t thisblock_numel;
    uint32_t elemsize;
    uint32_t numel;    
  };
  
  struct mmap_vector_intermediate_header {
    uint64_t thisblock_numel;
    mmap_allocator_offset_t nextblock;
  };
}


// forward declaration
class mmap_allocator;
  
/**
 * A mmaped vector as created by mmap_allocator.
 * This object cannot be created directly and can only be created by the 
 * mmap_allocator. 
 */
class mmap_allocator_vector {
 public:
  typedef uint32_t idx_t;
  /** Sets the value of an entry. 'element size' bytes will be read and stored
   * from the 'val' pointer. This function returns true on success and 
   * false if the entry idx does not exist.
   */
  bool set_entry(idx_t idx, const void* val);
  
  /** Reads an entry. 'element size' bytes will be written into oval.
   * User must ensure that oval has sufficient room to store the data.
   * This function returns true on success and if the entry idx does not exist.
   */
  bool get_entry(idx_t idx, void* oval) const;
  
  /// Returns the element size
  idx_t get_elem_size() const {
    return header.elemsize;
  }
  
  /// Returns the number of elements in the vector
  idx_t get_len() const {
    return header.numel;
  }
  
  /// Resizes the vector. This only extends.
  void resize(idx_t len);
  
  /// Reserves at least len elements
  void reserve(idx_t len);
  
  /** Performs a complete dump of the vector of up to 'numel' units.
   * ptr must be able to store at least numel * get_elem_size() bytes.
   * Returns the actual number of elements copied
   */
  idx_t get_all(void* ptr, idx_t numel) const;
  
  /** Reverse of get_all(). Dumps the data in the ptr pointer to disk,
   * updating the number of elements. 
   * numel * get_elem_size() bytes will be read from the pointer.
   */
  void set_all(const void* ptr, idx_t numel);
  
  /**
   * Inserts an additional element. get_elem_size() bytes will be read
   * from the pointer
   */
  void push_back(const void* ptr);
  
  /** Releases this vector. All vectors must either be released
   * or go out of scope before the mmap wrapper can be closed
  */
  void release();
  
  void print_map();
  
  ~mmap_allocator_vector() {
    release();
  }
private:
  mmap_allocator_vector() { }
  mmap_allocator_vector(mmap_allocator* allocator, 
                        mmap_allocator_offset_t offset, uint32_t elemsize);
  mutable mmap_allocator* allocator; /// the allocator this is associated with
  mmap_allocator_offset_t offset; /// the offset address of this vector
  
  mmap_allocator_impl::mmap_vector_header header; /// the vector header
  mutable idx_t last_known_index; /// The last known index in the header map
  
  mutable std::vector<std::pair<idx_t, mmap_allocator_impl::mmap_vector_intermediate_header> > idx_header_map;
  
  /// Caches the index map up to an including the entry up_to_index_or_end
  void cache_header_index(idx_t up_to_index_or_end) const;
  
  /// Get the offset of an index. Returns -1 if does not exist
  mmap_allocator_offset_t find_index_pos(idx_t idx) const;
  
  rwlock lock;
  bool released;
  friend class mmap_allocator;
};


/**
 * A file storage manager built around storing vectors of arbitrary 
 * shape/length in an mmaped file.
 * get_vector reads an existing vector in a particular offset in the file.
 * There is always a vector of element size 8 at offset 0 which can be used
 * for indexing.
 * create_vector creates a vector of a particular element size returning the
 * offset.
 * 
 * This object is parallel and permits multiple vectors in the same file to be
 * created and accessed simultaneously. However, there should be at most one
 * vector object created per offset.
 */
class mmap_allocator {
public:
  /// Opens or creates an storage file
  mmap_allocator(const std::string fname);
  
  /// Closes the file.
  void close();
  
  /// Obtains a vector stored at a particular file offset
  inline mmap_allocator_vector get_vector(mmap_allocator_offset_t offset, 
                                          uint32_t elemsize) const {
    return mmap_allocator_vector(const_cast<mmap_allocator*>(this), offset, elemsize);
  }
  
  /// Creates a vector with a particular element size. Returns the offset
  mmap_allocator_offset_t create_vector(uint32_t elemsize, uint64_t start_numel);
  
  inline ~mmap_allocator() {
    close();
  }
private:
  mutable mmap_wrapper* mapped_file;
  
  uint64_t utilized_bytes;
  uint64_t file_length;
  // to be acquired when the root pointer is accessed
  mutable rwlock mapped_file_lock;
  
  bool closed;
  
  /// Requests an allocation of 'len' bytes. Returns the offset.
  mmap_allocator_offset_t mem_alloc(uint64_t len);
  
  /** Get a ptr at this offset. Pointer must be released with release_ptr
   * when done. TODO: Current implementation only permits each thread
   * to have only one active pointer. 
   */
  inline void* get_ptr(uint64_t offset) const {
    ASSERT_FALSE(closed);
    mapped_file_lock.readlock();
    return (char*)(mapped_file->mapped_ptr()) + 
            sizeof(mmap_allocator_impl::mmap_file_header) + offset;
  }
  
  inline void release_ptr(void* ptr) const {
    ASSERT_FALSE(closed);
    mapped_file_lock.unlock();
  }
  
  friend class mmap_allocator_vector;
  
  
};

}
#endif