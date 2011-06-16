#include <graphlab/util/mmap_allocator.hpp>
#include <graphlab/util/mmap_wrapper.hpp>
#include <graphlab/logger/assertions.hpp>
#include <limits>
#include <algorithm>
#define PAGE_SIZE 4096
#define PAGE_SIZE_LOG_2 12
namespace graphlab {


mmap_allocator::mmap_allocator(const std::string fname):closed(false) {
  // create mmap file. minimum length of 1 page.
  mapped_file = new mmap_wrapper(fname, PAGE_SIZE);
  // get the header
  mmap_allocator_impl::mmap_file_header* header = 
            (mmap_allocator_impl::mmap_file_header*) mapped_file->mapped_ptr();
  // create the header if non-existant
  bool newfile = false;
  if (header->utilized_bytes == 0) {
    header->utilized_bytes = 8;
    newfile = true;
  }
  // read the header
  utilized_bytes = header->utilized_bytes;
  file_length = mapped_file->file_length();
  
  // if this is a new file, create the base vector
  if (newfile) {
    ASSERT_EQ(create_vector(8, 128), 0);
  }
}


void mmap_allocator::close() {
  if (closed) return;
  closed = true;
  // update the header
  mmap_allocator_impl::mmap_file_header* header = 
                (mmap_allocator_impl::mmap_file_header*) mapped_file->mapped_ptr();
  header->utilized_bytes = utilized_bytes;
  mapped_file->close();
  delete mapped_file;
}


mmap_allocator_offset_t mmap_allocator::create_vector(uint32_t elemsize, 
                                                      uint64_t start_numel) {
  ASSERT_FALSE(closed);
  mmap_allocator_offset_t ret = mem_alloc(sizeof(mmap_allocator_impl::mmap_vector_header) + elemsize * start_numel);
  void* ptr = get_ptr(ret);
  // create the header
  mmap_allocator_impl::mmap_vector_header* vec_header = 
                                (mmap_allocator_impl::mmap_vector_header*)ptr;

  vec_header->numel = 0;
  vec_header->nextblock = 0;
  vec_header->thisblock_numel = start_numel;
  vec_header->elemsize = elemsize;
  release_ptr(ptr);
  return ret;
}

mmap_allocator_offset_t mmap_allocator::mem_alloc(uint64_t len) {
  ASSERT_FALSE(closed);
  // round len up to nearest multiple of 8 bytes
  mmap_allocator_offset_t ret;
  
  mapped_file_lock.writelock();
  
  len = ((len >> 3) + ((len & 7) != 0)) << 3;
  // do we need to extend the file?
  if (utilized_bytes + len < file_length) {
    // nope!
    ret = utilized_bytes;
    utilized_bytes += len;    
  }
  else {
    // extend file. by how much?
    uint64_t extension = utilized_bytes + len - file_length;
    // round extension up to nearest page size
    extension = ((extension >> PAGE_SIZE_LOG_2) + 
                ((extension & (PAGE_SIZE - 1)) != 0)) << PAGE_SIZE_LOG_2;
    
    mapped_file->extend_file_and_remap(extension);
    file_length += extension;
    ret = utilized_bytes;
    utilized_bytes += len; 
  }
  mapped_file_lock.unlock();
  return ret - sizeof(mmap_allocator_impl::mmap_file_header);
}



mmap_allocator_vector::mmap_allocator_vector(mmap_allocator* allocator, 
                      mmap_allocator_offset_t offset, uint32_t elemsize):allocator(allocator), 
                                                                       offset(offset),
                                                                       released(false){
                                                                         
  // read the vector header
  void* ptr = allocator->get_ptr(offset);
  // create the header
  mmap_allocator_impl::mmap_vector_header* vec_header = 
                      (mmap_allocator_impl::mmap_vector_header*)ptr;
  header = (*vec_header);
  allocator->release_ptr(ptr);
  // make sure the element size matches
  ASSERT_EQ(header.elemsize, elemsize);
  // begin caching the index map
  idx_header_map.resize(1);
  idx_header_map[0].second.thisblock_numel = header.thisblock_numel;
  idx_header_map[0].second.nextblock = header.nextblock;
  
  last_known_index = idx_header_map[0].first + 
                     idx_header_map[0].second.thisblock_numel - 1;
  
}

void mmap_allocator_vector::cache_header_index(idx_t up_to_index_or_end) const {
  ASSERT_FALSE(released);
  lock.writelock();
  size_t last_idx_header_entry = idx_header_map.size() - 1;
   
  // while I have not reached the desired index, and there are still un read blocks
  while(last_known_index <= up_to_index_or_end &&
    idx_header_map[last_idx_header_entry].second.nextblock != 0) {
    // add a new entry to the index
    last_idx_header_entry++;
    idx_header_map.resize(last_idx_header_entry + 1);
  
    idx_header_map[last_idx_header_entry].first = last_known_index + 1;
    // read the intermediate header
    void* ptr = allocator->get_ptr(idx_header_map[last_idx_header_entry - 1].second.nextblock);
    mmap_allocator_impl::mmap_vector_intermediate_header* int_header = 
                  (mmap_allocator_impl::mmap_vector_intermediate_header*)(ptr);
    // copy the header
    idx_header_map[last_idx_header_entry].second = (*int_header);
    allocator->release_ptr(ptr);
    
    // update last index
    last_known_index = idx_header_map[last_idx_header_entry].first + 
                        idx_header_map[last_idx_header_entry].second.thisblock_numel - 1;
  } 
  lock.unlock();
}

mmap_allocator_offset_t mmap_allocator_vector::find_index_pos(idx_t idx) const {
  ASSERT_FALSE(released);
  // if it is a valid index, but I know nothing about it
  if (idx > last_known_index &&
      idx < header.numel) {
    cache_header_index(idx);
    ASSERT_LE(idx, last_known_index);
  }

  if (idx <= last_known_index) {
    lock.readlock();
    // binary search for it in the idx header map
    size_t low = 0; // low may include
    size_t high = idx_header_map.size();  // high must always exclude
    while (high > low) {
      size_t center = (high + low) / 2;
      ASSERT_NE(center, high);  // center must never equals to high. 
      if (idx_header_map[center].first <= idx &&
          idx_header_map[center].first + idx_header_map[center].second.thisblock_numel > idx) {
        // it is in this block!
        // shift by the header
        mmap_allocator_offset_t ret;
        if (center == 0) {
          ret = offset + sizeof(mmap_allocator_impl::mmap_vector_header);
        }
        else {
          ret = idx_header_map[center - 1].second.nextblock + sizeof(mmap_allocator_impl::mmap_vector_intermediate_header);
        }
        
        ret = ret + header.elemsize * (idx - idx_header_map[center].first);
        lock.unlock();
        return ret;
      }
      else if (idx_header_map[center].first >= idx) {
        high = center;
      }
      else if (idx_header_map[center].first + idx_header_map[center].second.thisblock_numel <= idx) {
        low = center + 1;
      }
    }
    lock.unlock();
    ASSERT_MSG(false, "binary search in index map failed! Should be impossible!");
    // unreachable
    return (mmap_allocator_offset_t)(-1);
  }

  if (idx >= header.numel) {
    return (mmap_allocator_offset_t)(-1);
  }
  ASSERT_MSG(false, "Index corruption");
  return (mmap_allocator_offset_t)(-1);;
}

void mmap_allocator_vector::resize(idx_t len) {
  ASSERT_FALSE(released);
  if (len > header.numel) {
    reserve(len);
    lock.writelock();
    header.numel = len;
    lock.unlock();
  }
}

/// Resizes the vector. This only extends.
void mmap_allocator_vector::reserve(idx_t len) {
  ASSERT_FALSE(released);
  // ok. we need to extend. Complete the cache
  cache_header_index(std::numeric_limits<idx_t>::max());
  if (len <= last_known_index + 1) return;
  lock.writelock();  
  // now, allocate a new block
  size_t new_block_len = len - last_known_index - 1;
  mmap_allocator_offset_t ret = allocator->mem_alloc(new_block_len * header.elemsize + sizeof(mmap_allocator_impl::mmap_vector_intermediate_header));
  // create the new header and update the previous header
  idx_header_map[idx_header_map.size() - 1].second.nextblock = ret;
  void* ptr = allocator->get_ptr(idx_header_map[idx_header_map.size() - 1].second.nextblock);
  mmap_allocator_impl::mmap_vector_intermediate_header* new_header = 
                  (mmap_allocator_impl::mmap_vector_intermediate_header*)(ptr);
  new_header->thisblock_numel = new_block_len;
  new_header->nextblock = 0;
  // add this new entry into the index map
  idx_header_map.push_back(std::make_pair(last_known_index + 1, *new_header));
  allocator->release_ptr(ptr);


  
  ASSERT_GE(idx_header_map.size(), 2);
  // update the previous block as well as the index map
  if (idx_header_map.size() == 2) {
  // if size is 2, the previous block is the vector header
    header.nextblock = ret;
    idx_header_map[idx_header_map.size() - 2].second.nextblock = ret;
  }
  else {
    // otherwise size is at least 3. The previous block header is held in the
    // previous previous block's header
    mmap_allocator_offset_t prevblock = idx_header_map[idx_header_map.size() - 3].second.nextblock;
    void* ptr = allocator->get_ptr(prevblock);
    mmap_allocator_impl::mmap_vector_intermediate_header* prev_int_header = 
                  (mmap_allocator_impl::mmap_vector_intermediate_header*)(ptr);
    prev_int_header->nextblock = ret;
    allocator->release_ptr(ptr);
    idx_header_map[idx_header_map.size() - 2].second.nextblock = ret;
  }
  
  // update the last known entry
  last_known_index = idx_header_map[idx_header_map.size() - 1].first + 
                        idx_header_map[idx_header_map.size() - 1].second.thisblock_numel - 1;
  lock.unlock();
}



bool mmap_allocator_vector::set_entry(idx_t idx, const void* val) {
  ASSERT_FALSE(released);
  if (idx >= header.numel) return false;
  mmap_allocator_offset_t ret = mmap_allocator_vector::find_index_pos(idx);

  void* ptr = allocator->get_ptr(ret);
  memcpy(ptr, val, header.elemsize);
  allocator->release_ptr(ptr);
  return true;
}
  
bool mmap_allocator_vector::get_entry(idx_t idx, void* oval) const {
  ASSERT_FALSE(released);
  if (idx >= header.numel) return false;
  mmap_allocator_offset_t ret = mmap_allocator_vector::find_index_pos(idx);
  void* ptr = allocator->get_ptr(ret);
  memcpy(oval, ptr, header.elemsize);
  allocator->release_ptr(ptr);
  return true;
}
  
mmap_allocator_vector::idx_t mmap_allocator_vector::get_all(void* _ptr, idx_t len) const {
  ASSERT_FALSE(released);
  char* ptr = (char*)(_ptr);
  idx_t startlen = len;
  // loop through all the blocks and copy everything
  // do the first block first
  mmap_allocator_offset_t block = offset;
  void* srcptr = allocator->get_ptr(block + sizeof(mmap_allocator_impl::mmap_vector_header));
  size_t numel_to_copy = std::min<uint64_t>(len, header.thisblock_numel);
  len -= numel_to_copy;
  memcpy(ptr, srcptr, numel_to_copy * header.elemsize);
  allocator->release_ptr(srcptr);
  ptr += numel_to_copy * header.elemsize;
  block = header.nextblock;
  
  // walk the linked list
  while(block != 0 && len > 0){
    srcptr = allocator->get_ptr(block);
    mmap_allocator_impl::mmap_vector_intermediate_header* int_header = (mmap_allocator_impl::mmap_vector_intermediate_header*)(srcptr);
    numel_to_copy = std::min<uint64_t>(len, int_header->thisblock_numel);
    memcpy(ptr, 
           (char*)srcptr + sizeof(mmap_allocator_impl::mmap_vector_intermediate_header), 
           numel_to_copy * header.elemsize);
    len -= numel_to_copy;
    ptr += numel_to_copy * header.elemsize;
    block = int_header->nextblock;
    allocator->release_ptr(srcptr);
  }
  return startlen - len;
}

void mmap_allocator_vector::set_all(const void* _ptr, idx_t len) {
  ASSERT_FALSE(released);
  // resize if necessary
  resize(len + 1);
  const char* ptr = (const char*)(_ptr);
  // loop through all the blocks and copy everything
  // do the first block first
  mmap_allocator_offset_t block = offset;
  void* destptr = allocator->get_ptr(block + sizeof(mmap_allocator_impl::mmap_vector_header));
  size_t numel_to_copy = std::min<uint64_t>(len, header.thisblock_numel);
  len -= numel_to_copy;
  memcpy(destptr, ptr, numel_to_copy * header.elemsize);
  allocator->release_ptr(destptr);
  ptr += numel_to_copy * header.elemsize;
  block = header.nextblock;
  
  // walk the linked list
  while(block != 0 && len > 0){
    destptr = allocator->get_ptr(block);
    mmap_allocator_impl::mmap_vector_intermediate_header* int_header = (mmap_allocator_impl::mmap_vector_intermediate_header*)(destptr);
    numel_to_copy = std::min<uint64_t>(len, int_header->thisblock_numel);
    memcpy((char*)destptr + sizeof(mmap_allocator_impl::mmap_vector_intermediate_header),
           ptr, 
           numel_to_copy * header.elemsize);
    len -= numel_to_copy;
    ptr += numel_to_copy * header.elemsize;
    block = int_header->nextblock;
    allocator->release_ptr(destptr);
  }
}
  
  
void mmap_allocator_vector::push_back(const void* ptr) {
  ASSERT_FALSE(released);
  cache_header_index(std::numeric_limits<idx_t>::max());
  // if we are out of room
  if (header.numel == last_known_index + 1) {
    reserve(header.numel * 2);
  }
  header.numel++;
  set_entry(header.numel - 1, ptr);
}


void mmap_allocator_vector::release() {
  if (released) return;
  released = true;
  // write back the header
  void* ptr = allocator->get_ptr(offset);
  // create the header
  mmap_allocator_impl::mmap_vector_header* vec_header = 
                      (mmap_allocator_impl::mmap_vector_header*)ptr;
  (*vec_header) = header;
  allocator->release_ptr(ptr);
}

void mmap_allocator_vector::print_map() {
  std::cout << "vec size: " << header.numel << std::endl;
  std::cout << "element size: " << header.elemsize << std::endl;
  mmap_allocator_offset_t block = offset;
  for (size_t i = 0;i < idx_header_map.size(); ++i) {
    std::cout << block << ": " << idx_header_map[i].second.thisblock_numel << " elements" << std::endl;
    block = idx_header_map[i].second.nextblock;
  }
  
}
}

