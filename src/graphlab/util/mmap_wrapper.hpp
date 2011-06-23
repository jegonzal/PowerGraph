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


#ifndef GRAPHLAB_MMAP_WRAPPER_HPP
#define GRAPHLAB_MMAP_WRAPPER_HPP
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <string>
#include <graphlab/logger/assertions.hpp>

namespace graphlab {
  class mmap_wrapper{ 
  public:
    /**
     * mmaps a file into memory if pad > 0, the file will be padded so
     * that it is at least "pad" bytes if diskuncached is true, all
     * writes to the mmapped location will not be cached but will
     * write to physical disk immediately.
     */
    inline mmap_wrapper(std::string file, 
			size_t pad = 0, 
			bool diskuncached = false):
      fname(file), fd(0), ptr(NULL), ptrlen(0) {
      advisetype = MADV_NORMAL;
      fd = 0;
      if (diskuncached) {
#ifdef __APPLE__
        fd = open(file.c_str(), O_RDWR | O_CREAT | O_SYNC, 
            S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH); // chmod 644 (u+rw g+r o+r)
#else
        fd = open(file.c_str(), O_RDWR | O_CREAT | O_DSYNC, 
         		  S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH); // chmod 644 (u+rw g+r o+r)
#endif
      }
      else {
	      fd = open(file.c_str(), O_RDWR | O_CREAT, 
		    S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH); // chmod 644 (u+rw g+r o+r)
      }
      ASSERT_MSG(fd >= 0, strerror(errno));    
    
      extend_file(pad);
    
      // mmap it into memory
      ptrlen = file_length();
      ptr = mmap(0, ptrlen, PROT_READ | PROT_WRITE, MAP_SHARED , fd, 0);
      ASSERT_MSG(ptr != MAP_FAILED, strerror(errno));
    }
  
    void extend_file_and_remap(uint64_t pad) {
      extend_file(file_length() + pad);
      if (remap_nomove() == false) remap();
    }
  
    uint64_t file_length() {
      struct stat statbuf;
      int ret = fstat(fd, &statbuf); ASSERT_EQ(ret, 0);
      return (uint64_t)statbuf.st_size;
    }
  
    inline void* mapped_ptr() {
      return ptr;
    }

    inline bool remap_nomove() {
#ifdef __APPLE__
      // OS X has no mremap? Can we make this better?
      return false;
#else
      size_t newptrlen = file_length();
      void* ret = mremap(ptr, ptrlen, newptrlen,  0);
      if (ret != MAP_FAILED) {
        ptrlen = newptrlen;
        return true;
      }
      else {
        return false;
      }
#endif   
    }
    
    inline void remap() {
#ifdef __APPLE__
      // OS X has no mremap? Can we make this better?
      munmap(ptr, ptrlen);
      size_t newptrlen = file_length();
      ptr = mmap(0, ptrlen, PROT_READ | PROT_WRITE, MAP_SHARED , fd, 0);
      ASSERT_MSG(ptr != MAP_FAILED, strerror(errno));
      ptrlen = newptrlen;
#else
      size_t newptrlen = file_length();
      ptr = mremap(ptr, ptrlen, newptrlen, MREMAP_MAYMOVE);
      ASSERT_MSG(ptr != MAP_FAILED, strerror(errno));
      ptrlen = newptrlen;
#endif
    }
  
    inline void sync(void* start, size_t length) {
      msync(start, length, MS_SYNC);
    }
  
    inline void sync_all() {
      sync(ptr, ptrlen);
    }
  
    inline void background_sync(void* start, size_t length) {
      msync(start, length, MS_ASYNC);
    }
  
    inline void background_sync_all() {
      background_sync(ptr, ptrlen);
    }
  
    inline void close() {
      if (ptr != NULL) {
        sync_all();
        munmap(ptr, ptrlen);
        ::close(fd);
        ptr = NULL;
        fd = 0;
        ptrlen = 0;
      }
    }
  
    inline void prefer_seq_access() {
      madvise(ptr, ptrlen, MADV_SEQUENTIAL);
      advisetype = MADV_SEQUENTIAL;
    }
  
    inline void prefer_random_access() {
      madvise(ptr, ptrlen, MADV_RANDOM);
      advisetype = MADV_RANDOM;
    }
  
    inline void prefer_reset() {
      madvise(ptr, ptrlen, MADV_NORMAL);
      advisetype = MADV_NORMAL;
    }

    inline void prefetch(void* loc, size_t len) {
      madvise(loc, len, MADV_WILLNEED);
    }
  
    inline ~mmap_wrapper() {
      close();
    }
  private:
    
    void extend_file(uint64_t pad) {
      // get the file length
      if (pad > 0 && file_length() < pad) {
        lseek(fd, pad - 1, SEEK_SET);
        const ssize_t error(write(fd, "\0", 1));
        ASSERT_NE(error, ssize_t(-1));        
        ASSERT_GE(file_length(), pad);
      }
    }
    
    std::string fname;

    int fd;
    void* ptr;
    size_t ptrlen;
    int advisetype;
  };

} // end namespace graphlab
#endif

