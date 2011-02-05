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
    
      // get the file length
      struct stat statbuf;
      int ret = fstat(fd, &statbuf); ASSERT_EQ(ret, 0);
    
      // seek to the padding point and write a byte
      if (pad > 0 && statbuf.st_size < (int)pad) {
        lseek(fd, pad, SEEK_SET);
        write(fd, " ", 1);
        ret = fstat(fd, &statbuf); ASSERT_EQ(ret, 0);
        ASSERT_GE(statbuf.st_size, pad);
      }
      // mmap it into memory
      ptrlen = statbuf.st_size;
      ptr = mmap(0, ptrlen, PROT_READ | PROT_WRITE, MAP_SHARED , fd, 0);
      ASSERT_MSG(ptr != MAP_FAILED, strerror(errno));
    }
  
    inline void* mapped_ptr() {
      return ptr;
    }

    inline void remap() {
      munmap(ptr, ptrlen);
      ptr = mmap(0, ptrlen, PROT_READ | PROT_WRITE, MAP_SHARED , fd, 0);
      ASSERT_MSG(ptr != MAP_FAILED, strerror(errno));
      madvise(ptr, ptrlen, advisetype);
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
    std::string fname;

    int fd;
    void* ptr;
    size_t ptrlen;
    int advisetype;
  };

} // end namespace graphlab
#endif
