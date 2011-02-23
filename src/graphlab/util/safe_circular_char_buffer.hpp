#ifndef SAFE_CIRCULAR_CHAR_BUFFER_HPP
#define SAFE_CIRCULAR_CHAR_BUFFER_HPP
#include <graphlab/rpc/circular_char_buffer.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/assertions.hpp>

namespace graphlab {

/**
A non-resizing circular char buffer
with thread-safe write operations and a single reader 
*/
class safe_circular_char_buffer {
 public:
  safe_circular_char_buffer(std::streamsize bufsize = 1048576 /*1 MB */);

  ~safe_circular_char_buffer();
  
  void stop_reader();

  std::streamsize size();
  
  /**
  Returns 0 if the write doesn't fit
  */
  std::streamsize write(const char* c, std::streamsize clen);


  /**
   * Returns a pointer (through s) and a length of the read.
   * This pointer is a direct pointer into the internal buffer 
   * of this datastructure. The pointer is valid as long as no other operations
   * are performed on this structure.
   * The length of the introspective_read may be less than the number 
   * of bytes requested. Multiple calls to introspective_read may be 
   * necessary to read all data in the buffer. If the function returns 0,
   * the buffer is empty.
   */  
  std::streamsize introspective_read(char* &s, std::streamsize clen);
  
  
  /**
   * Same as introspective read. But blocks until there is something to read
   */
  std::streamsize blocking_introspective_read(char* &s, std::streamsize clen);
  
  
 private:
  char* buffer;
  std::streamsize bufsize; // current size of the buffer

  /** 
   * points to the head of the queue. 
   * Reader reads from here
   */
  std::streamsize head;  
  
  /** 
   * points to one past the end of the queue. 
   * writer writes to here. if tail == head, buffer must be empty
   */
  std::streamsize tail;  

  mutex mut;
  conditional cond;
  
  bool done;
  bool iswaiting;
};

}

#endif

