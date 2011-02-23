#include <graphlab/util/safe_circular_char_buffer.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/assertions.hpp>

namespace graphlab {

safe_circular_char_buffer::safe_circular_char_buffer(std::streamsize bufsize) 
    :bufsize(bufsize), head(0), tail(0), done(false), iswaiting(false){
  ASSERT_GT(bufsize, 0);
  buffer = (char*)malloc(bufsize);
}

safe_circular_char_buffer::~safe_circular_char_buffer() {
  free(buffer);
}

void safe_circular_char_buffer::stop_reader() {
  mut.lock();
  done = true;
  cond.signal();
  mut.unlock();
}

std::streamsize safe_circular_char_buffer::size() {
  if (tail >= head) return tail - head;
  else if (head < tail) return head + bufsize - tail;
  return 0;
}

std::streamsize safe_circular_char_buffer::write(const char* c, 
                                          std::streamsize clen) {
  if (clen >= bufsize) return 0;
  // is there enough room in the buffer?
  if (bufsize - size() - 1 < clen) return 0;
  
  mut.lock();

  std::streamsize firstcopy = std::min(clen, bufsize - tail);
  memcpy(buffer + tail, c, firstcopy);
  tail += firstcopy;
  if (tail == bufsize) tail = 0;
  if (firstcopy == clen) {
    if (iswaiting) cond.signal();
    mut.unlock();
    return clen;
  }
  
  std::streamsize secondcopy = clen - firstcopy;
  memcpy(buffer, c + firstcopy, secondcopy);
  tail += secondcopy;
  if (iswaiting) cond.signal();
  mut.unlock();
  return clen;
}

std::streamsize safe_circular_char_buffer::introspective_read(char* &s, std::streamsize clen) {
  s = buffer + head;
  // how much we do read?
  // we can go up to the end of the buffer, or until a looparound
  // case 1: no looparound
  // case 2: looparound
  std::streamsize readlen = 0;
  if (tail >= head) {
    readlen = tail - head;
  }
  else {
    readlen = bufsize - head;
  }
  
  // skip the readlen
  head += readlen;
  if (head >= bufsize) head -= bufsize;

  return readlen;
}

std::streamsize safe_circular_char_buffer::blocking_introspective_read(char* &s, std::streamsize clen) {
  // try to read
  size_t ret = introspective_read(s, clen);
  if (ret != 0) return ret;
  
  // if read failed. acquire the lock and try again
  while(1) {
    iswaiting = true;
    mut.lock();
    while (head == tail && !done) cond.wait(mut);
    mut.unlock();    
    iswaiting = false;
    if (done) return 0;
    ret = introspective_read(s, clen);
    if (ret != 0) return ret;
  }
}



}


