#include <cassert>
#include <iostream>
#include <string>
#include <graphlab/util/safe_circular_char_buffer.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;

int main(int argc, char** argv) {
  safe_circular_char_buffer cbuf(1024);

  
  // one more time with an introspective write, expand, read cycle
  size_t wrotebytes = 0;
  size_t readbytes = 0;
  for (size_t i = 0;i < 1024; ++i) {
    char v[128];
    char* tmp;
    size_t s = cbuf.write(v, 100);
    wrotebytes += 100;
    while (1) {
      size_t r = cbuf.introspective_read(tmp, 30);
      readbytes += r;
      if (r == 0) break;
    }
  }
  ASSERT_EQ(readbytes, wrotebytes);
 
}

