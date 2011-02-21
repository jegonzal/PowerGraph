#include <cassert>
#include <iostream>
#include <string>
#include <graphlab/rpc/circular_char_buffer.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;

int main(int argc, char** argv) {
  circular_char_buffer cbuf(16);
  
  // write peek read cycle
  for (size_t i = 0; i < 1024; ++i) {
    std::string s = "hello world";
    cbuf.write(s.c_str(), s.length());
    if (i % 5 == 0) cbuf.squeeze();
    ASSERT_EQ(cbuf.size(), s.length());
    
    std::string s2;
    ASSERT_EQ(cbuf.peek(s2, s.length()), s.length());
    ASSERT_EQ(s2, s);
    ASSERT_EQ(cbuf.size(), s.length());    
    cbuf.reserve(rand() % 100);
    ASSERT_EQ(cbuf.size(), s.length());    

    s2.clear();
    
    ASSERT_EQ(cbuf.read(s2, s.length()), s.length());
    ASSERT_EQ(s2, s);
    ASSERT_EQ(cbuf.size(), 0);
  }
  
  // do it again! but with a stream
  // write peek read cycle
  cbuf.clear();
  cbuf.squeeze();
  boost::iostreams::stream<circular_char_buffer_device> strm(cbuf);
  
  for (size_t i = 0; i < 1024; ++i) {
    std::string s = "hello world";
    strm << s << std::endl;

    std::string s2;
    std::getline(strm, s2);
    // get the \n
    ASSERT_EQ(s2, std::string("hello world"));
  }
  
  // one more time with an introspective write, expand, read cycle
  for (size_t i = 0;i < 1024; ++i) {
    char* tmp;
    size_t s = cbuf.introspective_write(tmp);
    cbuf.advance_write(s);
    s = cbuf.introspective_write(tmp);
    cbuf.advance_write(s < 25?s:25);
    
    cbuf.consistency_check();
    cbuf.reserve(cbuf.reserved_size() + 37);
    cbuf.align();
    cbuf.introspective_read(tmp, 31);
  }
 
}
