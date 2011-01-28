#include <cstring>
#include <iostream>
#include <graphlab/util/mmap_wrapper.hpp>

using namespace graphlab;


// mac has no strnlen??!!
size_t my_strnlen(const char* s, size_t n)
{
  const char* p = (const char*)memchr(s, 0, n);
  if (p == NULL) return n;
  else return p - s;
}


int main(int argc, char** argv) { 
  mmap_wrapper wrap("testfile.txt",1024);
  char* c = (char*)wrap.mapped_ptr();
  if (my_strnlen(c,1024)  < 128) {
    std::cout << c << "\n";
  }

  strcpy(c, "helloworld\n");
  std::cout << c << "\n";
  getchar();
  wrap.sync_all();
  wrap.close();
}
