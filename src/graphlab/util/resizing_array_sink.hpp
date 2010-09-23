#ifndef RESIZING_COUNTING_SINK
#define RESIZING_COUNTING_SINK
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/categories.hpp>

namespace graphlab {
struct resizing_array_sink{
  resizing_array_sink(size_t initial) { 
    str = (char*)(malloc(initial));
    len = 0;
    size = initial;
  }
  char *str;
  size_t len;
  size_t size;
  typedef char        char_type;
  struct category: public boost::iostreams::device_tag,
                   public boost::iostreams::output,
                   public boost::iostreams::multichar_tag { };

 /** the optimal buffer size is 0. */
  inline std::streamsize optimal_buffer_size() const { return 0; }

  inline std::streamsize write(const char* s, std::streamsize n) {
    if (len + n > size) {
      // double in length if we need more buffer
      size = 2 * (len + n);
      str = (char*)realloc(str, size);
    }
    memcpy(str + len, s, n);
    len += n;
    return n;
  }
};
}
#endif

