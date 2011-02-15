#ifndef CHAR_COUNTING_SINK
#define CHAR_COUNTING_SINK
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/categories.hpp>

namespace graphlab {
  
/**
 \ingroup util_internal
A boost sink device which counts the number of characters written
*/
struct char_counting_sink {
  char_counting_sink(size_t initial = 0):count(initial) { }
  
  char_counting_sink(const char_counting_sink &buf):count(buf.count) { }
  size_t count;
  typedef char        char_type;
  struct category: public boost::iostreams::device_tag,
                       public boost::iostreams::output,
                       public boost::iostreams::multichar_tag,
                       public boost::iostreams::optimally_buffered_tag { };

 /** the optimal buffer size is 0. */
  inline std::streamsize optimal_buffer_size() const { return 0; }

  inline std::streamsize write(const char* s, std::streamsize n) {
    count += n;
    return n;
  }
};
}
#endif
