#ifndef GRAPHLAB_CHARSTREAM
#define GRAPHLAB_CHARSTREAM

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/categories.hpp>

namespace graphlab {

  /// \ingroup util_internal
  namespace charstream_impl {
    /// \ingroup util_internal
    template <bool self_deleting>
    struct resizing_array_sink {


      resizing_array_sink(size_t initial = 0) : str(NULL) { 
        if(initial > 0) {
          str = (char*)(malloc(initial));
          assert(str != NULL);
        } 
        len = 0;
        buffer_size = initial;
      }

      resizing_array_sink(const resizing_array_sink& other) :
        len(other.len), buffer_size(other.buffer_size) {
        if(self_deleting) {
          str = (char*)(malloc(other.buffer_size));
          assert(str != NULL);
          memcpy(str, other.str, len);
        } else {
          str = other.str;
        }
      }

      ~resizing_array_sink() {
        if( self_deleting && str != NULL) {
          free(str);
          str = NULL;
        }        
      }


      size_t size() const { return len; }
      char* c_str() { return str; }

      void clear() {
        len = 0;
      }

      char *str;
      size_t len;
      size_t buffer_size;
      typedef char        char_type;
      struct category: public boost::iostreams::device_tag,
                       public boost::iostreams::output,
                       public boost::iostreams::multichar_tag { };
      
      /** the optimal buffer size is 0. */
      inline std::streamsize optimal_buffer_size() const { return 0; }
      
      inline std::streamsize write(const char* s, std::streamsize n) {
        if (len + n > buffer_size) {
          // double in length if we need more buffer
          buffer_size = 2 * (len + n);
          str = (char*)realloc(str, buffer_size);
          assert(str != NULL);
        }
        memcpy(str + len, s, n);
        len += n;
        return n;
      }
    };

  }; // end of impl;

  typedef boost::iostreams::stream< charstream_impl::resizing_array_sink<true> > 
  charstream;


}; // end of namespace graphlab
#endif
