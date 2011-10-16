#ifndef SERIALIZE_TO_FROM_STRING_HPP
#define SERIALIZE_TO_FROM_STRING_HPP
#include <sstream>
#include <boost/iostreams/stream.hpp>

namespace graphlab {
  template <typename T>
  inline std::string serialize_to_string(const T &t) {
    std::stringstream strm;
    oarchive oarc(strm);
    oarc << t;
    strm.flush();
    return strm.str();
  }
  
  template <typename T>
  inline void deserialize_from_string(const std::string &s, T &t) {
    boost::iostreams::stream<boost::iostreams::array_source> 
          istrm(s.c_str(), s.length());   
    iarchive iarc(istrm);
    iarc >> t;
  }
}

#endif
