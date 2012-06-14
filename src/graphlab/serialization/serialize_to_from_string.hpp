#ifndef SERIALIZE_TO_FROM_STRING_HPP
#define SERIALIZE_TO_FROM_STRING_HPP
#include <sstream>
#include <boost/iostreams/stream.hpp>

namespace graphlab {
  /**
   * \ingroup group_serialization
   * \brief Serializes a object to a string
   * 
   * Converts a \ref serializable object t to a string
   * using the serializer.
   * 
   * \tparam T the type of object to serialize. Typically
   *           will be inferred by the compiler. 
   *
   * \param t The object to serializer
   * \returns A string containing a serialized form of t 
   *
   * \see deserialize_from_string()
   */
  template <typename T>
  inline std::string serialize_to_string(const T &t) {
    std::stringstream strm;
    oarchive oarc(strm);
    oarc << t;
    strm.flush();
    return strm.str();
  }


  /**
   * \ingroup group_serialization
   * \brief Deserializes a object from a string
   * 
   * Deserializes a \ref serializable object t from a string
   * using the deserializer.
   * 
   * \tparam T the type of object to deserialize. Typically
   *           will be inferred by the compiler. 
   *
   * \param s The string to deserialize 
   * \param t A reference to the object which will contain 
   *          the deserialized object when the function returns
   *
   * \see serialize_from_string()
   */
  template <typename T>
  inline void deserialize_from_string(const std::string &s, T &t) {
    boost::iostreams::stream<boost::iostreams::array_source> 
          istrm(s.c_str(), s.length());   
    iarchive iarc(istrm);
    iarc >> t;
  }
}

#endif
