#ifndef GRAPHLAB_IARCHIVE_HPP
#define GRAPHLAB_IARCHIVE_HPP

#include <iostream>
#include <cassert>
#include <string>
#include <utility>
#include <stdint.h>
#include <typeinfo>
#include <boost/static_assert.hpp>
#include <boost/utility.hpp>
#include <graphlab/logger/logger.hpp>

#ifdef _MSC_VER
#include <itpp/base/ittypes.h> // for int32_t etc.
#endif

namespace graphlab {

  class iarchive {
  public:
    std::istream* i;

    iarchive(std::istream& is)
      :i(&is) { }
  };


  /** Deserializes a single character. 
      Assertion fault on failure. */
  iarchive& operator>>(iarchive& a, char& i);

  /** Deserializes a 64 bit integer. 
      Assertion fault on failure. */


  /** Deserializes a boolean. Assertion fault on failure. */
  iarchive& operator>>(iarchive& a, bool& i);

  /** Deserializes a unsigned char.
      Assertion fault on failure. */
  iarchive& operator>>(iarchive& a, unsigned char& i);

  iarchive& operator>>(iarchive& a, short& i);
  iarchive& operator>>(iarchive& a, unsigned short& i);
  iarchive& operator>>(iarchive& a, int& i);
  iarchive& operator>>(iarchive& a, long& i);
  iarchive& operator>>(iarchive& a, long long& i);
  iarchive& operator>>(iarchive& a, unsigned long& i);
  iarchive& operator>>(iarchive& a, unsigned int& i);
  iarchive& operator>>(iarchive& a, unsigned long long& i);
  iarchive& deserialize_64bit_integer(iarchive& a, int64_t& i);


  /** Deserializes a floating point number. 
      Assertion fault on failure. */
  iarchive& operator>>(iarchive& a, float& i);

  /** Serializes a double precisition floating point number. 
      Assertion fault on failure. */
  iarchive& operator>>(iarchive& a, double& i);

  /** Deserializes a generic pointer object of known length.
      This call must match the corresponding serialize call : 
      \see{serialize(std::ostream &o, const iarchive&* i,const int length)}
      The length of the object is read from the file and checked against the 
      length parameter. If they do not match, the function returns with a failure.
      Otherwise, an additional (length) bytes will be read from the file stream
      into (*i). (*i) must contain at least (length) bytes of memory. Otherwise
      there will be a buffer overflow. 
      Assertion fault on failure. */
  iarchive& deserialize(iarchive& a, void* const i,const size_t length);

  /** Loads a C string. If s is NULL, it will allocate it */
  iarchive& operator>>(iarchive& a, char*& s);

  /** Loads a string.
      Assertion fault on failure. */
  iarchive& operator>>(iarchive& a, std::string& s);
  

  /** void*'s */
  iarchive& operator>>(iarchive& a, void* &t);

  /** Deserializes a pair
      Assertion fault on failure.  */
  template <typename T,typename U>
  iarchive& operator>>(iarchive& a, std::pair<T,U>& p){
    a >> p.first;
    a >> p.second;
    return a;
  }



namespace iarchive_detail {

 // SFINAE method derived from
 // http://stackoverflow.com/questions/87372/is-there-a-technique-in-c-to-know-if-a-class-has-a-member-function-of-a-given-s/87846#87846
  template<typename T>
  struct has_load_method
  {
	  template<typename U, void (U::*)(iarchive&)> struct SFINAE {};
	  template<typename U> static char Test(SFINAE<U, &U::load>*);
	  template<typename U> static int Test(...);
	  static const bool value = sizeof(Test<T>(0)) == sizeof(char);
  };


  template <typename ValueType>
  typename boost::enable_if_c<has_load_method<ValueType>::value, void>::type 
  load_or_fail(iarchive& o, ValueType &t) { 
    t.load(o);
  }
  
  template <typename ValueType>
  typename boost::disable_if_c<has_load_method<ValueType>::value, void>::type 
  load_or_fail(iarchive& o, ValueType &t) { 
    ASSERT_MSG(false,"Trying to deserializable type %s without valid load method.", typeid(ValueType).name()); 
  }
}
   

  /** Catch all serializer as member of iarchive */
  template <typename T>
  inline iarchive& operator>>(iarchive& a, T& t) {
    iarchive_detail::load_or_fail(a, t);
    return a;
  }
} // namespace prl

#endif //PRL_IARCHIVE_HPP
