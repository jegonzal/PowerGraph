#ifndef GRAPHLAB_OARCHIVE_HPP
#define GRAPHLAB_OARCHIVE_HPP

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
  
  class oarchive{
   public:
    std::ostream* o;

    oarchive(std::ostream& os)
      : o(&os) {}

    ~oarchive() {
      if (o)
        o->flush();
    }
  };

  /** Serializes a single character. 
      Assertion fault on failure. */
  oarchive& operator<<(oarchive& a, const char i);

  /** Serializes a floating point number. 
      Assertion fault on failure. */
  oarchive& operator<<(oarchive& a, const float i);

  /** Serializes a double precisition floating point number. 
      Assertion fault on failure. */
  oarchive& operator<<(oarchive& a, const double i);

  oarchive& operator<<(oarchive& a, const bool i);
  oarchive& operator<<(oarchive& a, const unsigned char i);
  
  
  /** Serializes a integer. 
      Assertion fault on failure. */
  oarchive& operator<<(oarchive& a, const short i);
  oarchive& operator<<(oarchive& a, const unsigned short i);
  oarchive& operator<<(oarchive& a, const int i);
  oarchive& operator<<(oarchive& a, const long i);
  oarchive& operator<<(oarchive& a, const long long i);
  oarchive& operator<<(oarchive& a, const unsigned long i);
  oarchive& operator<<(oarchive& a, const unsigned int i);
  oarchive& operator<<(oarchive& a, const unsigned long long  i);
  oarchive& serialize_64bit_integer(oarchive& a, const int64_t i);


  /** Serializes a generic pointer object. bytes from (i) to (i + length - 1) 
      inclusive will be written to the file stream. The length will also be
      Assertion fault on failure.  */
  oarchive& serialize(oarchive& a, const void* i,const size_t length);

  /** Serializes a C string 
      Assertion fault on failure. */
  oarchive& operator<<(oarchive& a, const char* s);

  /** Serializes a string. 
      Assertion fault on failure.  */
  oarchive& operator<<(oarchive  &a, const std::string& s);

  oarchive& operator<<(oarchive& a, void* const t);
  
  /** Serializes a pair
      Assertion fault on failure.   */
  template <typename T,typename U>
  oarchive& operator<<(oarchive& a, const std::pair<T,U>& p) {
    a << p.first;
    a << p.second;
    return a;
  }
  
 
namespace oarchive_detail {

 // SFINAE method derived from
 // http://stackoverflow.com/questions/87372/is-there-a-technique-in-c-to-know-if-a-class-has-a-member-function-of-a-given-s/87846#87846
  template<typename T>
  struct has_save_method
  {
	  template<typename U, void (U::*)(oarchive&) const> struct SFINAE {};
	  template<typename U> static char Test(SFINAE<U, &U::save>*);
	  template<typename U> static int Test(...);
	  static const bool value = sizeof(Test<T>(0)) == sizeof(char);
  };


  template <typename ValueType>
  typename boost::enable_if_c<has_save_method<ValueType>::value, void>::type 
  save_or_fail(oarchive& o, const ValueType &t) { 
    t.save(o);
  }
  
  template <typename ValueType>
  typename boost::disable_if_c<has_save_method<ValueType>::value, void>::type 
  save_or_fail(oarchive& o, const ValueType &t) { 
    ASSERT_MSG(false,"Trying to serializable type %s without valid save method.", typeid(ValueType).name()); 
  }
}
 
 

  /** catch all operator<< as member of iarchive */
  template <typename T>
  inline oarchive& operator<<(oarchive& a, const T& t) {
    oarchive_detail::save_or_fail(a, t);
    return a;
  }
} // namespace prl

#endif  //PRL_OARCHIVE_HPP
