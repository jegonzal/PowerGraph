#ifndef HAS_SAVE_HPP
#define HAS_SAVE_HPP
#include <typeinfo>
/**
Detects if a class has a save function implemented
*/
namespace graphlab {
namespace archive_detail {

 // SFINAE method derived from
 // http://stackoverflow.com/questions/87372/is-there-a-technique-in-c-to-know-if-a-class-has-a-member-function-of-a-given-s/87846#87846
  template<typename ArcType, typename T>
  struct has_save_method
  {
	  template<typename U, void (U::*)(ArcType&) const> struct SFINAE {};
	  template<typename U> static char Test(SFINAE<U, &U::save>*);
	  template<typename U> static int Test(...);
	  static const bool value = sizeof(Test<T>(0)) == sizeof(char);
  };


  template <typename ArcType, typename ValueType>
  typename boost::enable_if_c<has_save_method<ArcType, ValueType>::value, void>::type 
  save_or_fail(ArcType& o, const ValueType &t) { 
    t.save(o);
  }
  
  template <typename ArcType, typename ValueType>
  typename boost::disable_if_c<has_save_method<ArcType, ValueType>::value, void>::type 
  save_or_fail(ArcType& o, const ValueType &t) { 
    ASSERT_MSG(false,"Trying to serializable type %s without valid save method.", typeid(ValueType).name()); 
  }
 
}  // archive_detail
}  // graphlab

#endif
