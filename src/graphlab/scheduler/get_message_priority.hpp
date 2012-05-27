#ifndef GRAPHLAB_SCHEDULER_GET_MESSAGE_PRIORITY_HPP
#define GRAPHLAB_SCHEDULER_GET_MESSAGE_PRIORITY_HPP

#include <boost/type_traits.hpp>
#include <typeinfo>

namespace graphlab {
  
namespace scheduler_impl {
  
template <typename T>
class implements_priority_member {
  template<typename U, double (U::*)() const> struct SFINAE {};
  template <typename U> char test(SFINAE<U, &U::priority>*);
  template <typename U> int test(...);
  static const bool value = sizeof(test<T>(0)) == sizeof(char);
};

template <typename MessageType>
typename boost::enable_if_c<implements_priority_member<MessageType>::value, double>::type 
get_message_priority(const MessageType &m) { 
  return m.priority();
}

template <typename MessageType>
typename boost::disable_if_c<implements_priority_member<MessageType>::value, double>::type 
get_message_priority(const MessageType &m) { 
  return 1.0;
}

} //namespace scheduler_impl
} //namespace graphlab


#endif