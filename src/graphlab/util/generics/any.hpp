// Modified from boost 1.37 boost::any
// Extended to handle graphlab graphlab/serialization/deserialization functions
// See http://www.boost.org/libs/any for Documentation.

#ifndef GRAPHLAB_ANY_INCLUDED
#define GRAPHLAB_ANY_INCLUDED

// what:  variant type boost::any
// who:   contributed by Kevlin Henney,
//        with features contributed and bugs found by
//        Ed Brey, Mark Rodgers, Peter Dimov, and James Curran
// when:  July 2001
// where: tested with BCC 5.5, MSVC 6.0, and g++ 2.95

#include <algorithm>
#include <typeinfo>
#include <map>
#include <iostream>
#include <stdint.h>

#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/throw_exception.hpp>
#include <boost/static_assert.hpp>
#include <boost/utility.hpp>
#include <boost/exception/detail/is_output_streamable.hpp>

#include <graphlab/logger/assertions.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
  
  namespace any_detail {
    template <typename ValueType>
    typename boost::disable_if_c<boost::is_output_streamable<ValueType>::value, void>::type 
    print_type_or_content(std::ostream& out, const ValueType &h) { 
      out << "Type " << typeid(ValueType).name() ; 
    }

    template <typename ValueType>
    typename boost::enable_if_c<boost::is_output_streamable<ValueType>::value, void>::type 
    print_type_or_content(std::ostream& out, const ValueType &h) { 
      out << h; 
    }
  }

  class __any_placeholder {
  public: // structors

    virtual ~__any_placeholder() { }

  public: // queries

    virtual const std::type_info & type() const = 0;

    virtual __any_placeholder * clone() const = 0;

    virtual uint64_t get_deserializer_id() const = 0;

    static __any_placeholder* base_load(iarchive &arc);

    virtual void deep_op_equal(__any_placeholder* c) = 0;

    void base_save(oarchive &arc) const;

    virtual void load(iarchive &arc) = 0;
    virtual void save(oarchive &arc) const = 0;
    virtual std::ostream& print(std::ostream& out) const = 0;

  };

  typedef __any_placeholder* (*__any_registration_deserializer_type)(iarchive &arc);
  typedef std::map<uint64_t, __any_registration_deserializer_type>
  __any_registration_map_type;

  __any_registration_map_type& __get_registration_map();


  class any {
  public: // structors

    any() : content(NULL) { }

    template<typename ValueType>
    any(const ValueType & value)
      : content(new holder<ValueType>(value)) {
      // force instantiation of the registration template (without running the constructor)
      any::holder<ValueType>::__registration.get_deserializer_id();
    }

    any(const any & other)
      : content(other.empty() ? NULL : other.content->clone()) { }

    ~any() {
      delete content;
    }

  public: // queries

    bool empty() const {
      return content == NULL;
    }

    const std::type_info & type() const {
      return empty() ? typeid(void) : content->type();
    }

    // serialization
    void load(iarchive &arc) {
      if(content != NULL) delete content;
      bool isempty;
      arc >> isempty;
      if (isempty == false) {
        content = __any_placeholder::base_load(arc);
      }
    }
    void save(oarchive &arc) const {
      bool isempty = empty();
      arc << isempty;
      if (isempty == false) content->base_save(arc);
    }

  private:
    template <typename T>
    class any_registration {
    public:
      any_registration() {
        __get_registration_map()[get_deserializer_id()] =
          any_registration<T>::deserialize;
        //std::cout << "registered " << typeid(T).name() << " to " << get_deserializer_id() << std::endl;
      }
      
      static bool inited; // whether localid has been created
      static uint64_t localid;  // cached id of this type. Avoids rehashing everything I serialize
      
      static uint64_t get_deserializer_id() {
        if (inited == false) compute_deserializer_id();
        return localid;
      }
      
      static void compute_deserializer_id() {
        // FNV hash function
        const char *p = typeid(T).name();
        uint64_t h = 2166136261u;
        while((*p) != 0) {
          h = ( h * 16777619 ) ^ ((unsigned char)(*p));
          ++p;
        }
        inited = true;
        localid = h;
      }
      static __any_placeholder* deserialize(iarchive &arc) {
        any::holder<T> *newholder = new any::holder<T>(arc);
        return newholder;
      }
    };

    template<typename ValueType>
    class holder : public __any_placeholder {
    public: // structors
      typedef ValueType value_type;
    
      holder(const ValueType & value)
        : held(value) { }

      holder(iarchive &arc) { arc >> held; }

      static any_registration<ValueType> __registration;

    public: // queries

      const std::type_info & type() const {
        return typeid(ValueType);
      }

      __any_placeholder * clone() const {
        return new holder(held);
      }

      void deep_op_equal(__any_placeholder* c) {
        held = static_cast<holder<ValueType>*>(c)->held;
      }

      uint64_t get_deserializer_id() const{
        return __registration.get_deserializer_id();
      }

      void load(iarchive &arc) {
        arc >> held;
      }

      void save(oarchive &arc) const {
        arc << held;
      }
      
      std::ostream& print(std::ostream& out) const {
        any_detail::print_type_or_content(out, held);
        return out;
      }


    public: // representation

      ValueType held;

    private: // intentionally left unimplemented
   
      holder & operator=(const holder &);
    };
  

  public:
    // reading functions
    template<typename ValueType>
    ValueType& as() {
      // force instantiation of the registration template (without running the constructor)
      any::holder<ValueType>::__registration.get_deserializer_id();
      DASSERT_TRUE(type() == typeid(ValueType));
      DASSERT_FALSE(empty());
      return static_cast<any::holder<ValueType> *>(content)->held;
    }

    template<typename ValueType>
    inline const ValueType& as() const{
      // force instantiation of the registration template (without running the constructor)
      any::holder<ValueType>::__registration.get_deserializer_id();
      DASSERT_TRUE(type() == typeid(ValueType));
      DASSERT_FALSE(empty());
      return static_cast<any::holder<ValueType> *>(content)->held;
    }

    // modifiers
    any& swap(any & rhs) {
      std::swap(content, rhs.content);
      return *this;
    }

    template<typename ValueType>
    any& operator=(const ValueType & rhs) {
      // force instantiation of the registration template (without running the constructor)
      any::holder<ValueType>::__registration.get_deserializer_id();
      if (content != NULL && content->type() == typeid(ValueType)) {
        as<ValueType>() = rhs;
      }
      else {
        any(rhs).swap(*this);
      }
      return *this;
    }    

    any& operator=(const any & rhs) {
      if (rhs.empty()) {
        if (content) delete content;
        content = NULL;
      }
      else {
        if (content != NULL && content->type() == rhs.content->type()) {
          content->deep_op_equal(rhs.content);
        }
        else {
          any(rhs).swap(*this);
        }
      }
      return *this;
    }

    std::ostream& print(std::ostream& out) const {
      DASSERT_FALSE(empty());
      return content->print(out);        
    }
    
  private: // representation

    __any_placeholder* content;

  };


  template<typename T> any::any_registration<T> any::holder<T>::__registration;
  template<typename T> bool any::any_registration<T>::inited = false; 
  template<typename T> uint64_t any::any_registration<T>::localid; 
  //template<typename T> size_t any::any_registration<T>::unused;

}



// std::ostream& operator<<(std::ostream& out, const graphlab::any& any);


// Copyright Kevlin Henney, 2000, 2001, 2002. All rights reserved.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#endif
