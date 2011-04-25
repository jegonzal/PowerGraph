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

    static __any_placeholder* base_load(iarchive_soft_fail &arc);

    virtual void deep_op_equal(__any_placeholder* c) = 0;

    void base_save(oarchive_soft_fail &arc) const;

    virtual void load(iarchive_soft_fail &arc) = 0;
    virtual void save(oarchive_soft_fail &arc) const = 0;
    virtual std::ostream& print(std::ostream& out) const = 0;

  };

  typedef __any_placeholder* (*__any_registration_deserializer_type)(iarchive_soft_fail &arc);
  typedef std::map<uint64_t, __any_registration_deserializer_type>
  __any_registration_map_type;

  __any_registration_map_type& __get_registration_map();

  /**
   A generic "variant" object obtained from Boost::Any and modified
   to be serializable. A variable of type "any" can store any datatype 
   (even dynamically changeable at runtime), but the caveat is that you must
   know the exact stored type to be able to extract the data safely.
   
   To serialize/deserialize the any, regular serialization procedures apply.
   However, since a statically initialized type registration system is used to 
   identify the type of the deserialized object, so the user must pay attention
   to a couple of minor issues.
   
   On serialization: 
   
   \li \b a) If an any contains a serializable type, the any can be serialized.
   \li \b b) If an any contains an unserializable type, the serialization will fail at runtime.
   
   On deserialization:
   
   \li \b c) An empty any can be constructed with no type information and it can be
             deserialized from an archive. 
   \li \b d) However, the deserialization will fail at runtime if the true type of the any
             is never accessed / instantiated anywhere in the code.
             
   Condition \b d) is particular unusual so I will illustrate with an example.
   
   Given a simple user struct:
   \code
   struct UserStruct {
     int i;
     void save (graphlab::oarchive& oarc) const {
       oarc << i;
     }
     void load (graphlab::iarchive& iarc) {
       iarc << i;
     }
   }
  \endcode
  
   If an any object contains the struct, it will be serializable.
   
   \code
   UserStruct us;
   us.i = 10;
   any a = us;
   // output file
   std::ofstream fout("test.bin");
   graphlab::oarchive oarc(fout);
   oarc << a;    // write the any
   \endcode
   
   To deserialize, I will open an input archive and stream into an any.
   
   \code
   // open input 
   std::ifstream fin("test.bin");
   graphlab::iarchive iarc(fin);
   // create an any and read it
   any a;
   iarc >> a; 
   \endcode
   
   Now, unusually, the above code will fail, while the following code
   will succeed
   
   \code
   // open input 
   std::ifstream fin("test.bin");
   graphlab::iarchive iarc(fin);
   // create an any and read it
   any a;
   iarc >> a; 
   std::cout << a.as<UserStruct>().i;
   \endcode
   
   The <tt> a.as<UserStruct>() </tt> forces the instantiation of static functions
   which allow the any deserialization to identify the UserStruct type.
  */
  class any {
  public: // structors
    /// default constructor. Creates an empty any
    any() : content(NULL) { }

    /// Creates an any which stores the value
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

    /// Returns true if the object does not contain any stored data
    bool empty() const {
      return content == NULL;
    }

    /// Returns the type information of the stored data.
    const std::type_info & type() const {
      return empty() ? typeid(void) : content->type();
    }

    /// loads the any from a file.
    void load(iarchive &arc) {
      iarchive_soft_fail isoftarc(arc);
      if(content != NULL) delete content;
      bool isempty;
      isoftarc >> isempty;
      if (isempty == false) {
        content = __any_placeholder::base_load(isoftarc);
      }
    }
    
    /// Saves the any to a file. Caveats apply. See the main any docs.
    void save(oarchive &arc) const {
      oarchive_soft_fail osoftarc(arc);
      bool isempty = empty();
      osoftarc << isempty;
      if (isempty == false) content->base_save(osoftarc);
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
      static __any_placeholder* deserialize(iarchive_soft_fail &arc) {
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

      holder(iarchive_soft_fail &arc) { arc >> held; }

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

      void load(iarchive_soft_fail &arc) {
        arc >> held;
      }

      void save(oarchive_soft_fail &arc) const {
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
    /// Extracts a reference to the contents of the any as a type of ValueType
    template<typename ValueType>
    ValueType& as() {
      // force instantiation of the registration template (without running the constructor)
      any::holder<ValueType>::__registration.get_deserializer_id();
      DASSERT_TRUE(type() == typeid(ValueType));
      DASSERT_FALSE(empty());
      return static_cast<any::holder<ValueType> *>(content)->held;
    }

    /// Extracts a constant reference to the contents of the any as a type of ValueType
    template<typename ValueType>
    inline const ValueType& as() const{
      // force instantiation of the registration template (without running the constructor)
      any::holder<ValueType>::__registration.get_deserializer_id();
      DASSERT_TRUE(type() == typeid(ValueType));
      DASSERT_FALSE(empty());
      return static_cast<any::holder<ValueType> *>(content)->held;
    }

    /// Exchanges the contents of two any's
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




} // namespace graphlab

// std::ostream& operator<<(std::ostream& out, const graphlab::any& any);


// Copyright Kevlin Henney, 2000, 2001, 2002. All rights reserved.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#endif

