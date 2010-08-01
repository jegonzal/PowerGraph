#include <iostream>
#include <graphlab/serialization/iarchive.hpp>

namespace graphlab {

  /***************************************************************************
   *                        Basic Deserializers                              * 
   ***************************************************************************/

  // generate the serialize call for a whole bunch of integer types
#define GENCASTDESERIALIZE(typesrc, typedest)           \
  iarchive& operator>>(iarchive &a, typesrc &i) {       \
    typedest c;                                         \
    operator>>(a, c);                                   \
    i = static_cast<typesrc>(c);                        \
    return a;                                           \
  } 

#define GENINTDESERIALIZE(typesrc)                      \
  iarchive& operator>>(iarchive &a, typesrc &i) {       \
    int64_t c;                                          \
    deserialize_64bit_integer(a, c);                    \
    i = static_cast<typesrc>(c);                        \
    return a;                                           \
  } 
  
  iarchive& operator>>(iarchive& a, char& i) {
    a.i->get(i);
    assert(!a.i->fail());
    return a;
  }

  iarchive& deserialize_64bit_integer(iarchive& a, int64_t& i) {
    a.i->read(reinterpret_cast<char*>(&i), sizeof(int64_t));
    assert(!a.i->fail());
    return a;
  }


  GENCASTDESERIALIZE(bool, char);
  GENCASTDESERIALIZE(unsigned char, char);
  // serializes a bunch of integer types
  GENINTDESERIALIZE(short);
  GENINTDESERIALIZE(unsigned short);
  GENINTDESERIALIZE(int);
  GENINTDESERIALIZE(long);
  GENINTDESERIALIZE(long long);
  GENINTDESERIALIZE(unsigned long);
  GENINTDESERIALIZE(unsigned int);
  GENINTDESERIALIZE(unsigned long long);


  iarchive& operator>>(iarchive& a, float& i) {
    a.i->read(reinterpret_cast<char*>(&i), sizeof(float));
    assert(!a.i->fail());
    return a;
  }

  iarchive& operator>>(iarchive& a, double& i) {
    a.i->read(reinterpret_cast<char*>(&i), sizeof(double));
    assert(!a.i->fail());
    return a;
  }

  iarchive& deserialize(iarchive& a, void* const i,const size_t length) {
    // Save the length and check if lengths match
    size_t length2;
    operator>>(a, length2);
    assert(length == length2);

    //operator>> the rest
    a.i->read(reinterpret_cast<char*>(i), length);
    assert(!a.i->fail());
    return a;
  }

  iarchive& operator>>(iarchive& a, char*& s) {
    // Save the length and check if lengths match
    size_t length;
    operator>>(a, length);
    if (s == NULL) s = new char[length+1];
    //operator>> the rest
    a.i->read(reinterpret_cast<char*>(s), length);
    s[length] = 0;
    assert(!a.i->fail());
    return a;
  }

  iarchive& operator>>(iarchive& a, std::string& s){
    //read the length
    size_t length;
    a >> length;
    //resize the string and read the characters
    s.resize(length);
    a.i->read(const_cast<char*>(s.c_str()), length);
    assert(!a.i->fail());
    return a;
  }


  iarchive& operator>>(iarchive& a, void* &t) {
    size_t tmp;
    a >> tmp;
    t = (void*)tmp;
    return a;
  }
} // namespace prl

