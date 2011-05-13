/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PODIFY_HPP
#define PODIFY_HPP
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
namespace graphlab {

 /**
This is a macro which allows you to temporarily define
an object as a POD. 
for instance if you have a struct that has no save/load
function defined, but you want to treat it as a POD.
Now you can write:

SOmePODStruct s;
oarc << PODIFY(s);
and
iarc >> PODIFY(s);
*/

#define PODIFY(X) (__podify__<typeof(X)>(X))

/**
This wraps and adds temporary support for serialization
of a POD object
*/
template <typename T>
class __podify__{
 public:
   
  T alt;
  const T* ref;
  // if created without constructor, assign 
  // the pointer to the built-in version
  __podify__():ref(&alt) { }
  
  __podify__(const T &a):ref(&a) { }
  
  operator const T&() const{
    return *ref;
  }
  
  operator T&() {
    return *(const_cast<T*>(ref));
  }
  
  void save(oarchive &a) const{
    serialize(a, ref, sizeof(T));
  }
  // unfortunately load will not work 
  // because we this class is only created as a temporary
  // we are unable to get a reference to it
  // so we make the pass by value version below
};


template <typename T>
iarchive& operator>>(iarchive& a, __podify__<T> i) {
    T* refptr = (const_cast<T*>(i.ref));
    deserialize(a, refptr, sizeof(T));
    return a;
}

}
#endif

