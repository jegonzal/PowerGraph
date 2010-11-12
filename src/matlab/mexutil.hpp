#ifndef MEXUTIL_HPP
#define MEXUTIL_HPP
#include <iostream>
#include <typeinfo>
#include <cstring>
#include "mex_classid_traits.hpp"

template<class T, class EMXType>
bool clear_array(EMXType &out) {
  out.data = (T*)malloc(sizeof(T) * 1);
  out.size = (int32_t*)malloc(sizeof(int32_t) * 2);
  out.size[0] = 1; out.size[1] = 1;
  out.allocatedSize = 1;
  out.numDimensions = 2;
  out.canFreeData = 1;
}

/**
 * Converts an mxArray to an emxArray(EMXType ) of a particular type (T)
 * out MUST be an unallocated emxArray structure. 
 * Returns true on success, false on failure
 */
template<class T, class EMXType>
bool read_array(const mxArray *array, EMXType &out) {
  // reset out
  static bool printed = false;
  if (array == NULL) {
    clear_array<T, EMXType>(out);
    return true;
  }
  // check for type compatibility
  // element sizes must be the same
  // and it must pass the class compatibility check
  if (mxGetElementSize(array) != sizeof(T) ||
     compatible_classid<T>(mxGetClassID(array)) == false) {
    if (!printed) {
      std::cerr << "Parse type incompatibility when reading mxArray." << std::endl;
      std::cerr << "Failed conversion from " << mxGetClassName(array) << " to " << typeid(T).name() << std::endl;
    }
    printed = true;
    return false;
  }
  // get the number of dimensions
  out.numDimensions = mxGetNumberOfDimensions(array);
  // get the dimensions. Allocate the size array and copy it
  // we have to do this element by element since mwSize may not be the same
  // type as int32_t
  const mwSize* dimptr = mxGetDimensions(array);
  out.size = (int32_t *)malloc((sizeof(int32_t) * out.numDimensions));
  // print if allocation failed
  if (out.size == NULL) {
    if (!printed) std::cerr << "Malloc Failure allocating " << out.numDimensions << " * 4" << std::endl;
    printed = true;
    return false;
  }
  // count the number of elements
  size_t numel = 1;
  for (size_t i = 0; i < out.numDimensions; ++i) {
    out.size[i] = dimptr[i];
    numel = numel * dimptr[i];
  }
  
  // fill in the allocated size and allocate the data
  out.allocatedSize = numel;
  out.data = (T*)malloc((sizeof(T) * numel));
  out.canFreeData = 1;
  // print if allocation failed
  if (out.data== NULL) {
    if (!printed) std::cerr << "Malloc Failure allocating " << numel << " * " << sizeof(T) << std::endl;
    printed = true;
    return false;
  }
  // copy the data. Since we have guaranteed type storage compatibility, 
  // we can do a direct memcpy
  void* ptr = NULL;
  if (mxIsComplex(array)) {
    ptr = mxGetImagData(array);
  }
  else {
    ptr = mxGetData(array);
  }
  memcpy(out.data, ptr, sizeof(T) * numel);
  return true;
}










/**
 * Converts an emxArray(EMXType ) of a particular type T to an mxArray
 * array must not be allocated.
 * Returns true on success, false on failure
 */
template<class T, class EMXType>
bool write_array(EMXType &in, mxArray * &array) {
  static bool printed = false;
  // empty array
  if (in.numDimensions == 0) {
    array = NULL;
    return true;
  }
  
  // get the dimensions. 
  mwSize dimptr[in.numDimensions];
  for (size_t i = 0; i < in.numDimensions; ++i) dimptr[i] = in.size[i];
  
  
  // check if I can get a valid class id from the emx type
  if (prefered_classid<T>::cid() == mxUNKNOWN_CLASS) {
    if (!printed) {
      std::cerr << "Unrecognized type " << typeid(T).name() << std::endl;
      std::cerr << "Bad! Notify the developers!" << std::endl;
    }
    printed = true;
    return false;
  }
  
  // allocate the array
  array = mxCreateNumericArray(in.numDimensions, 
                               dimptr, 
                               prefered_classid<T>::cid(), 
                               prefered_classid<T>::complex()?mxCOMPLEX:mxREAL);
  
  // check for type storage compatibility
  if (mxGetElementSize(array) != sizeof(T)) {
    if (!printed) {
      std::cerr << "Failed conversion from " << typeid(T).name() << " to " << mxGetClassName(array) << std::endl;
      std::cerr << "Bad! Notify the developers!" << std::endl;
    }
    printed = true;
    return false;
  }
  // copy the data. Since we have guaranteed type storage compatibility, 
  // we can do a direct memcpy
  
  void* ptr = NULL;
  if (prefered_classid<T>::complex()) {
    ptr = mxGetImagData(array);
  }
  else {
    ptr = mxGetData(array);
  }
  memcpy(ptr, in.data, sizeof(T) * mxGetNumberOfElements(array));
  return true;
}

/**
 * If the mxArray is a struct and has a field, this is set to the field's mxArray
 * Returns NULL otherwise
 */
inline const mxArray* struct_has_field(const mxArray *array, const char *fieldname) {
  if (mxGetClassID(array) == mxSTRUCT_CLASS) {
      return mxGetField(array,0,fieldname);
  }
  return NULL;
}
#endif