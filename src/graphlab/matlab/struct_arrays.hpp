#ifndef mex_h
#error "mex.h must be included before this file"
#endif

#ifndef STRUCT_ARRAYS_HPP
#define STRUCT_ARRAYS_HPP
#include <iostream>

template<class T, class EMXType>
bool clear_struct_array(EMXType &out) {
  out.data = (T*)malloc(sizeof(T) * 1);
  converter<T>::clearemx(out.data[0]);
  out.size = (int32_t*)malloc(sizeof(int32_t) * 2);
  out.size[0] = 1; out.size[1] = 1;
  out.allocatedSize = 1;
  out.numDimensions = 2;
  out.canFreeData = 1;
  return true;
}


template<class T, class EMXType>
void free_struct_array(EMXType &in) {
  if (!in.canFreeData) return;
  for (size_t i = 0 ;i < in.allocatedSize; ++i) {
    freeemx(in.data[i]);
  }
  free(in.data);
  free(in.size);
  in.canFreeData = 0;
}



/**
 * Converts an mxArray to an emxArray(EMXType ) of a particular struct type (T)
 * out MUST be an unallocated emxArray structure. 
 * The annoying part about struct arrays is that I am unable to slice an
 * struct mxArray easily.
 * Returns true on success, false on failure 
 */
template<class T, class EMXType>
bool read_struct_array(const mxArray *array, EMXType &out) {
  // make sure the array is the right type
  if (array == NULL) {
    clear_struct_array<T, EMXType>(out);
    return true;
  }
  
  static bool printed = false;
  if (mxGetClassID(array) != mxSTRUCT_CLASS) {
    if (!printed) std::cerr << "Trying to convert a non-struct array into a struct" << std::endl;
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
  // perform the copy
  // now this is annoyingly complicated since I can't slice the mxarray easily
  // get the field names
  size_t numfields = mxGetNumberOfFields(array);
  const char* fieldnames[numfields];
  for (size_t i = 0;i < numfields; ++i) {
    fieldnames[i] = mxGetFieldNameByNumber(array, i);
  }
  // create a temporary struct array. this will store the slice
  mwSize dims[1];
  dims[0] = 1;
  mxArray* tmp = mxCreateStructArray(1,dims,numfields,fieldnames);
  // acceleration technique: So that we can use shallow copies
  // remember the old mxarrays
  mxArray* oldmxarrays[numfields];
  for (size_t j = 0;j < numfields; ++j) {
    oldmxarrays[j] = mxGetFieldByNumber(tmp,0,j);
  }
  
  // loop through all the elements
  for (size_t i = 0;i < numel; ++i) {
    // slice
    // we use a shallow slice here. Passing pointers to the inner mxarrays
    for (size_t j = 0;j < numfields; ++j) {
      mxSetFieldByNumber(tmp, 0, j, mxGetFieldByNumber(array,i,j));
    }
    // store
    converter<T>::mxarray2emx(tmp, out.data[i]);
  }
  // before we destroy, we restore the old mxarrays
  for (size_t j = 0;j < numfields; ++j) {
    mxSetFieldByNumber(tmp,0,j, oldmxarrays[j]);
  }
  
  mxDestroyArray(tmp);
  return true;
}

/**
 * Converts an emxArray(EMXType ) of a particular struct type T to an mxArray
 * array must not be allocated.
 * Returns true on success, false on failure 
 */
template<class T, class EMXType>
bool write_struct_array(EMXType &in, mxArray * &array) {
  static bool printed = false;
  // empty array
  if (in.numDimensions == 0) {
    array = NULL;
    return true;
  }
  // ok. I need a sample of a conversion so I can get the field names
  mxArray* tmp = NULL;
  converter<T>::emx2mxarray(in.data[0], tmp);
  
  // get the field names
  size_t numfields = mxGetNumberOfFields(tmp);
  const char* fieldnames[numfields];
  for (size_t i = 0;i < numfields; ++i) {
    fieldnames[i] = mxGetFieldNameByNumber(tmp, i);
  }
  mxDestroyArray(tmp);
  // get the dimensions. 
  mwSize dimptr[in.numDimensions];
  size_t numel = 1;
  for (size_t i = 0; i < in.numDimensions; ++i) {
    dimptr[i] = in.size[i];
    numel = numel * in.size[i];
  }
  // construct the array
  array = mxCreateStructArray(in.numDimensions,dimptr,numfields,fieldnames);
  
  // store the data
  
  // loop through all the elements
  for (size_t i = 0;i < numel; ++i) {
    // slice
    // read an element
    tmp = NULL;
    converter<T>::emx2mxarray(in.data[i], tmp);
    for (size_t j = 0;j < numfields; ++j) {
        // destroy the old mxarray here if any
        mxArray* oldmxarray = mxGetFieldByNumber(array,i,j);
        if (oldmxarray != NULL) mxDestroyArray(tmp);
        mxSetFieldByNumber(array, i, j, mxDuplicateArray(mxGetFieldByNumber(tmp,0,j)));
    }
    mxDestroyArray(tmp);
  }
  return true;
}
#endif  