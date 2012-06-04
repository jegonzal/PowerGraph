/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef MATWRAP
#define MATWRAP

#include <cstring>
#include <cstdio>

#include "mex.h"



 
/**
 * A safe assertion for matlab
 */
void safe_assert(const bool value, const char* msg) {
  if(!value) { mexErrMsgTxt(msg); }
  //  if( __builtin_expect(!value, 0) ) { mexErrMsgTxt(msg); }
};


void flush_screen() {
  mexEvalString("drawnow");
}


/**
 * A convenient wrapper around mxArray objects
 */
struct matwrap {

  mxArray* array;
  
  matwrap(mxArray* array = NULL) : array(array) {  }



  bool is_null() const { return array == NULL; }

  matwrap get_property(const char* property) const { 
    safe_assert(array != NULL, "dereferenced null mxArray");
    mxArray* result(mxGetProperty(array, 0, property));
    if(result == NULL) {
      char buffer[256];
      sprintf(buffer, "Invalid property %s\n", property);
      mexErrMsgTxt(buffer);
    }
    return matwrap(result);
  } // end of get property

  
  template<typename T>  T* get_data() {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return reinterpret_cast<T*>(mxGetData(array));
  }

  double* get_double_array() {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetPr(array);
  }

  double& mat_index2d(const size_t i, const size_t j) {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetPr(array)[ i + j*rows()];    
  }


  void set_cell_index2d(const size_t i, const size_t j,
                        matwrap contents) {
    safe_assert(array != NULL, "dereferenced null mxArray");
    mxSetCell(array, i + j*rows(), contents.array);    
  }

  matwrap get_cell_index2d(const size_t i, const size_t j) {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetCell(array, i + j*rows());
  }


  matwrap get_field(const char* fieldname) const {
    safe_assert(is_struct(),
                "Attempted to access field of a non-struct element.");
    return mxGetField(array,0,fieldname);
  }

  int get_number_of_fields() const { return  mxGetNumberOfFields(array); }

  matwrap get_field(const int field_id) const {
    safe_assert(is_struct(),
                "Attempted to access field of a non-struct element.");
    safe_assert(field_id < get_number_of_fields(), "Invalid field id!");
    return mxGetFieldByNumber(array,0,field_id);
  }

  int get_field_number(const char* field_name) const {
    safe_assert(is_struct(),
                "Attempted to access field of a non-struct element.");
    return mxGetFieldNumber(array, field_name);
  }

  bool is_class(const char* classname) const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    safe_assert(classname != NULL, "Invalid classname argument");
    return mxIsClass(array, classname);     
  } // end of is class

  const char* get_classname() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetClassName(array);
  } // end of is class


  size_t size() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetNumberOfElements(array);
  }

  matwrap get_cell(size_t index) const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    safe_assert(mxGetClassID(array) == mxCELL_CLASS,
                "Attempted to access a cell in a non-cell array.");
    return mxGetCell(array, index);
  }  

  bool is_cell() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxIsCell(array);
  }

  bool is_struct() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxIsStruct(array);
  }

  bool is_double() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetClassID(array) == mxDOUBLE_CLASS;
  }


  bool is_uint32() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetClassID(array) == mxUINT32_CLASS;
  }

  bool is_string() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetClassID(array) == mxCHAR_CLASS;
  }


  void as_string(char* str_buffer, size_t buffer_len) {
    safe_assert(str_buffer != NULL, "NULL string buffer");
    int error = mxGetString(array, str_buffer, buffer_len);
    safe_assert(!error, "Error processing string!");
  }

  const mwSize* get_dimensions() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetDimensions(array);
  }

  size_t rows() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetM(array);
  }

  size_t cols() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetN(array);
  }


  const mwSize get_num_dimensions() const {
    safe_assert(array != NULL, "dereferenced null mxArray");
    return mxGetNumberOfDimensions(array);

  }











  static matwrap create_matrix(size_t m, size_t n) {
    return matwrap(mxCreateDoubleMatrix(m, n, mxREAL));
  }

  static matwrap create_cell(size_t m, size_t n) {
    return matwrap(mxCreateCellMatrix(m, n));
  }




  
};



#endif
