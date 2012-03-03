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


#include <execinfo.h>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <cxxabi.h>


/** Code from http://mykospark.net/2009/09/runtime-backtrace-in-c-with-name-demangling/ */
std::string demangle(const char* symbol) {
  size_t size;
  int status;
  char temp[1024];
  char* demangled;
  //first, try to demangle a c++ name
  if (1 == sscanf(symbol, "%*[^(]%*[^_]%127[^)+]", temp)) {
    if (NULL != (demangled = abi::__cxa_demangle(temp, NULL, &size, &status))) {
      std::string result(demangled);
      free(demangled);
      return result;
    }
  }
  //if that didn't work, try to get a regular c symbol
  if (1 == sscanf(symbol, "%127s", temp)) {
    return temp;
  }
 
  //if all else fails, just return the symbol
  return symbol;
}

/* Obtain a backtrace and print it to stderr. */
void __print_back_trace() {
    void    *array[1024];
    size_t  size, i;
    char    **strings;

    size = backtrace(array, 1024);
    strings = backtrace_symbols(array, size);

    fprintf(stderr, "Pointers\n");
    fprintf(stderr, "------------\n");
    for (i = 0; i < size; ++i) {
        fprintf(stderr, "%p\n", array[i]);
    }
 

    fprintf(stderr, "Raw\n");
    fprintf(stderr, "------------\n");
    for (i = 0; i < size; ++i) {
        fprintf(stderr, "%s\n", strings[i]);
    }
    fprintf(stderr, "\nDemangled\n");
    fprintf(stderr, "------------\n");
 
    for (i = 0; i < size; ++i) {
        std::string ret = demangle(strings[i]);
        fprintf(stderr, "%s\n", ret.c_str());
    }
    free(strings);
}

