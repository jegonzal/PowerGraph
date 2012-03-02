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

/* Obtain a backtrace and print it to stderr. */
void __print_back_trace() {
    void    *array[1024];
    size_t  size, i;
    char    **strings;

    size = backtrace(array, 1024);
    strings = backtrace_symbols(array, size);

    char demangled_name[10240];
    size_t length;
    int status;
    for (i = 0; i < size; ++i) {
        char* ret = abi::__cxa_demangle(strings[i], demangled_name, &length, &status);
        if (ret != NULL) {
          fprintf(stderr, "%s\n", ret);
        }
        else {
          fprintf(stderr, "%s\n", strings[i]);
        }
    }
    free(strings);
}

