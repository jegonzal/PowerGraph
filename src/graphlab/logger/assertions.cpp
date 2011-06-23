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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Obtain a backtrace and print it to stderr. */
void __print_back_trace() {
  const size_t array_size(1024);
  void *array[array_size];
  int size;
  
  size = backtrace(array, array_size);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  // backtrace_symbols_fd(array, size, STDOUT_FILENO);
}

