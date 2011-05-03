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

#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Obtain a backtrace and print it to stderr. */
void __print_back_trace() {
  const size_t array_size(1024);
  void *array[array_size];
  size_t size;
  
  size = backtrace(array, array_size);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  // backtrace_symbols_fd(array, size, STDOUT_FILENO);
}
