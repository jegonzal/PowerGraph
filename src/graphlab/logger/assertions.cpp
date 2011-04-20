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
