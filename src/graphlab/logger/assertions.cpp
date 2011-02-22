#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Obtain a backtrace and print it to stderr. */
void __print_back_trace() {
  void *array[200];
  size_t size;
  
  size = backtrace(array, 200);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
}
