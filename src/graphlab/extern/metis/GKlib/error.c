/*!
\file  error.c
\brief Various error-handling functions

This file contains functions dealing with error reporting and termination

\author George
\date 1/1/2007
\version\verbatim $Id: error.c 1065 2007-03-06 23:11:52Z karypis $ \endverbatim
*/


#define _GK_ERROR_C_  /* this is needed to properly declare the gk_return_to_entry
                         as an extern function in GKlib.h */

#include <GKlib.h>



/* This is the jmp_buf for the graceful exit in case of severy error */
#ifdef __linux__
__thread jmp_buf gk_return_to_entry;
#else
jmp_buf gk_return_to_entry;
#endif

/* These are the holders of the old singal handlers for the trapped signals */
typedef void (*sighandler_t)(int);  /* this should be in signals.h, but is not there */


#ifdef __linux__
static __thread sighandler_t old_SIGFPE_handler;
static __thread sighandler_t old_SIGILL_handler;
static __thread sighandler_t old_SIGSEGV_handler;
#ifndef WIN32
static __thread sighandler_t old_SIGBUS_handler;
#endif
static __thread sighandler_t old_SIGABRT_handler;
static __thread sighandler_t old_SIGMEM_handler;  /* Custom signal */
static __thread sighandler_t old_SIGERR_handler;  /* Custom signal */

#else

static  sighandler_t old_SIGFPE_handler;
static  sighandler_t old_SIGILL_handler;
static  sighandler_t old_SIGSEGV_handler;
#ifndef WIN32
static  sighandler_t old_SIGBUS_handler;
#endif
static  sighandler_t old_SIGABRT_handler;
static  sighandler_t old_SIGMEM_handler;  /* Custom signal */
static  sighandler_t old_SIGERR_handler;  /* Custom signal */

#endif




/*************************************************************************
* This function prints an error message and exits 
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;

  va_start(argp, f_str);
  vfprintf(stderr, f_str, argp);
  va_end(argp);

  fprintf(stderr,"\n");
  fflush(stderr);

  exit(-2);
}


/*************************************************************************
* This function prints an error message and raises a signum signal
**************************************************************************/
void gk_errexit(int signum, char *f_str,...)
{
  va_list argp;

  va_start(argp, f_str);
  vfprintf(stderr, f_str, argp);
  va_end(argp);

  fprintf(stderr,"\n");
  fflush(stderr);

  raise(signum);
}


/***************************************************************************
* This function sets a number of signal handlers and sets the return point 
* of a longjmp
****************************************************************************/
void gk_SetSignalHandlers() 
{
  old_SIGFPE_handler  = signal(SIGFPE,  gk_NonLocalExit_Handler);
  old_SIGILL_handler  = signal(SIGILL,  gk_NonLocalExit_Handler);
  old_SIGSEGV_handler = signal(SIGSEGV, gk_NonLocalExit_Handler);
#ifndef WIN32
  old_SIGBUS_handler  = signal(SIGBUS,  gk_NonLocalExit_Handler);
#endif
  old_SIGABRT_handler = signal(SIGABRT, gk_NonLocalExit_Handler);
  if (SIGMEM != SIGABRT) 
    old_SIGMEM_handler  = signal(SIGMEM,  gk_NonLocalExit_Handler);
  if (SIGERR != SIGABRT) 
    old_SIGERR_handler  = signal(SIGERR,  gk_NonLocalExit_Handler);
}
  

/***************************************************************************
* This function sets the handlers for the signals to their default handlers
****************************************************************************/
void gk_UnsetSignalHandlers() 
{
  signal(SIGFPE,  old_SIGFPE_handler);
  signal(SIGILL,  old_SIGILL_handler);
  signal(SIGSEGV, old_SIGSEGV_handler);
#ifndef WIN32
  signal(SIGBUS,  old_SIGBUS_handler);
#endif
  signal(SIGABRT, old_SIGABRT_handler);
  if (SIGMEM != SIGABRT) 
    signal(SIGMEM,  old_SIGMEM_handler);
  if (SIGERR != SIGABRT) 
    signal(SIGERR,  old_SIGERR_handler);
}
  

/*************************************************************************
* This function is the handler for SIGUSR1 that implements the cleaning up 
* process prior to a non-local exit.
**************************************************************************/
void gk_NonLocalExit_Handler(int signum)
{
  gk_malloc_cleanup();

  //printf("Calling longjmp...\n");
  longjmp(gk_return_to_entry, signum);
}
  

/*************************************************************************/
/*! \brief Thread-safe implementation of strerror() */
/**************************************************************************/
char *gk_strerror(int errnum)
{
#ifdef __linux__
  static __thread char buf[1024];
#else
  static char buf[1024];
#endif

  strerror_r(errnum, buf, 1024);

  buf[1023] = '\0';
  return buf;
}



/*************************************************************************
* This function prints a backtrace of callinf functions
**************************************************************************/
void PrintBackTrace()
{
#ifdef LINUX
  void *array[10];
  int i, size;
  char **strings;

  size = backtrace(array, 10);
  strings = backtrace_symbols(array, size);
  
  printf("Obtained %d stack frames.\n", size);
  for (i=0; i<size; i++) {
    printf("%s\n", strings[i]);
  }
  free(strings);
#endif
}
