/*!
\file gk_arch.h
\brief This file contains various architecture-specific declerations

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_arch.h 1430 2007-04-07 17:53:07Z karypis $ \endverbatim
*/

#ifndef _GK_ARCH_H_
#define _GK_ARCH_H_

/*************************************************************************
* Architecture-specific differences in header files
**************************************************************************/
#ifdef LINUX
#if !defined(__USE_XOPEN)
#define __USE_XOPEN
#endif
#if !defined(_XOPEN_SOURCE)
#define _XOPEN_SOURCE 600
#endif
#if !defined(__USE_XOPEN2K)
#define __USE_XOPEN2K 
#endif
#include <execinfo.h>
#endif


#ifdef __MSC__ 
  #include <crtdefs.h>
#else
  #include <stdint.h>
  #include <inttypes.h>
  #include <sys/types.h>
  #include <sys/resource.h>
  #include <sys/time.h>
#endif



/*************************************************************************
* Architecture-specific modifications
**************************************************************************/
#ifdef WIN32
#define __thread __declspec( thread )

typedef __int32                 int32_t;
typedef __int64                 int64_t;
typedef unsigned __int32        uint32_t;
typedef unsigned __int64        uint64_t;
#endif


#ifdef SUNOS
#define __thread 
#endif


#ifdef DARWIN
#define __thread 
#endif


#ifdef __MSC__
#define rint(x) ((int)((x)+0.5))  /* MSC does not have rint() function */
#endif


#ifdef __GNUC__
#if !defined(strdup)
extern char* strdup (const char *);
#endif
#endif


#endif
