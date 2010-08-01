/*
 * GKlib.h
 * 
 * George's library of most frequently used routines
 *
 * $Id: GKlib.h 1430 2007-04-07 17:53:07Z karypis $
 *
 */

#ifndef _GKLIB_H_
#define _GKLIB_H_ 1

#define GKMSPACE

#if defined(_MSC_VER)
#define __MSC__
#endif
#if defined(__ICC)
#define __ICC__
#endif
#include <signal.h>

#include <gk_arch.h> /*!< This should be here, prior to the includes */


/*************************************************************************
* Header file inclusion section
**************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <limits.h>

#include <setjmp.h>
#include <assert.h>
#include <sys/stat.h>

#if defined(__WITHPCRE__)
#include <pcreposix.h>
#else
#include <regex.h>
#endif


#if defined(__OPENMP__) 
#include <omp.h>
#endif




#include <gk_dlmalloc.h>
#include <gk_types.h>
#include <gk_struct.h>
#include <gk_externs.h>
#include <gk_defs.h>
#include <gk_macros.h>
#include <gk_getopt.h>

#include <gk_sort.h>
#include <gk_blas.h>
#include <gk_memory.h>

#include <gk_proto.h>

#endif  /* GKlib.h */


