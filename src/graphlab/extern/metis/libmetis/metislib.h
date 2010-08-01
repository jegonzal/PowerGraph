/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metislib.h,v 1.1 2002/08/10 06:29:31 karypis Exp $
 */

#include <GKlib.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <signal.h>
#include <setjmp.h>
#include <assert.h>

#if defined(ENABLE_OPENMP)
  #include <omp.h>
#endif


#include <metis.h>

#include <defs.h>
#include <struct.h>
#include <macros.h>
#include <rename.h>
#include <proto.h>


#if defined(COMPILER_MSC)
#define rint(x) ((idxtype)((x)+0.5))  /* MSC does not have rint() function */
#endif


#if defined(COMPILER_GCC)
//extern char* strdup (const char *);
#endif

