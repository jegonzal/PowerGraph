/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * time.c
 *
 * This file contains various functions that do with time
 *
 * $Id: timers.c 934 2007-02-21 18:31:42Z karypis $
 *
 */

#include <GKlib.h>




/*************************************************************************
* This function returns the CPU seconds
**************************************************************************/
gk_wclock_t gk_WClockSeconds(void)
{
  return time(NULL);
}


/*************************************************************************
* This function returns the CPU seconds
**************************************************************************/
double gk_CPUSeconds(void)
{
#ifdef __OPENMP__
  return omp_get_wtime();
#else
  #ifdef WIN32
    return((double) clock()/CLOCKS_PER_SEC);
  #else
    struct rusage r;

    getrusage(RUSAGE_SELF, &r);
    return ((r.ru_utime.tv_sec + r.ru_stime.tv_sec) + 1.0e-6*(r.ru_utime.tv_usec + r.ru_stime.tv_usec));
  #endif
#endif
}

