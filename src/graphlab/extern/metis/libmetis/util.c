/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id: util.c,v 1.4 2003/04/13 04:45:12 karypis Exp $
 */

#include <metislib.h>




/*************************************************************************
* The following are utility functions defined for idxtype using GKlib's
* function generating macros
**************************************************************************/
GK_XMALLOC(idxmalloc, idxtype)
GK_XREALLOC(idxrealloc, idxtype)
GK_XSMALLOC(idxsmalloc, idxtype, idxset)
GK_SET(idxset, idxtype)
GK_ARGMAX(idxargmax, idxtype)
GK_ARGMIN(idxargmin, idxtype)
GK_SUM(idxsum, idxtype, idxtype)
GK_AXPY(idxaxpy, idxtype)




/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
idxtype idxargmax_strd(size_t n, idxtype *x, idxtype incx)
{
  size_t i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}




/*************************************************************************
* These functions return the index of the almost maximum element in a vector
**************************************************************************/
idxtype famax2(size_t n, float *x)
{
  size_t i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}







/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void RandomPermute(size_t n, idxtype *p, idxtype flag)
{
  size_t i, u, v;
  idxtype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  if (n <= 4)
    return;

  for (i=0; i<n; i+=16) {
    u = RandomInRange(n-4);
    v = RandomInRange(n-4);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  }
}


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
void InitRandom(idxtype seed)
{
  srand((seed == -1 ? 4321 : seed)); 
}


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
idxtype strtoidx(const char *nptr, char **endptr, int base)
{
  if (sizeof(idxtype) == sizeof(long int))
    return (idxtype)strtol(nptr, endptr, base);

  if (sizeof(idxtype) == sizeof(long long int))
    return (idxtype)strtoll(nptr, endptr, base);

  return (idxtype)strtoimax(nptr, endptr, base);
}
