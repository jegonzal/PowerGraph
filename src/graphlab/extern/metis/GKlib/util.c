/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This file contains utility functions
 *
 * $Id: util.c 1421 2007-04-06 14:37:41Z karypis $
 *
 */

#include <GKlib.h>



/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void gk_RandomPermute(size_t n, int *p, int flag)
{
  gk_idx_t i, u, v;
  int tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  for (i=0; i<n/2; i++) {
    v = RandomInRange(n);
    u = RandomInRange(n);
    SWAP(p[v], p[u], tmp);
  }
}


/*************************************************************************
* This function converts an element-based set membership to a set-based
* element membership.
* That is, it takes an array such as part[] that stores where each element
* belongs to, and returns a pair of arrays (pptr[], pind[]) that store in
* CSF format the list of elements belonging in each partition.
* Parameters:
*    n      - the number of elements in the array (e.g., # of vertices)
*    range  - the cardinality of the set (e.g., # of partitions)
**************************************************************************/
void gk_array2csr(size_t n, size_t range, int *array, int *ptr, int *ind)
{
  gk_idx_t i;

  gk_iset(range+1, 0, ptr);

  for (i=0; i<n; i++) 
    ptr[array[i]]++;

  /* Compute the ptr, ind structure */
  MAKECSR(i, range, ptr);
  for (i=0; i<n; i++)
    ind[ptr[array[i]]++] = i;
  SHIFTCSR(i, range, ptr);
}



/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
int gk_log2(int a)
{
  gk_idx_t i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}


/*************************************************************************
* This function checks if the argument is a power of 2
**************************************************************************/
int gk_ispow2(int a)
{
  return (a == (1<<gk_log2(a)));
}


/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
float gk_flog2(float a)
{
  return log(a)/log(2.0);
}

