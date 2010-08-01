/*
 * mutil.c 
 *
 * This file contains various utility functions for the MOC portion of the
 * code
 *
 * Started 2/15/98
 * George
 *
 * $Id: mutil.c,v 1.2 2002/08/10 06:29:33 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function checks if the vertex weights of two vertices are below 
* a given set of values
**************************************************************************/
idxtype AreAllVwgtsBelow(idxtype ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float limit)
{
  idxtype i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] > limit)
      return 0;

  return 1;
}


/*************************************************************************
* This function checks if the vertex weights of two vertices are below 
* a given set of values
**************************************************************************/
idxtype AreAnyVwgtsBelow(idxtype ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float limit)
{
  idxtype i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] < limit)
      return 1;

  return 0;
}



/*************************************************************************
* This function checks if the vertex weights of two vertices are above 
* a given set of values
**************************************************************************/
idxtype AreAllVwgtsAbove(idxtype ncon, float alpha, float *vwgt1, float beta, float *vwgt2, float limit)
{
  idxtype i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] < limit)
      return 0;

  return 1;
}


/*************************************************************************
* This function computes the load imbalance over all the constrains
* For now agk_fsume that we just want balanced partitionings
**************************************************************************/ 
float ComputeLoadImbalance(idxtype ncon, idxtype nparts, float *npwgts, float *tpwgts)
{
  idxtype i, j;
  float max, lb=0.0;

  for (i=0; i<ncon; i++) {
    max = 0.0;
    for (j=0; j<nparts; j++) {
      if (npwgts[j*ncon+i] > max)
        max = npwgts[j*ncon+i];
    }
    if (max*nparts > lb)
      lb = max*nparts;
  }

  return lb;
}

/*************************************************************************
* This function checks if the vertex weights of two vertices are below 
* a given set of values
**************************************************************************/
idxtype AreAllBelow(idxtype ncon, float *v1, float *v2)
{
  idxtype i;

  for (i=0; i<ncon; i++)
    if (v1[i] > v2[i])
      return 0;

  return 1;
}
