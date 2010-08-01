/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * fortran.c
 *
 * This file contains code for the fortran to C interface
 *
 * Started 8/19/97
 * George
 *
 */

#include <metislib.h>


/*************************************************************************
* This function changes the numbering to start from 0 instead of 1
**************************************************************************/
void Change2CNumbering(idxtype nvtxs, idxtype *xadj, idxtype *adjncy)
{
  idxtype i, nedges;

  for (i=0; i<=nvtxs; i++)
    xadj[i]--;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]--;
}

/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumbering(idxtype nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vector)
{
  idxtype i, nedges;

  for (i=0; i<nvtxs; i++)
    vector[i]++;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;
}

/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumbering2(idxtype nvtxs, idxtype *xadj, idxtype *adjncy)
{
  idxtype i, nedges;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;
}



/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumberingOrder(idxtype nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *v1, idxtype *v2)
{
  idxtype i, nedges;

  for (i=0; i<nvtxs; i++) {
    v1[i]++;
    v2[i]++;
  }

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;

}



/*************************************************************************
* This function changes the numbering to start from 0 instead of 1
**************************************************************************/
void ChangeMesh2CNumbering(idxtype n, idxtype *mesh)
{
  idxtype i;

  for (i=0; i<n; i++)
    mesh[i]--;

}


/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void ChangeMesh2FNumbering(idxtype n, idxtype *mesh, idxtype nvtxs, idxtype *xadj, idxtype *adjncy)
{
  idxtype i, nedges;

  for (i=0; i<n; i++)
    mesh[i]++;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;

}


/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void ChangeMesh2FNumbering2(idxtype n, idxtype *mesh, idxtype ne, idxtype nn, idxtype *epart, idxtype *npart)
{
  idxtype i, nedges;

  for (i=0; i<n; i++)
    mesh[i]++;

  for (i=0; i<ne; i++)
    epart[i]++;

  for (i=0; i<nn; i++)
    npart[i]++;

}


/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void ChangeMesh2FNumbering3(idxtype n, idxtype *mesh)
{
  idxtype i;

  for (i=0; i<n; i++)
    mesh[i]++;
}
