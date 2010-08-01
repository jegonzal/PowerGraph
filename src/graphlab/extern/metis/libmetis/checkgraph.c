/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * checkgraph.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 */

#include <metislib.h>



/*************************************************************************
* This function checks if a graph is valid
**************************************************************************/
idxtype CheckGraph(GraphType *graph)
{
  idxtype i, j, k, l;
  idxtype nvtxs, ncon, err=0;
  idxtype minedge, maxedge, minewgt, maxewgt;
  float minvwgt[MAXNCON], maxvwgt[MAXNCON];
  idxtype *xadj, *adjncy, *adjwgt, *htable;
  float *nvwgt, ntvwgts[MAXNCON];

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  htable = idxsmalloc(nvtxs, 0, "htable");

  if (ncon > 1) {
    for (j=0; j<ncon; j++) { 
      minvwgt[j] = maxvwgt[j] = nvwgt[j];
      ntvwgts[j] = 0.0;
    }
  }

  minedge = maxedge = adjncy[0];
  minewgt = maxewgt = adjwgt[0];

  for (i=0; i<nvtxs; i++) {
    if (ncon > 1) {
      for (j=0; j<ncon; j++) {
        ntvwgts[j] += nvwgt[i*ncon+j];
        minvwgt[j] = (nvwgt[i*ncon+j] < minvwgt[j]) ? nvwgt[i*ncon+j] : minvwgt[j];
        maxvwgt[j] = (nvwgt[i*ncon+j] > maxvwgt[j]) ? nvwgt[i*ncon+j] : maxvwgt[j];
      }
    }

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];

      minedge = (k < minedge) ? k : minedge;
      maxedge = (k > maxedge) ? k : maxedge;
      minewgt = (adjwgt[j] < minewgt) ? adjwgt[j] : minewgt;
      maxewgt = (adjwgt[j] > maxewgt) ? adjwgt[j] : maxewgt;

      if (i == k) {
        mprintf("Vertex %D contains a self-loop (i.e., diagonal entry in the matrix)!\n", i);
        err++;
      }
      else {
        for (l=xadj[k]; l<xadj[k+1]; l++) {
          if (adjncy[l] == i) {
            if (adjwgt != NULL && adjwgt[l] != adjwgt[j]) {
              mprintf("Edges (%D %D) and (%D %D) do not have the same weight! %D %D\n", i,k,k,i, adjwgt[l], adjwgt[j]);
              err++;
            }
            break;
          }
        }
        if (l == xadj[k+1]) {
          mprintf("Missing edge: (%D %D)!\n", k, i);
          err++;
        }
      }

      if (htable[k] == 0) {
        htable[k]++;
      }
      else {
        mprintf("Edge %D from vertex %D is repeated %D times\n", k, i, htable[k]++);
        err++;
      }
    }

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      htable[adjncy[j]] = 0;
    }
  }

  if (ncon > 1) {
    for (j=0; j<ncon; j++) {
      if (fabs(ntvwgts[j] - 1.0) > 0.0001) {
        mprintf("Normalized vwgts don't sum to one.  Weight %D = %.8f.\n", j, ntvwgts[j]);
        err++;
      }
    }
  }

/*
  mprintf("errs: %D, adjncy: [%D %D], adjwgt: [%D %D]\n",
  err, minedge, maxedge, minewgt, maxewgt);
  if (ncon > 1) {
    for (j=0; j<ncon; j++)
      mprintf("[%.5f %.5f] ", minvwgt[j], maxvwgt[j]);
    mprintf("\n");
  }
*/
 
  if (err > 0) { 
    mprintf("A total of %D errors exist in the input file. Correct them, and run again!\n", err);
  }

  gk_free((void **)&htable, LTERM);
  return (err == 0 ? 1 : 0);
}

