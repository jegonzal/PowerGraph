/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ccgraph.c
 *
 * This file contains the functions that create the coarse graph
 *
 * Started 8/11/97
 * George
 *
 */

#include <metislib.h>



/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, idxtype cnvtxs, idxtype *match, idxtype *perm)
{
  idxtype i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask, dovsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;

  dovsize = (ctrl->optype == OP_KVMETIS ? 1 : 0);

  mask = HTLENGTH;
  if (cnvtxs < 8*mask || graph->nedges/graph->nvtxs > 15) { 
    CreateCoarseGraphNoMask(ctrl, graph, cnvtxs, match, perm);
    return;
  }

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs, dovsize);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cvsize = cgraph->vsize;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore; 
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1)); 

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      gk_fcopy(ncon, nvwgt+v*ncon, cnvwgt+cnvtxs*ncon);

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m] += adjwgt[j];
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj] += adjwgt[j];
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = adjwgt[j];
        }
      }
    }

    if (v != u) { 
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        gk_faxpy(ncon, 1.0, nvwgt+u*ncon, 1, cnvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m] += adjwgt[j];
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj] += adjwgt[j];
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = adjwgt[j];
          }
        }
      }

      /* Remove the contracted adjacency weight */
      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs) 
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) { /* This 2nd check is needed for non-adjacent matchings */
        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraphNoMask(CtrlType *ctrl, GraphType *graph, idxtype cnvtxs, idxtype *match, idxtype *perm)
{
  idxtype i, j, k, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, dovsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;

  dovsize = (ctrl->optype == OP_KVMETIS ? 1 : 0);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;


  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs, dovsize);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cvsize = cgraph->vsize;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  htable = idxset(cnvtxs, -1, idxwspacemalloc(ctrl, cnvtxs));

  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore; 
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      gk_fcopy(ncon, nvwgt+v*ncon, cnvwgt+cnvtxs*ncon);

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      if ((m = htable[k]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[k] = nedges++;
      }
      else {
        cadjwgt[m] += adjwgt[j];
      }
    }

    if (v != u) { 
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        gk_faxpy(ncon, 1.0, nvwgt+u*ncon, 1, cnvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        if ((m = htable[k]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[k] = nedges++;
        }
        else {
          cadjwgt[m] += adjwgt[j];
        }
      }

      /* Remove the contracted adjacency weight */
      if ((j = htable[cnvtxs]) != -1) {
        ASSERT(cadjncy[j] == cnvtxs);
        cadjwgtsum[cnvtxs] -= cadjwgt[j];
        cadjncy[j] = cadjncy[--nedges];
        cadjwgt[j] = cadjwgt[nedges];
        htable[cnvtxs] = -1;
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d\n", cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt)));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]] = -1;  /* Zero out the htable */

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, cnvtxs);
}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph_NVW(CtrlType *ctrl, GraphType *graph, idxtype cnvtxs, idxtype *match, idxtype *perm)
{
  idxtype i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask;
  idxtype *xadj, *adjncy, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;


  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs, 0);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore; 
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  mask = HTLENGTH;
  htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1)); 

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    cvwgt[cnvtxs] = 1;
    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = 1;
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m]++;
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj]++;
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = 1;
        }
      }
    }

    if (v != u) { 
      cvwgt[cnvtxs]++;
      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = 1;
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m]++;
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj]++;
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = 1;
          }
        }
      }

      /* Remove the contracted adjacency weight */
      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs) 
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) { /* This 2nd check is needed for non-adjacent matchings */
        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, 0);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}


/*************************************************************************
* Setup the various arrays for the coarse graph
**************************************************************************/
GraphType *SetUpCoarseGraph(GraphType *graph, idxtype cnvtxs, idxtype dovsize)
{
  GraphType *cgraph;

  cgraph = CreateGraph();

  cgraph->nvtxs = cnvtxs;
  cgraph->ncon  = graph->ncon;

  cgraph->finer  = graph;
  graph->coarser = cgraph;


  /* Allocate memory for the coarser graph */
  cgraph->xadj       = idxmalloc(cnvtxs+1, "SetUpCoarseGraph: xadj");
  cgraph->adjwgtsum  = idxmalloc(cnvtxs,   "SetUpCoarseGraph: adjwgtsum");
  cgraph->cmap       = idxmalloc(cnvtxs,   "SetUpCoarseGraph: cmap");
  cgraph->adjncy     = idxmalloc(graph->nedges,   "SetUpCoarseGraph: adjncy");
  cgraph->adjwgt     = idxmalloc(graph->nedges,   "SetUpCoarseGraph: adjwgt");

  if (graph->ncon == 1)
    cgraph->vwgt     = idxmalloc(cnvtxs,   "SetUpCoarseGraph: vwgt");
  else
    cgraph->nvwgt    = gk_fmalloc(graph->ncon*cnvtxs, "SetUpCoarseGraph: nvwgt");

  if (dovsize)
    cgraph->vsize    = idxmalloc(cnvtxs,   "SetUpCoarseGraph: vsize");

  return cgraph;
}



/*************************************************************************
* This function re-adjusts the amount of memory that was allocated if
* it will lead to significant savings
**************************************************************************/
void ReAdjustMemory(GraphType *graph, GraphType *cgraph, idxtype dovsize) 
{
  if (cgraph->nedges > 10000 && cgraph->nedges < 0.8*graph->nedges) {
    cgraph->adjncy = idxrealloc(cgraph->adjncy, cgraph->nedges, "ReAdjustMemory: adjncy");
    cgraph->adjwgt = idxrealloc(cgraph->adjwgt, cgraph->nedges, "ReAdjustMemory: adjwgt");
  }
}
