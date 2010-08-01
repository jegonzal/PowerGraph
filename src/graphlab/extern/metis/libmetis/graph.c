/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * graph.c
 *
 * This file contains functions that deal with setting up the graphs
 * for METIS.
 *
 * Started 7/25/97
 * George
 *
 */

#include <metislib.h>

/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void SetUpGraph(GraphType *graph, idxtype OpType, idxtype nvtxs, idxtype ncon,
       idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype wgtflag)
{
  idxtype i, j, k, sum;
  float *nvwgt;
  idxtype tvwgt[MAXNCON];


  InitGraph(graph);

  graph->nvtxs  = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon   = ncon;

  graph->xadj        = xadj;
  graph->free_xadj   = 0;

  graph->adjncy      = adjncy;
  graph->free_adjncy = 0;

  /* setup the vertex weights */
  if (ncon == 1) { /* We are in the non mC mode */
    if ((wgtflag&2) == 0) {
      vwgt = graph->vwgt = idxsmalloc(nvtxs, 1, "SetUpGraph: vwgt");
    }
    else {
      graph->vwgt      = vwgt;
      graph->free_vwgt = 0;
    }
  }
  else {  /* Set up the graph in MOC mode */
    for (i=0; i<ncon; i++) 
      tvwgt[i] = idxsum(nvtxs, vwgt+i, ncon);
    
    nvwgt = graph->nvwgt = gk_fmalloc(ncon*nvtxs, "SetUpGraph: nvwgt");

    for (i=0; i<nvtxs; i++) {
      for (j=0; j<ncon; j++) 
        nvwgt[i*ncon+j] = (1.0*vwgt[i*ncon+j])/(1.0*tvwgt[j]);
    }
  }

  /* setup the edge weights */
  if ((wgtflag&1) == 0) {
    adjwgt = graph->adjwgt = idxsmalloc(graph->nedges, 1, "SetUpGraph: adjwgt");
  }
  else {
    graph->adjwgt      = adjwgt;
    graph->free_adjwgt = 0;
  }

  /* Compute the initial values of the adjwgtsum */
  graph->adjwgtsum = idxmalloc(nvtxs, "SetUpGraph: adjwgtsum");

  for (i=0; i<nvtxs; i++) {
    for (sum=0, j=xadj[i]; j<xadj[i+1]; j++)
      sum += adjwgt[j];
    graph->adjwgtsum[i] = sum;
  }

  graph->cmap = idxmalloc(nvtxs, "SetUpGraph: cmap");


  if (OpType != OP_KMETIS && OpType != OP_KVMETIS) {
    graph->label = idxmalloc(nvtxs, "SetUpGraph: label");

    for (i=0; i<nvtxs; i++)
      graph->label[i] = i;
  }

}




/*************************************************************************
* This function sets up the graph from the user input. The difference 
* from the previous routine is that the vertex weights came already in
* a normalized fashion.
**************************************************************************/
void SetUpGraph2(GraphType *graph, idxtype nvtxs, idxtype ncon, idxtype *xadj, 
       idxtype *adjncy, float *nvwgt, idxtype *adjwgt)
{
  idxtype i, j, sum;

  InitGraph(graph);

  graph->nvtxs  = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon   = ncon;

  graph->xadj        = xadj;
  graph->free_xadj   = 0;

  graph->adjncy      = adjncy;
  graph->free_adjncy = 0;

  graph->adjwgt      = adjwgt;
  graph->free_adjwgt = 0;

  graph->nvwgt = gk_fmalloc(nvtxs*ncon, "SetUpGraph2: graph->nvwgt");
  gk_fcopy(nvtxs*ncon, nvwgt, graph->nvwgt);


  /* Compute the initial values of the adjwgtsum */
  graph->adjwgtsum = idxmalloc(nvtxs, "SetUpGraph2: adjwgtsum");
  for (i=0; i<nvtxs; i++) {
    for (sum=0, j=xadj[i]; j<xadj[i+1]; j++)
      sum += adjwgt[j];
    graph->adjwgtsum[i] = sum;
  }

  graph->cmap = idxmalloc(nvtxs, "SetUpGraph2: cmap");

  graph->label = idxmalloc(nvtxs, "SetUpGraph: label");
  for (i=0; i<nvtxs; i++)
    graph->label[i] = i;

}


/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void VolSetUpGraph(GraphType *graph, idxtype OpType, idxtype nvtxs, idxtype ncon, idxtype *xadj, 
                   idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype wgtflag)
{
  idxtype i, j, k, sum;
  idxtype *adjwgt;
  float *nvwgt;
  idxtype tvwgt[MAXNCON];

  InitGraph(graph);

  graph->nvtxs  = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->ncon   = ncon;

  graph->xadj      = xadj;
  graph->free_xadj = 0;

  graph->adjncy      = adjncy;
  graph->free_adjncy = 0;

  /* Setup the vwgt/vwgt */
  if (ncon == 1) { /* We are in the non mC mode */
    if ((wgtflag&2) == 0) {
      vwgt = graph->vwgt = idxsmalloc(nvtxs, 1, "VolSetUpGraph: vwgt");
    }
    else {
      graph->vwgt      = vwgt;
      graph->free_vwgt = 0;
    }
  }
  else {  /* Set up the graph in MOC mode */
    /* Create the normalized vertex weights along each constraint */
    for (i=0; i<ncon; i++) 
      tvwgt[i] = idxsum(nvtxs, vwgt+i, ncon);
    
    nvwgt = graph->nvwgt = gk_fmalloc(ncon*nvtxs, "SetUpGraph: nvwgt");

    for (i=0; i<nvtxs; i++) {
      for (j=0; j<ncon; j++) 
        nvwgt[i*ncon+j] = (1.0*vwgt[i*ncon+j])/(1.0*tvwgt[j]);
    }
  }


  /* Setup the vsize */
  if ((wgtflag&1) == 0) {
    vsize = graph->vsize = idxsmalloc(nvtxs, 1, "VolSetUpGraph: vsize");
  }
  else {
    graph->vsize      = vsize;
    graph->free_vsize = 0;
  }

  /* Allocate memory for edge weights and initialize them to the sum of the vsize */
  adjwgt = graph->adjwgt = idxmalloc(graph->nedges, "VolSetUpGraph: adjwgt");
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++)
      adjwgt[j] = 1+vsize[i]+vsize[adjncy[j]];
  }


  /* Compute the initial values of the adjwgtsum */
  graph->adjwgtsum = idxmalloc(nvtxs, "VolSetUpGraph: adjwgtsum");
  for (i=0; i<nvtxs; i++) {
    for (sum=0, j=xadj[i]; j<xadj[i+1]; j++)
      sum += adjwgt[j];
    graph->adjwgtsum[i] = sum;
  }

  graph->cmap = idxmalloc(nvtxs, "VolSetUpGraph: cmap");


  if (OpType != OP_KVMETIS) {
    graph->label = idxmalloc(nvtxs, "SetUpGraph: label");

    for (i=0; i<nvtxs; i++)
      graph->label[i] = i;
  }

}


/*************************************************************************
* This function randomly permutes the adjacency lists of a graph
**************************************************************************/
void RandomizeGraph(GraphType *graph)
{
  idxtype i, j, k, l, tmp, nvtxs;
  idxtype *xadj, *adjncy, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  for (i=0; i<nvtxs; i++) {
    l = xadj[i+1]-xadj[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = xadj[i] + RandomInRange(l);
      SWAP(adjncy[j], adjncy[k], tmp);
      SWAP(adjwgt[j], adjwgt[k], tmp);
    }
  }
}


/*************************************************************************
* This function checks whether or not partition pid is contigous
**************************************************************************/
idxtype IsConnectedSubdomain(CtrlType *ctrl, GraphType *graph, idxtype pid, idxtype report)
{
  idxtype i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  idxtype *xadj, *adjncy, *where, *touched, *queue;
  idxtype *cptr;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue   = idxmalloc(nvtxs, "IsConnected: queue");
  cptr    = idxmalloc(nvtxs+1, "IsConnected: cptr");

  nleft = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] == pid) 
      nleft++;
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] == pid) 
      break;
  }

  touched[i] = 1;
  queue[0] = i;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (where[i] == pid && !touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] == pid && !touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  if (ncmps > 1 && report) {
    mprintf("The graph has %D connected components in partition %D:\t", ncmps, pid);
    for (i=0; i<ncmps; i++) {
      wgt = 0;
      for (j=cptr[i]; j<cptr[i+1]; j++)
        wgt += graph->vwgt[queue[j]];
      mprintf("[%5D %5D] ", cptr[i+1]-cptr[i], wgt);
      /*
      if (cptr[i+1]-cptr[i] == 1)
        mprintf("[%D %D] ", queue[cptr[i]], xadj[queue[cptr[i]]+1]-xadj[queue[cptr[i]]]);
      */
    }
    mprintf("\n");
  }

  gk_free((void **)&touched, &queue, &cptr, LTERM);

  return (ncmps == 1 ? 1 : 0);
}


/*************************************************************************
* This function checks whether a graph is contigous or not
**************************************************************************/
idxtype IsConnected(CtrlType *ctrl, GraphType *graph, idxtype report)
{
  idxtype i, j, k, nvtxs, first, last;
  idxtype *xadj, *adjncy, *touched, *queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue = idxmalloc(nvtxs, "IsConnected: queue");

  touched[0] = 1;
  queue[0] = 0;
  first = 0; last = 1;

  while (first < last) {
    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }

  if (first != nvtxs && report)
    mprintf("The graph is not connected. It has %D disconnected vertices!\n", nvtxs-first);

  return (first == nvtxs ? 1 : 0);
}


/*************************************************************************
* This function checks whether or not partition pid is contigous
**************************************************************************/
idxtype IsConnected2(GraphType *graph, idxtype report)
{
  idxtype i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  idxtype *xadj, *adjncy, *where, *touched, *queue;
  idxtype *cptr;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue   = idxmalloc(nvtxs, "IsConnected: queue");
  cptr    = idxmalloc(nvtxs+1, "IsConnected: cptr");

  nleft = nvtxs;
  touched[0] = 1;
  queue[0] = 0;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (!touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  if (ncmps > 1 && report) {
    mprintf("%D connected components:\t", ncmps);
    for (i=0; i<ncmps; i++) {
      if (cptr[i+1]-cptr[i] > 200)
        mprintf("[%5D] ", cptr[i+1]-cptr[i]);
    }
    mprintf("\n");
  }

  gk_free((void **)&touched, &queue, &cptr, LTERM);

  return (ncmps == 1 ? 1 : 0);
}


/*************************************************************************
* This function returns the number of connected components in cptr,cind
* The separator of the graph is used to split it and then find its components.
**************************************************************************/
idxtype FindComponents(CtrlType *ctrl, GraphType *graph, idxtype *cptr, idxtype *cind)
{
  idxtype i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  idxtype *xadj, *adjncy, *where, *touched, *queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: queue");

  for (i=0; i<graph->nbnd; i++)
    touched[graph->bndind[i]] = 1;

  queue = cind;

  nleft = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2) 
      nleft++;
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2)
      break;
  }

  touched[i] = 1;
  queue[0] = i;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (!touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  gk_free((void **)&touched, LTERM);

  return ncmps;
}



