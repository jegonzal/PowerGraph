/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 2/24/96
 * George
 *
 * $Id: memory.c,v 1.2 2002/08/10 06:29:31 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void AllocateWorkSpace(CtrlType *ctrl, GraphType *graph, idxtype nparts)
{
  ctrl->wspace.pmat = NULL;

  if (ctrl->optype == OP_KMETIS) {
    ctrl->wspace.edegrees = (EDegreeType *)gk_malloc(graph->nedges*sizeof(EDegreeType), "AllocateWorkSpace: edegrees");
    ctrl->wspace.vedegrees = NULL;
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.edegrees;

    ctrl->wspace.pmat = idxmalloc(nparts*nparts, "AllocateWorkSpace: pmat");

    /* Memory requirements for different phases
          Coarsening
                    Matching: 4*nvtxs vectors
                 Contraction: 2*nvtxs vectors (from the above 4), 1*nparts, 1*Nedges
            Total = MAX(4*nvtxs, 2*nvtxs+nparts+nedges)

          Refinement
                Random Refinement/Balance: 5*nparts + 1*nvtxs + 2*nedges
                Greedy Refinement/Balance: 5*nparts + 2*nvtxs + 2*nedges + 1*PQueue(==Nvtxs)
            Total = 5*nparts + 3*nvtxs + 2*nedges

         Total = 5*nparts + 3*nvtxs + 2*nedges 
    */
    ctrl->wspace.maxcore = 3*(graph->nvtxs+1) +                 /* Match/Refinement vectors */
                           5*(nparts+1) +                       /* Partition weights etc */
                           graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* Greedy k-way balance/refine */
                           20  /* padding for 64 bit machines */
                           ;
  }
  else if (ctrl->optype == OP_KVMETIS) {
    ctrl->wspace.edegrees = NULL;
    ctrl->wspace.vedegrees = (VEDegreeType *)gk_malloc(graph->nedges*sizeof(VEDegreeType), "AllocateWorkSpace: vedegrees");
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.vedegrees;

    ctrl->wspace.pmat = idxmalloc(nparts*nparts, "AllocateWorkSpace: pmat");

    /* Memory requirements for different phases are identical to KMETIS */
    ctrl->wspace.maxcore = 3*(graph->nvtxs+1) +                 /* Match/Refinement vectors */
                           3*(nparts+1) +                       /* Partition weights etc */
                           graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* Greedy k-way balance/refine */
                           20  /* padding for 64 bit machines */
                           ;
  }
  else {
    ctrl->wspace.edegrees = (EDegreeType *)idxmalloc(graph->nedges, "AllocateWorkSpace: edegrees");
    ctrl->wspace.vedegrees = NULL;
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.edegrees;

    ctrl->wspace.maxcore = 5*(graph->nvtxs+1) +                 /* Refinement vectors */
                           4*(nparts+1) +                       /* Partition weights etc */
                           2*graph->ncon*graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* 2-way refinement */
                           2*graph->ncon*(NEG_GAINSPAN+PLUS_GAINSPAN+1)*(sizeof(ListNodeType *)/sizeof(idxtype)) + /* 2-way refinement */
                           20  /* padding for 64 bit machines */
                           ;
  }

  ctrl->wspace.maxcore += HTLENGTH;
  ctrl->wspace.core = idxmalloc(ctrl->wspace.maxcore, "AllocateWorkSpace: maxcore");
  ctrl->wspace.ccore = 0;
}


/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void FreeWorkSpace(CtrlType *ctrl, GraphType *graph)
{
  gk_free((void **)&ctrl->wspace.edegrees, &ctrl->wspace.vedegrees, &ctrl->wspace.core, &ctrl->wspace.pmat, LTERM);
}

/*************************************************************************
* This function returns how may words are left in the workspace
**************************************************************************/
idxtype WspaceAvail(CtrlType *ctrl)
{
  return ctrl->wspace.maxcore - ctrl->wspace.ccore;
}


/*************************************************************************
* This function allocate space from the core 
**************************************************************************/
idxtype *idxwspacemalloc(CtrlType *ctrl, idxtype n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore += n;
  ASSERT(ctrl->wspace.ccore <= ctrl->wspace.maxcore);
  return ctrl->wspace.core + ctrl->wspace.ccore - n;
}

/*************************************************************************
* This function frees space from the core 
**************************************************************************/
void idxwspacefree(CtrlType *ctrl, idxtype n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore -= n;
  ASSERT(ctrl->wspace.ccore >= 0);
}


/*************************************************************************
* This function allocate space from the core 
**************************************************************************/
float *fwspacemalloc(CtrlType *ctrl, idxtype n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore += n;
  ASSERT(ctrl->wspace.ccore <= ctrl->wspace.maxcore);
  return (float *) (ctrl->wspace.core + ctrl->wspace.ccore - n);
}

/*************************************************************************
* This function frees space from the core 
**************************************************************************/
void fwspacefree(CtrlType *ctrl, idxtype n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore -= n;
  ASSERT(ctrl->wspace.ccore >= 0);
}



/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
GraphType *CreateGraph(void)
{
  GraphType *graph;

  graph = (GraphType *)gk_malloc(sizeof(GraphType), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
void InitGraph(GraphType *graph) 
{
  /* graph size constants */
  graph->nvtxs     = -1;
  graph->nedges    = -1;
  graph->ncon      = -1;
  graph->mincut    = -1;
  graph->minvol    = -1;
  graph->nbnd      = -1;

  /* memory for the graph structure */
  graph->xadj      = NULL;
  graph->vwgt      = NULL;
  graph->nvwgt     = NULL;
  graph->vsize     = NULL;
  graph->adjncy    = NULL;
  graph->adjwgt    = NULL;
  graph->adjwgtsum = NULL;
  graph->label     = NULL;
  graph->cmap      = NULL;
  graph->coords    = NULL;

  /* by default these are set to true, but the can be explicitly changed afterwards */
  graph->free_xadj   = 1;
  graph->free_vwgt   = 1;
  graph->free_vsize  = 1;
  graph->free_adjncy = 1;
  graph->free_adjwgt = 1;


  /* memory for the partition/refinement structure */
  graph->where     = NULL;
  graph->pwgts     = NULL;
  graph->npwgts    = NULL;
  graph->id        = NULL;
  graph->ed        = NULL;
  graph->bndptr    = NULL;
  graph->bndind    = NULL;
  graph->rinfo     = NULL;
  graph->vrinfo    = NULL;
  graph->nrinfo    = NULL;

  /* linked-list structure */
  graph->coarser   = NULL;
  graph->finer     = NULL;
}



/*************************************************************************
* This function allocates memory for the requested fields of the graph's
* structure (i.e., not partition/refinement structure). 
* The size of these fields is determined by the size of the graph itself
* The fields should be provided as a comma-separated string with no spaces
* For example 
*   AllocGraphFields(&graph, "xadj,adjncy,vwgt")
**************************************************************************/
void AllocGraphFields(GraphType *graph, char *fields)
{
  char *sptr, *eptr;

}


/*************************************************************************
* This function frees the refinement/partition memory
* any memory stored in a graph
**************************************************************************/
void FreeRData(GraphType *graph) 
{
  /* free partition/refinement structure */
  gk_free((void **)&graph->where, &graph->pwgts, &graph->npwgts, &graph->id,
    &graph->ed, &graph->bndptr, &graph->bndind, &graph->rinfo, &graph->vrinfo,
    &graph->nrinfo, LTERM);
}


/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeGraph(GraphType *graph, int flag) 
{
  
  /* free graph structure */
  if (graph->free_xadj)
    gk_free((void **)&graph->xadj, LTERM);
  if (graph->free_vwgt)
    gk_free((void **)&graph->vwgt, LTERM);
  if (graph->free_vsize)
    gk_free((void **)&graph->vsize, LTERM);
  if (graph->free_adjncy)
    gk_free((void **)&graph->adjncy, LTERM);
  if (graph->free_adjwgt)
    gk_free((void **)&graph->adjwgt, LTERM);
    
  gk_free((void **)&graph->nvwgt, &graph->adjwgtsum, &graph->label, &graph->cmap, 
    &graph->coords, LTERM);

  /* free partition/refinement structure */
  gk_free((void **)&graph->where, &graph->pwgts, &graph->npwgts, &graph->id,
    &graph->ed, &graph->bndptr, &graph->bndind, &graph->rinfo, &graph->vrinfo,
    &graph->nrinfo, LTERM);

  if (flag)
    gk_free((void **)&graph, LTERM);
}

