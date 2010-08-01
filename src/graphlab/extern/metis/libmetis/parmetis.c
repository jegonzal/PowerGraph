/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * parmetis.c
 *
 * This file contains top level routines that are used by ParMETIS
 *
 * Started 10/14/97
 * George
 *
 * $Id: parmetis.c,v 1.5 2002/08/13 07:16:43 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function is the entry point for KMETIS with seed specification
* in options[7] 
**************************************************************************/
void METIS_PartGraphKway2(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                         idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, 
                         idxtype *options, idxtype *edgecut, idxtype *part)
{
  idxtype i;
  float *tpwgts;

  tpwgts = gk_fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++) 
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphKway2(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, 
                       tpwgts, options, edgecut, part);

  gk_free((void **)&tpwgts, LTERM);
}


/*************************************************************************
* This function is the entry point for KWMETIS with seed specification
* in options[7] 
**************************************************************************/
void METIS_WPartGraphKway2(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                          idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, 
                          float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  idxtype i, j;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_KMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = KMETIS_CTYPE;
    ctrl.IType = KMETIS_ITYPE;
    ctrl.RType = KMETIS_RTYPE;
    ctrl.dbglvl = KMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = 20*(*nparts);
  ctrl.maxvwgt = 1.5*((graph.vwgt ? idxsum(*nvtxs, graph.vwgt, 1) : (*nvtxs))/ctrl.CoarsenTo);

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

  *edgecut = MlevelKWayPartitioning(&ctrl, &graph, *nparts, part, tpwgts, 1.03);

  IFSET(ctrl.dbglvl, DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function is the entry point for the node ND code for ParMETIS
**************************************************************************/
void METIS_NodeNDP(idxtype nvtxs, idxtype *xadj, idxtype *adjncy, idxtype npes, 
                   idxtype *options, idxtype *perm, idxtype *iperm, idxtype *sizes) 
{
  idxtype i, ii, j, l, wflag, nflag;
  GraphType graph;
  CtrlType ctrl;
  idxtype *cptr, *cind;

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType   = ONMETIS_CTYPE;
    ctrl.IType   = ONMETIS_ITYPE;
    ctrl.RType   = ONMETIS_RTYPE;
    ctrl.dbglvl  = ONMETIS_DBGLVL;
    ctrl.oflags  = ONMETIS_OFLAGS;
    ctrl.pfactor = ONMETIS_PFACTOR;
    ctrl.nseps   = ONMETIS_NSEPS;
  }
  else {
    ctrl.CType   = options[OPTION_CTYPE];
    ctrl.IType   = options[OPTION_ITYPE];
    ctrl.RType   = options[OPTION_RTYPE];
    ctrl.dbglvl  = options[OPTION_DBGLVL];
    ctrl.oflags  = options[OPTION_OFLAGS];
    ctrl.pfactor = options[OPTION_PFACTOR];
    ctrl.nseps   = options[OPTION_NSEPS];
  }
  if (ctrl.nseps < 1)
    ctrl.nseps = 1;

  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

  InitRandom(-1);

  if (ctrl.oflags&OFLAG_COMPRESS) {
    /*============================================================
    * Compress the graph 
    ==============================================================*/
    cptr = idxmalloc(nvtxs+1, "ONMETIS: cptr");
    cind = idxmalloc(nvtxs, "ONMETIS: cind");

    CompressGraph(&ctrl, &graph, nvtxs, xadj, adjncy, cptr, cind);

    if (graph.nvtxs >= COMPRESSION_FRACTION*(nvtxs)) {
      ctrl.oflags--; /* We actually performed no compression */
      gk_free((void **)&cptr, &cind, LTERM);
    }
    else if (2*graph.nvtxs < nvtxs && ctrl.nseps == 1)
      ctrl.nseps = 2;
  }
  else {
    SetUpGraph(&graph, OP_ONMETIS, nvtxs, 1, xadj, adjncy, NULL, NULL, 0);
  }


  /*=============================================================
  * Do the nested dissection ordering 
  --=============================================================*/
  ctrl.maxvwgt = 1.5*(idxsum(graph.nvtxs, graph.vwgt, 1)/ctrl.CoarsenTo);
  AllocateWorkSpace(&ctrl, &graph, 2);

  idxset(2*npes-1, 0, sizes);
  MlevelNestedDissectionP(&ctrl, &graph, iperm, graph.nvtxs, npes, 0, sizes);

  FreeWorkSpace(&ctrl, &graph);

  if (ctrl.oflags&OFLAG_COMPRESS) { /* Uncompress the ordering */
    if (graph.nvtxs < COMPRESSION_FRACTION*(nvtxs)) { 
      /* construct perm from iperm */
      for (i=0; i<graph.nvtxs; i++)
        perm[iperm[i]] = i; 
      for (l=ii=0; ii<graph.nvtxs; ii++) {
        i = perm[ii];
        for (j=cptr[i]; j<cptr[i+1]; j++)
          iperm[cind[j]] = l++;
      }
    }

    gk_free((void **)&cptr, &cind, LTERM);
  }


  for (i=0; i<nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl.dbglvl, DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

}



/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void MlevelNestedDissectionP(CtrlType *ctrl, GraphType *graph, idxtype *order, idxtype lastvtx, 
                             idxtype npes, idxtype cpos, idxtype *sizes)
{
  idxtype i, j, nvtxs, nbnd, tvwgt, tpwgts2[2];
  idxtype *label, *bndind;
  GraphType lgraph, rgraph;
  float ubfactor;

  nvtxs = graph->nvtxs;

  if (nvtxs == 0) {
    FreeGraph(graph, 0);
    return;
  }

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt, 1);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  if (cpos >= npes-1) 
    ubfactor = ORDER_UNBALANCE_FRACTION;
  else 
    ubfactor = 1.05;


  MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);

  IFSET(ctrl->dbglvl, DBG_SEPINFO, mprintf("Nvtxs: %6D, [%6D %6D %6D]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

  if (cpos < npes-1) {
    sizes[2*npes-2-cpos] = graph->pwgts[2];
    sizes[2*npes-2-(2*cpos+1)] = graph->pwgts[1];
    sizes[2*npes-2-(2*cpos+2)] = graph->pwgts[0];
  }

  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;
  for (i=0; i<nbnd; i++) 
    order[label[bndind[i]]] = --lastvtx;

  SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  FreeGraph(graph, 0);

  if (rgraph.nvtxs > MMDSWITCH || 2*cpos+1 < npes-1) 
    MlevelNestedDissectionP(ctrl, &rgraph, order, lastvtx, npes, 2*cpos+1, sizes);
  else {
    MMDOrder(ctrl, &rgraph, order, lastvtx); 
    FreeGraph(&rgraph, 0);
  }
  if (lgraph.nvtxs > MMDSWITCH || 2*cpos+2 < npes-1) 
    MlevelNestedDissectionP(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs, npes, 2*cpos+2, sizes);
  else {
    MMDOrder(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs); 
    FreeGraph(&lgraph, 0);
  }
}




/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_NodeComputeSeparator(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
           idxtype *adjwgt, idxtype *options, idxtype *sepsize, idxtype *part) 
{
  idxtype i, j, tvwgt, tpwgts[2];
  GraphType graph;
  CtrlType ctrl;

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, 3);
  tvwgt = idxsum(*nvtxs, graph.vwgt, 1);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags  = 0;
  ctrl.pfactor = 0;
  ctrl.nseps = 3;
  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = amin(100, *nvtxs-1);
  ctrl.maxvwgt = 1.5*tvwgt/ctrl.CoarsenTo;

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, 2);

  /*============================================================
   * Perform the bisection
   *============================================================*/ 
  tpwgts[0] = tvwgt/2;
  tpwgts[1] = tvwgt-tpwgts[0];

  MlevelNodeBisectionMultiple(&ctrl, &graph, tpwgts, 1.02);

  *sepsize = graph.pwgts[2];
  idxcopy(*nvtxs, graph.where, part);

  FreeGraph(&graph, 0);

  FreeWorkSpace(&ctrl, &graph);

}



/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_EdgeComputeSeparator(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
           idxtype *adjwgt, idxtype *options, idxtype *sepsize, idxtype *part) 
{
  idxtype i, j, tvwgt, tpwgts[2];
  GraphType graph;
  CtrlType ctrl;

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, 3);
  tvwgt = idxsum(*nvtxs, graph.vwgt, 1);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags  = 0;
  ctrl.pfactor = 0;
  ctrl.nseps = 1;
  ctrl.optype = OP_OEMETIS;
  ctrl.CoarsenTo = amin(100, *nvtxs-1);
  ctrl.maxvwgt = 1.5*tvwgt/ctrl.CoarsenTo;

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, 2);

  /*============================================================
   * Perform the bisection
   *============================================================*/ 
  tpwgts[0] = tvwgt/2;
  tpwgts[1] = tvwgt-tpwgts[0];

  MlevelEdgeBisection(&ctrl, &graph, tpwgts, 1.05);
  ConstructMinCoverSeparator(&ctrl, &graph, 1.05);

  *sepsize = graph.pwgts[2];
  idxcopy(*nvtxs, graph.where, part);

  FreeGraph(&graph, 0);

  FreeWorkSpace(&ctrl, &graph);
}
