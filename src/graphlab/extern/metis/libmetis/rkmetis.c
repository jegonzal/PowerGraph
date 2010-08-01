/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This file contains the top level routines for the multilevel k-way partitioning
 * algorithm KMETIS.
 *
 * Started 7/28/97
 * George
 *
 * $Id: rkmetis.c,v 1.3 2003/04/01 05:09:37 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void METIS_RefineGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                           idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, 
                           idxtype *options, idxtype *edgecut, idxtype *part)
{
  idxtype i;
  float *tpwgts;

  tpwgts = gk_fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++) 
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WRefineGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, 
                         tpwgts, options, edgecut, part);

  gk_free((void **)&tpwgts, LTERM);
}


/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void METIS_WRefineGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
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
    ctrl.dbglvl = DBG_REFINE;
    ctrl.dbglvl = KMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = amax((*nvtxs)/(40*gk_log2(*nparts)), 20*(*nparts));
  ctrl.maxvwgt = 0; /* GK-MOD: Ensure that no coarsening will take place */

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

  *edgecut = MlevelKWayRefinement(&ctrl, &graph, *nparts, part, tpwgts, 1.03);

  IFSET(ctrl.dbglvl, DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
idxtype MlevelKWayRefinement(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, 
          float *tpwgts, float ubfactor)
{
  idxtype i, j, nvtxs, tvwgt, tpwgts2[2];
  GraphType *cgraph;
  idxtype wgtflag=3, numflag=0, options[10], edgecut;

  cgraph = Coarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
  AllocateKWayPartitionMemory(ctrl, cgraph, nparts);

  if (cgraph->nvtxs != graph->nvtxs)
    errexit("GK-MOD Failed: %d %d\n", cgraph->nvtxs, graph->nvtxs);

  for (i=0; i<graph->nvtxs; i++)
    cgraph->where[graph->cmap[i]] = part[i];

  RefineKWayRefinement(ctrl, graph, cgraph, nparts, tpwgts, ubfactor);

  idxcopy(graph->nvtxs, graph->where, part);

  FreeGraph(graph, 0);

  return graph->mincut;

}

