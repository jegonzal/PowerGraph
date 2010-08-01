/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mrkmetis.c
 *
 * This file contains the top level routines for the multilevel k-way partitioning
 * algorithm KMETIS.
 *
 * Started 7/28/97
 * George
 *
 * $Id: mrkmetis.c,v 1.2 2003/04/04 23:25:10 karypis Exp $
 *
 */

#include <metislib.h>



/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void METIS_mCRefineGraphKway(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, 
                          idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, 
                          idxtype *nparts, float *rubvec, idxtype *options, idxtype *edgecut, 
                          idxtype *part)
{
  idxtype i, j;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_KMETIS, *nvtxs, *ncon, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType  = McKMETIS_CTYPE;
    ctrl.IType  = McKMETIS_ITYPE;
    ctrl.RType  = McKMETIS_RTYPE;
    ctrl.dbglvl = McKMETIS_DBGLVL;
  }
  else {
    ctrl.CType  = options[OPTION_CTYPE];
    ctrl.IType  = options[OPTION_ITYPE];
    ctrl.RType  = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype    = OP_KMETIS;
  ctrl.CoarsenTo = amax((*nvtxs)/(20*gk_log2(*nparts)), 30*(*nparts));
  ctrl.nmaxvwgt  = 0.0; /* GK-MOD: Ensure that no coarsening will take place */

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

  *edgecut = MCMlevelKWayRefinement(&ctrl, &graph, *nparts, part, rubvec);

  IFSET(ctrl.dbglvl, DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
idxtype MCMlevelKWayRefinement(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, 
      float *rubvec)
{
  idxtype i, j, nvtxs;
  GraphType *cgraph;
  idxtype options[10], edgecut;

  cgraph = MCCoarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
  MocAllocateKWayPartitionMemory(ctrl, cgraph, nparts);

  if (cgraph->nvtxs != graph->nvtxs)
    errexit("GK-MOD Failed: %d %d\n", cgraph->nvtxs, graph->nvtxs);

  for (i=0; i<graph->nvtxs; i++)
    cgraph->where[graph->cmap[i]] = part[i];

  MocRefineKWayHorizontal(ctrl, graph, cgraph, nparts, rubvec);

  idxcopy(graph->nvtxs, graph->where, part);

  FreeGraph(graph, 0);

  return graph->mincut;

}

