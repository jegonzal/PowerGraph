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
 */

#include <metislib.h>


/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void METIS_PartGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                         idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, 
                         idxtype *options, idxtype *edgecut, idxtype *part)
{
  idxtype i;
  float *tpwgts;

  tpwgts = gk_fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++) 
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, 
                       tpwgts, options, edgecut, part);

  gk_free((void **)&tpwgts, LTERM);
}


/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void METIS_WPartGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
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
  ctrl.CoarsenTo = amax((*nvtxs)/(40*gk_log2(*nparts)), 20*(*nparts));
  ctrl.maxvwgt = 1.5*((graph.vwgt ? idxsum(*nvtxs, graph.vwgt, 1) : (*nvtxs))/ctrl.CoarsenTo);

  InitRandom(-1);

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
* This function takes a graph and produces a bisection of it
**************************************************************************/
idxtype MlevelKWayPartitioning(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, float *tpwgts, float ubfactor)
{
  idxtype i, j, nvtxs, tvwgt, tpwgts2[2];
  GraphType *cgraph;
  idxtype wgtflag=3, numflag=0, options[10], edgecut;

  cgraph = Coarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
  AllocateKWayPartitionMemory(ctrl, cgraph, nparts);

  options[0] = 1; 
  options[OPTION_CTYPE] = MTYPE_SHEMKWAY;
  options[OPTION_ITYPE] = ITYPE_GGPKL;
  options[OPTION_RTYPE] = RTYPE_FM;
  options[OPTION_DBGLVL] = 0;

  METIS_WPartGraphRecursive(&cgraph->nvtxs, cgraph->xadj, cgraph->adjncy, cgraph->vwgt, 
                            cgraph->adjwgt, &wgtflag, &numflag, &nparts, tpwgts, options, 
                            &edgecut, cgraph->where);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));
  IFSET(ctrl->dbglvl, DBG_IPART, mprintf("Initial %D-way partitioning cut: %D\n", nparts, edgecut));

  IFSET(ctrl->dbglvl, DBG_KWAYPINFO, ComputePartitionInfo(cgraph, nparts, cgraph->where));

  RefineKWay(ctrl, graph, cgraph, nparts, tpwgts, ubfactor);

  idxcopy(graph->nvtxs, graph->where, part);

  FreeGraph(graph, 0);

  return graph->mincut;

}

