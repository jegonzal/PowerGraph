/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mkmetis.c
 *
 * This file contains the top level routines for the multilevel k-way partitioning
 * algorithm KMETIS.
 *
 * Started 7/28/97
 * George
 *
 * $Id: mkmetis.c,v 1.3 2002/08/10 07:07:27 karypis Exp $
 *
 */

#include <metislib.h>



/*************************************************************************
* This function is the entry point for KWMETIS
**************************************************************************/
void METIS_mCPartGraphKway(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, 
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
  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = amax((*nvtxs)/(20*gk_log2(*nparts)), 30*(*nparts));

  ctrl.nmaxvwgt = 1.5/(1.0*ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

  *edgecut = MCMlevelKWayPartitioning(&ctrl, &graph, *nparts, part, rubvec);

  IFSET(ctrl.dbglvl, DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
idxtype MCMlevelKWayPartitioning(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, 
      float *rubvec)
{
  idxtype i, j, nvtxs;
  GraphType *cgraph;
  idxtype options[10], edgecut;

  cgraph = MCCoarsen2Way(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
  MocAllocateKWayPartitionMemory(ctrl, cgraph, nparts);

  options[0] = 1; 
  options[OPTION_CTYPE] = MTYPE_SBHEM_INFNORM;
  options[OPTION_ITYPE] = ITYPE_RANDOM;
  options[OPTION_RTYPE] = RTYPE_FM;
  options[OPTION_DBGLVL] = 0;

  /* Determine what you will use as the initial partitioner, based on tolerances */
  for (i=0; i<graph->ncon; i++) {
    if (rubvec[i] > 1.2)
      break;
  }
  if (i == graph->ncon)
    METIS_mCPartGraphRecursiveInternal(&cgraph->nvtxs, &cgraph->ncon, 
          cgraph->xadj, cgraph->adjncy, cgraph->nvwgt, cgraph->adjwgt, &nparts, 
          options, &edgecut, cgraph->where);
  else
    METIS_mCHPartGraphRecursiveInternal(&cgraph->nvtxs, &cgraph->ncon, 
          cgraph->xadj, cgraph->adjncy, cgraph->nvwgt, cgraph->adjwgt, &nparts, 
          rubvec, options, &edgecut, cgraph->where);


  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));
  IFSET(ctrl->dbglvl, DBG_IPART, mprintf("Initial %D-way partitioning cut: %D\n", nparts, edgecut));

  IFSET(ctrl->dbglvl, DBG_KWAYPINFO, ComputePartitionInfo(cgraph, nparts, cgraph->where));

  MocRefineKWayHorizontal(ctrl, graph, cgraph, nparts, rubvec);

  idxcopy(graph->nvtxs, graph->where, part);

  FreeGraph(graph, 0);

  return graph->mincut;

}

