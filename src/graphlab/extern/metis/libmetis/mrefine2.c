/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mrefine2.c
 *
 * This file contains the driving routines for multilevel refinement
 *
 * Started 7/24/97
 * George
 *
 * $Id: mrefine2.c,v 1.2 2002/08/10 06:29:33 karypis Exp $
 */

#include <metislib.h>


/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void MocRefine2Way2(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, float *tpwgts, 
       float *ubvec)
{

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->UncoarsenTmr));

  /* Compute the parameters of the coarsest graph */
  MocCompute2WayPartitionParams(ctrl, graph);

  for (;;) {
    ASSERT(CheckBnd(graph));

    IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->RefTmr));
    switch (ctrl->RType) {
      case RTYPE_FM:
        MocBalance2Way2(ctrl, graph, tpwgts, ubvec);
        MocFM_2WayEdgeRefine2(ctrl, graph, tpwgts, ubvec, 8); 
        break;
      default:
        errexit("Unknown refinement type: %d\n", ctrl->RType);
    }
    IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    graph = graph->finer;
    IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->ProjectTmr));
    MocProject2WayPartition(ctrl, graph);
    IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->ProjectTmr));
  }

  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->UncoarsenTmr));
}


