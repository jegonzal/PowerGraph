/*
 * coarsen.c
 *
 * This file contains the driving routines for the coarsening process 
 *
 * Started 7/23/97
 * George
 *
 */

#include <metislib.h>


/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *Coarsen2Way(CtrlType *ctrl, GraphType *graph)
{
  idxtype clevel;
  GraphType *cgraph;

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->CoarsenTmr));

  cgraph = graph;

  /* The following is ahack to allow the multiple bisections to go through with correct
     coarsening */
  if (ctrl->CType > 20) {
    clevel = 1;
    ctrl->CType -= 20;
  }
  else
    clevel = 0;

  do {
    IFSET(ctrl->dbglvl, DBG_COARSEN, mprintf("%6D %7D %7D [%D] [%D %D]\n",
          cgraph->nvtxs, cgraph->nedges/2, idxsum(cgraph->nvtxs, cgraph->adjwgtsum, 1)/2,
          ctrl->CoarsenTo, ctrl->maxvwgt, 
          (cgraph->vwgt ? idxsum(cgraph->nvtxs, cgraph->vwgt, 1) : cgraph->nvtxs)));

    if (cgraph->adjwgt) {
      switch (ctrl->CType) {
        case MTYPE_RM:
          Match_RM(ctrl, cgraph);
          break;
        case MTYPE_HEM:
          if (clevel < 1 || cgraph->nedges == 0)
            Match_RM(ctrl, cgraph);
          else
            Match_HEM(ctrl, cgraph);
          break;
        case MTYPE_SHEM:
          if (clevel < 1 || cgraph->nedges == 0)
            Match_RM(ctrl, cgraph);
          else
            Match_SHEM(ctrl, cgraph);
          break;
        case MTYPE_SHEMKWAY:
          if (cgraph->nedges == 0)
            Match_RM(ctrl, cgraph);
          else
            Match_SHEM(ctrl, cgraph);
          break;
        default:
          errexit("Unknown CType: %d\n", ctrl->CType);
      }
    }
    else {
      Match_RM_NVW(ctrl, cgraph);
    }

    cgraph = cgraph->coarser;
    clevel++;

  } while (cgraph->nvtxs > ctrl->CoarsenTo && cgraph->nvtxs < COARSEN_FRACTION2*cgraph->finer->nvtxs && cgraph->nedges > cgraph->nvtxs/2); 

  IFSET(ctrl->dbglvl, DBG_COARSEN, mprintf("%6D %7D %7D [%D] [%D %D]\n",
        cgraph->nvtxs, cgraph->nedges/2, idxsum(cgraph->nvtxs, cgraph->adjwgtsum, 1)/2,
        ctrl->CoarsenTo, ctrl->maxvwgt, 
        (cgraph->vwgt ? idxsum(cgraph->nvtxs, cgraph->vwgt, 1) : cgraph->nvtxs)));

  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->CoarsenTmr));

  return cgraph;
}

