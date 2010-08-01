/*
 * mcoarsen.c
 *
 * This file contains the driving routines for the coarsening process 
 *
 * Started 7/23/97
 * George
 *
 * $Id: mcoarsen.c,v 1.3 2003/07/31 15:52:47 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *MCCoarsen2Way(CtrlType *ctrl, GraphType *graph)
{
  idxtype i, clevel;
  GraphType *cgraph;

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->CoarsenTmr));

  cgraph = graph;

  clevel = 0;
  do {
    if (ctrl->dbglvl&DBG_COARSEN) {
      mprintf("%6D %7D %10D [%D] [%6.4f", cgraph->nvtxs, cgraph->nedges, 
              idxsum(cgraph->nvtxs, cgraph->adjwgtsum, 1), ctrl->CoarsenTo, ctrl->nmaxvwgt);
      for (i=0; i<graph->ncon; i++)
        mprintf(" %5.3f", gk_fsum(cgraph->nvtxs, cgraph->nvwgt+i, cgraph->ncon));
      mprintf("]\n");
    }

    if (cgraph->nedges == 0) {
      MCMatch_RM(ctrl, cgraph);
    }
    else {
      switch (ctrl->CType) {
        case MTYPE_RM:
          MCMatch_RM(ctrl, cgraph);
          break;
        case MTYPE_HEM:
          if (clevel < 1)
            MCMatch_RM(ctrl, cgraph);
          else
            MCMatch_HEM(ctrl, cgraph);
          break;
        case MTYPE_SHEM:
          if (clevel < 1)
            MCMatch_RM(ctrl, cgraph);
          else
            MCMatch_SHEM(ctrl, cgraph);
          break;
        case MTYPE_SHEMKWAY:
          MCMatch_SHEM(ctrl, cgraph);
          break;
        case MTYPE_SHEBM_ONENORM:
          MCMatch_SHEBM(ctrl, cgraph, 1);
          break;
        case MTYPE_SHEBM_INFNORM:
          MCMatch_SHEBM(ctrl, cgraph, -1);
          break;
        case MTYPE_SBHEM_ONENORM:
          MCMatch_SBHEM(ctrl, cgraph, 1);
          break;
        case MTYPE_SBHEM_INFNORM:
          MCMatch_SBHEM(ctrl, cgraph, -1);
          break;
        default:
          errexit("Unknown CType: %d\n", ctrl->CType);
      }
    }

    cgraph = cgraph->coarser;
    clevel++;

  } while (cgraph->nvtxs > ctrl->CoarsenTo && cgraph->nvtxs < COARSEN_FRACTION2*cgraph->finer->nvtxs && cgraph->nedges > cgraph->nvtxs/2); 

  if (ctrl->dbglvl&DBG_COARSEN) {
    mprintf("%6D %7D %10D [%D] [%6.4f", cgraph->nvtxs, cgraph->nedges, 
            idxsum(cgraph->nvtxs, cgraph->adjwgtsum, 1), ctrl->CoarsenTo, ctrl->nmaxvwgt);
    for (i=0; i<graph->ncon; i++)
      mprintf(" %5.3f", gk_fsum(cgraph->nvtxs, cgraph->nvwgt+i, cgraph->ncon));
    mprintf("]\n");
  }


  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->CoarsenTmr));

  return cgraph;
}

