/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mbalance.c
 *
 * This file contains code that is used to forcefully balance either
 * bisections or k-sections
 *
 * Started 7/29/97
 * George
 *
 * $Id: mbalance.c,v 1.2 2002/08/10 06:29:31 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function is the entry point of the bisection balancing algorithms.
**************************************************************************/
void MocBalance2Way(CtrlType *ctrl, GraphType *graph, float *tpwgts, float lbfactor)
{

  if (Compute2WayHLoadImbalance(graph->ncon, graph->npwgts, tpwgts) < lbfactor)
    return;

  MocGeneral2WayBalance(ctrl, graph, tpwgts, lbfactor);

}


/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void MocGeneral2WayBalance(CtrlType *ctrl, GraphType *graph, float *tpwgts, float lbfactor)
{
  idxtype i, ii, j, k, l, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, pass, me, limit, tmp, cnum;
  idxtype *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idxtype *moved, *swaps, *perm, *qnum;
  float *nvwgt, *npwgts, mindiff[MAXNCON], origbal, minbal, newbal;
  PQueueType parts[MAXNCON][2];
  idxtype higain, oldgain, mincut, newcut, mincutorder;
  idxtype qsizes[MAXNCON][2];

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  id = graph->id;
  ed = graph->ed;
  npwgts = graph->npwgts;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  moved = idxwspacemalloc(ctrl, nvtxs);
  swaps = idxwspacemalloc(ctrl, nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  qnum = idxwspacemalloc(ctrl, nvtxs);

  limit = amin(amax(0.01*nvtxs, 15), 100);

  /* Initialize the queues */
  for (i=0; i<ncon; i++) {
    PQueueInit(ctrl, &parts[i][0], nvtxs, PLUS_GAINSPAN+1);
    PQueueInit(ctrl, &parts[i][1], nvtxs, PLUS_GAINSPAN+1);
    qsizes[i][0] = qsizes[i][1] = 0;
  }

  for (i=0; i<nvtxs; i++) {
    qnum[i] = gk_fargmax(ncon, nvwgt+i*ncon);
    qsizes[qnum[i]][where[i]]++;
  }

/*
  mprintf("Weight Distribution:    \t");
  for (i=0; i<ncon; i++) 
    mprintf(" [%D %D]", qsizes[i][0], qsizes[i][1]); 
  mprintf("\n");
*/

  for (from=0; from<2; from++) {
    for (j=0; j<ncon; j++) {
      if (qsizes[j][from] == 0) {
        for (i=0; i<nvtxs; i++) {
          if (where[i] != from)
            continue;

          k = famax2(ncon, nvwgt+i*ncon);
          if (k == j && qsizes[qnum[i]][from] > qsizes[j][from] && nvwgt[i*ncon+qnum[i]] < 1.3*nvwgt[i*ncon+j]) {
            qsizes[qnum[i]][from]--;
            qsizes[j][from]++;
            qnum[i] = j;
          }
        }
      }
    }
  }

/*
  mprintf("Weight Distribution (after):\t ");
  for (i=0; i<ncon; i++) 
    mprintf(" [%D %D]", qsizes[i][0], qsizes[i][1]); 
  mprintf("\n");
*/



  for (i=0; i<ncon; i++) 
    mindiff[i] = fabs(tpwgts[0]-npwgts[i]);
  minbal = origbal = Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);
  newcut = mincut = graph->mincut;
  mincutorder = -1;

  if (ctrl->dbglvl&DBG_REFINE) {
    mprintf("Parts: [");
    for (l=0; l<ncon; l++)
      mprintf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    mprintf("] T[%.3f %.3f], Nv-Nb[%5D, %5D]. ICut: %6D, LB: %.3f [B]\n", tpwgts[0], tpwgts[1], graph->nvtxs, graph->nbnd, graph->mincut, origbal);
  }

  idxset(nvtxs, -1, moved);

  ASSERT(ComputeCut(graph, where) == graph->mincut);
  ASSERT(CheckBnd(graph));

  /* Insert all nodes in the priority queues */
  nbnd = graph->nbnd;
  RandomPermute(nvtxs, perm, 1);
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    PQueueInsert(&parts[qnum[i]][where[i]], i, ed[i]-id[i]);
  }

  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if (minbal < lbfactor)
      break;

    SelectQueue(ncon, npwgts, tpwgts, &from, &cnum, parts);
    to = (from+1)%2;

    if (from == -1 || (higain = PQueueGetMax(&parts[cnum][from])) == -1)
      break;

    gk_faxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    gk_faxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
    newcut -= (ed[higain]-id[higain]);
    newbal = Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

    if (newbal < minbal || (newbal == minbal && 
        (newcut < mincut || (newcut == mincut && BetterBalance(ncon, npwgts, tpwgts, mindiff))))) {
      mincut = newcut;
      minbal = newbal;
      mincutorder = nswaps;
      for (i=0; i<ncon; i++)
        mindiff[i] = fabs(tpwgts[0]-npwgts[i]);
    }
    else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
      newcut += (ed[higain]-id[higain]);
      gk_faxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
      gk_faxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      break;
    }

    where[higain] = to;
    moved[higain] = nswaps;
    swaps[nswaps] = higain;

    if (ctrl->dbglvl&DBG_MOVEINFO) {
      mprintf("Moved %6D from %D(%D). Gain: %5D, Cut: %5D, NPwgts: ", higain, from, cnum, ed[higain]-id[higain], newcut);
      for (l=0; l<ncon; l++) 
        mprintf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
      mprintf(", %.3f LB: %.3f\n", minbal, newbal);
    }


    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1]) 
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];
      oldgain = ed[k]-id[k];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (moved[k] == -1)
        PQueueUpdate(&parts[qnum[k]][where[k]], k, oldgain, ed[k]-id[k]);

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1) 
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)  
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }



  /****************************************************************
  * Roll back computations
  *****************************************************************/
  for (nswaps--; nswaps>mincutorder; nswaps--) {
    higain = swaps[nswaps];

    to = where[higain] = (where[higain]+1)%2;
    SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    else if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    gk_faxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    gk_faxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      if (bndptr[k] != -1 && ed[k] == 0)
        BNDDelete(nbnd, bndind, bndptr, k);
      if (bndptr[k] == -1 && ed[k] > 0)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }

  if (ctrl->dbglvl&DBG_REFINE) {
    mprintf("\tMincut: %6D at %5D, NBND: %6D, NPwgts: [", mincut, mincutorder, nbnd);
    for (l=0; l<ncon; l++)
      mprintf("(%.3f, %.3f) ", npwgts[l], npwgts[ncon+l]);
    mprintf("], LB: %.3f\n", Compute2WayHLoadImbalance(ncon, npwgts, tpwgts));
  }

  graph->mincut = mincut;
  graph->nbnd = nbnd;


  for (i=0; i<ncon; i++) {
    PQueueFree(ctrl, &parts[i][0]);
    PQueueFree(ctrl, &parts[i][1]);
  }

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);

}

