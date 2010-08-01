/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initpart.c
 *
 * This file contains code that performs the initial partition of the
 * coarsest graph
 *
 * Started 7/23/97
 * George
 *
 */

#include <metislib.h>

/*************************************************************************
* This function computes the initial bisection of the coarsest graph
**************************************************************************/
void Init2WayPartition(CtrlType *ctrl, GraphType *graph, idxtype *tpwgts, float ubfactor) 
{
  idxtype dbglvl;

  dbglvl = ctrl->dbglvl;
  IFSET(ctrl->dbglvl, DBG_REFINE, ctrl->dbglvl -= DBG_REFINE);
  IFSET(ctrl->dbglvl, DBG_MOVEINFO, ctrl->dbglvl -= DBG_MOVEINFO);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));

  switch (ctrl->IType) {
    case ITYPE_GGPKL:
      if (graph->nedges == 0)
        RandomBisection(ctrl, graph, tpwgts, ubfactor);
      else
        GrowBisection(ctrl, graph, tpwgts, ubfactor);
      break;
    case ITYPE_RANDOM:
      RandomBisection(ctrl, graph, tpwgts, ubfactor);
      break;
    default:
      errexit("Unknown initial partition type: %d\n", ctrl->IType);
  }

  IFSET(ctrl->dbglvl, DBG_IPART, mprintf("Initial Cut: %D\n", graph->mincut));
  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));
  ctrl->dbglvl = dbglvl;

/*
  IsConnectedSubdomain(ctrl, graph, 0);
  IsConnectedSubdomain(ctrl, graph, 1);
*/
}

/*************************************************************************
* This function computes the initial bisection of the coarsest graph
**************************************************************************/
void InitSeparator(CtrlType *ctrl, GraphType *graph, float ubfactor) 
{
  idxtype dbglvl;

  dbglvl = ctrl->dbglvl;
  IFSET(ctrl->dbglvl, DBG_REFINE, ctrl->dbglvl -= DBG_REFINE);
  IFSET(ctrl->dbglvl, DBG_MOVEINFO, ctrl->dbglvl -= DBG_MOVEINFO);

  IFSET(ctrl->dbglvl, DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));

  GrowBisectionNode(ctrl, graph, ubfactor);
  Compute2WayNodePartitionParams(ctrl, graph);

  IFSET(ctrl->dbglvl, DBG_IPART, mprintf("Initial Sep: %D\n", graph->mincut));
  IFSET(ctrl->dbglvl, DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));

  ctrl->dbglvl = dbglvl;

}



/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void GrowBisection(CtrlType *ctrl, GraphType *graph, idxtype *tpwgts, float ubfactor)
{
  idxtype i, j, k, nvtxs, drain, nleft, first, last, pwgts[2], minpwgt[2], 
          maxpwgt[2], from, bestcut, icut, mincut, me, pass, nbfs, inbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where;
  idxtype *queue, *touched, *gain, *bestwhere;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  Allocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  queue = idxmalloc(nvtxs, "BisectGraph: queue");
  touched = idxmalloc(nvtxs, "BisectGraph: touched");

  ASSERTP(tpwgts[0]+tpwgts[1] == idxsum(nvtxs, vwgt, 1), ("%d %d\n", tpwgts[0]+tpwgts[1], idxsum(nvtxs, vwgt, 1)));

  maxpwgt[0] = ubfactor*tpwgts[0];
  maxpwgt[1] = ubfactor*tpwgts[1];
  minpwgt[0] = (1.0/ubfactor)*tpwgts[0];
  minpwgt[1] = (1.0/ubfactor)*tpwgts[1];

  nbfs = (nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  for (inbfs=0; inbfs<nbfs; inbfs++) {
    idxset(nvtxs, 0, touched);

    pwgts[1] = tpwgts[0]+tpwgts[1];
    pwgts[0] = 0;

    idxset(nvtxs, 1, where);

    queue[0] = RandomInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; last = 1;
    nleft = nvtxs-1;
    drain = 0;

    /* Start the BFS from queue to get a partition */
    for (;;) {
      if (first == last) { /* Empty. Disconnected graph! */
        if (nleft == 0 || drain)
          break;

        k = RandomInRange(nleft);
        for (i=0; i<nvtxs; i++) {
          if (touched[i] == 0) {
            if (k == 0)
              break;
            else
              k--;
          }
        }

        queue[0] = i;
        touched[i] = 1;
        first = 0; last = 1;;
        nleft--;
      }

      i = queue[first++];
      if (pwgts[0] > 0 && pwgts[1]-vwgt[i] < minpwgt[1]) {
        drain = 1;
        continue;
      }

      where[i] = 0;
      INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
      if (pwgts[1] <= maxpwgt[1])
        break;

      drain = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (touched[k] == 0) {
          queue[last++] = k;
          touched[k] = 1;
          nleft--;
        }
      }
    }

    /* Check to see if we hit any bad limiting cases */
    if (pwgts[1] == 0) { 
      i = RandomInRange(nvtxs);
      where[i] = 1;
      INC_DEC(pwgts[1], pwgts[0], vwgt[i]);
    }

    /*************************************************************
    * Do some partition refinement 
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    /*mprintf("IPART: %3D [%5D %5D] [%5D %5D] %5D\n", graph->nvtxs, pwgts[0], pwgts[1], graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    Balance2Way(ctrl, graph, tpwgts, ubfactor);
    /*mprintf("BPART: [%5D %5D] %5D\n", graph->pwgts[0], graph->pwgts[1], graph->mincut);*/

    FM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);
    /*mprintf("RPART: [%5D %5D] %5D\n", graph->pwgts[0], graph->pwgts[1], graph->mincut);*/

    if (inbfs == 0 || bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  gk_free((void **)&bestwhere, &queue, &touched, LTERM);
}




/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void GrowBisectionNode(CtrlType *ctrl, GraphType *graph, float ubfactor)
{
  idxtype i, j, k, nvtxs, drain, nleft, first, last, pwgts[2], tpwgts[2], 
          minpwgt[2], maxpwgt[2], from, bestcut, icut, mincut, me, pass, 
          nbfs, inbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where, *bndind;
  idxtype *queue, *touched, *gain, *bestwhere;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  queue = idxmalloc(nvtxs, "BisectGraph: queue");
  touched = idxmalloc(nvtxs, "BisectGraph: touched");

  tpwgts[0] = idxsum(nvtxs, vwgt, 1);
  tpwgts[1] = tpwgts[0]/2;
  tpwgts[0] -= tpwgts[1];

  maxpwgt[0] = ubfactor*tpwgts[0];
  maxpwgt[1] = ubfactor*tpwgts[1];
  minpwgt[0] = (1.0/ubfactor)*tpwgts[0];
  minpwgt[1] = (1.0/ubfactor)*tpwgts[1];

  /* Allocate refinement memory. Allocate sufficient memory for both edge and node */
  graph->pwgts  = idxmalloc(3, "GrowBisectionNode: pwgts");
  graph->where  = idxmalloc(nvtxs, "GrowBisectionNode: where");
  graph->bndptr = idxmalloc(nvtxs, "GrowBisectionNode: bndptr");
  graph->bndind = idxmalloc(nvtxs, "GrowBisectionNode: bndind");
  graph->id     = idxmalloc(nvtxs, "GrowBisectionNode: id");
  graph->ed     = idxmalloc(nvtxs, "GrowBisectionNode: ed");
  graph->nrinfo = (NRInfoType *)gk_malloc(nvtxs*sizeof(NRInfoType), "GrowBisectionNode: nrinfo");
  

  where = graph->where;
  bndind = graph->bndind;

  nbfs = (nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS) + 1;
  for (inbfs=0; inbfs<nbfs; inbfs++) {
    idxset(nvtxs, 0, touched);

    pwgts[1] = tpwgts[0]+tpwgts[1];
    pwgts[0] = 0;

    idxset(nvtxs, 1, where);

    queue[0] = RandomInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; last = 1;
    nleft = nvtxs-1;
    drain = 0;

    /* Start the BFS from queue to get a partition */
    if (nbfs >= 1) {
      for (;;) {
        if (first == last) { /* Empty. Disconnected graph! */
          if (nleft == 0 || drain)
            break;
  
          k = RandomInRange(nleft);
          for (i=0; i<nvtxs; i++) {
            if (touched[i] == 0) {
              if (k == 0)
                break;
              else
                k--;
            }
          }

          queue[0] = i;
          touched[i] = 1;
          first = 0; last = 1;;
          nleft--;
        }

        i = queue[first++];
        if (pwgts[1]-vwgt[i] < minpwgt[1]) {
          drain = 1;
          continue;
        }

        where[i] = 0;
        INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
        if (pwgts[1] <= maxpwgt[1])
          break;

        drain = 0;
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (touched[k] == 0) {
            queue[last++] = k;
            touched[k] = 1;
            nleft--;
          }
        }
      }
    }

    /*************************************************************
    * Do some partition refinement 
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    Balance2Way(ctrl, graph, tpwgts, ubfactor);
    FM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);

    /* Construct and refine the vertex separator */
    for (i=0; i<graph->nbnd; i++) 
      where[bndind[i]] = 2;

    Compute2WayNodePartitionParams(ctrl, graph); 
    FM_2WayNodeRefine(ctrl, graph, ubfactor, 6);

    /* mprintf("ISep: [%D %D %D] %D\n", graph->pwgts[0], graph->pwgts[1], graph->pwgts[2], bestcut); */

    if (inbfs == 0 || bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  Compute2WayNodePartitionParams(ctrl, graph); 

  gk_free((void **)&bestwhere, &queue, &touched, LTERM);
}


/*************************************************************************
* This function takes a graph and produces a bisection by using a region
* growing algorithm. The resulting partition is returned in
* graph->where
**************************************************************************/
void RandomBisection(CtrlType *ctrl, GraphType *graph, idxtype *tpwgts, float ubfactor)
{
  idxtype i, ii, j, k, nvtxs, pwgts[2], minpwgt[2], maxpwgt[2], from, bestcut, 
          icut, mincut, me, pass, nbfs, inbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where;
  idxtype *perm, *bestwhere;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  Allocate2WayPartitionMemory(ctrl, graph);
  where = graph->where;

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  perm = idxmalloc(nvtxs, "BisectGraph: queue");

  ASSERTP(tpwgts[0]+tpwgts[1] == idxsum(nvtxs, vwgt, 1), ("%d %d\n", tpwgts[0]+tpwgts[1], idxsum(nvtxs, vwgt, 1)));

  maxpwgt[0] = ubfactor*tpwgts[0];
  maxpwgt[1] = ubfactor*tpwgts[1];
  minpwgt[0] = (1.0/ubfactor)*tpwgts[0];
  minpwgt[1] = (1.0/ubfactor)*tpwgts[1];

  nbfs = (nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);
  for (inbfs=0; inbfs<nbfs; inbfs++) {
    RandomPermute(nvtxs, perm, 1);

    idxset(nvtxs, 1, where);
    pwgts[1] = tpwgts[0]+tpwgts[1];
    pwgts[0] = 0;


    if (nbfs != 1) {
      for (ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        if (pwgts[0]+vwgt[i] < maxpwgt[0]) {
          where[i] = 0;
          pwgts[0] += vwgt[i];
          pwgts[1] -= vwgt[i];
          if (pwgts[0] > minpwgt[0])
            break;
        }
      }
    }

    /*************************************************************
    * Do some partition refinement 
    **************************************************************/
    Compute2WayPartitionParams(ctrl, graph);
    /* mprintf("IPART: %3D [%5D %5D] [%5D %5D] %5D\n", graph->nvtxs, pwgts[0], pwgts[1], graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    Balance2Way(ctrl, graph, tpwgts, ubfactor);
    /* mprintf("BPART: [%5D %5D] %5D\n", graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    FM_2WayEdgeRefine(ctrl, graph, tpwgts, 4);
    /* mprintf("RPART: [%5D %5D] %5D\n", graph->pwgts[0], graph->pwgts[1], graph->mincut); */

    if (inbfs==0 || bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
      if (bestcut == 0)
        break;
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  gk_free((void **)&bestwhere, &perm, LTERM);
}




