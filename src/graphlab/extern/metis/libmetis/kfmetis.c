/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kfmetis.c
 *
 * This file contains the top level routines for computing a k-way partitioning
 * that also balances the fill-in.
 *
 * Started 8/12/02
 * George
 *
 */

#include <metislib.h>


/*************************************************************************
* This function is the entry point for KFMETIS
**************************************************************************/
void METIS_PartFillGraph(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                         idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, 
                         idxtype *options, idxtype *edgecut, idxtype *part)
{
  idxtype i, j;
  GraphType graph;
  CtrlType ctrl;

  /* Compute an initial k-way partitioning of the graph */
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, 
                      options, edgecut, part);

  /*-----------------------------------------------------------------------
   * Internalize the graph and allocate memory for the workspace 
   *-----------------------------------------------------------------------*/
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


  /*-----------------------------------------------------------------------
   * Balance the fill-in 
   *-----------------------------------------------------------------------*/
  BalanceFillIn(&ctrl, &graph, *nparts, part);


  *edgecut = ComputeCut(&graph, part);

  /*-----------------------------------------------------------------------
   * Post-process the results
   *-----------------------------------------------------------------------*/
  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}




/**************************************************************************
* This function is the top-level routine of the fill-balancing code 
***************************************************************************/
void BalanceFillIn(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part)
{
  idxtype i, j, k, nvtxs, options[10];
  idxtype *fills, *subfills, *ssize, *spart, *vmap, *ppart, *vimap;
  idxtype oldavg, oldmax, oldcut, oldmaxcut, newavg, newmax, newcut, newmaxcut;
  idxtype *tssize;
  float diff, alpha;
  GraphType *pgraph;

  nvtxs = graph->nvtxs;

  fills    = idxsmalloc(nparts, 0, "BalanceFillIn: fills");
  subfills = idxsmalloc(2*nparts, 0, "BalanceFillIn: subfills");
  ssize    = idxsmalloc(nparts, 0, "BalanceFillIn: ssize");
  tssize   = idxsmalloc(nparts, 0, "BalanceFillIn: tssize");

  spart    = idxmalloc(nvtxs, "BalanceFillIn: spart");
  ppart    = idxmalloc(nvtxs, "BalanceFillIn: ppart");
  vmap     = idxmalloc(nvtxs, "BalanceFillIn: vmap");
  vimap    = idxmalloc(nvtxs, "BalanceFillIn: vimap");

  /* Go and compute a top-level vertex separator for each partition */
  for (i=0; i<nparts; i++) {
    pgraph = ExtractPartitionGraph(graph, part, i, vmap, vimap);

    options[0] = 0;
    options[7] = -1;
    METIS_NodeComputeSeparator(&(pgraph->nvtxs), pgraph->xadj, pgraph->adjncy, pgraph->vwgt, 
                               pgraph->adjwgt, options, &ssize[i], ppart);


    /* Copy the separator and partition information to the global graph's spart vector */
    for (j=0; j<pgraph->nvtxs; j++) 
      spart[vimap[j]] = ppart[j];

    FreeGraph(pgraph, 1);
  }
    
    
  /* Go and compute the fill-ins of each partition */
  for (i=0; i<nparts; i++) 
    ComputePartitionFillIn(ctrl, graph, nparts, part, spart, i, &fills[i], &subfills[2*i]);

  oldavg    = idxsum(nparts, fills, 1)/nparts;
  oldmax    = fills[idxargmax(nparts, fills)];
  oldcut    = ComputeCut(graph, part);
  oldmaxcut = ComputeMaxCut(graph, nparts, part);

  mprintf("AverageFill: %10D, MaxFill: %10D, LoadImbalance: %.3f, EdgeCut: %5D, MaxCut: %5D\n", 
         oldavg, oldmax, 1.0*oldmax/oldavg, oldcut, oldmaxcut);

  /* Determine the target separator sizes of the domains that are above average */
  for (i=0; i<nparts; i++) {
    if (fills[i] < 1.05*oldavg) {
      tssize[i] = ssize[i];
    }
    else {
      if (oldavg - (subfills[2*i] +subfills[2*i+1])< 0) {
        tssize[i] = ssize[i]/5.0;
        mprintf("* %2D %7D %7D %6D %6D\n", i, fills[i], subfills[2*i]+subfills[2*i+1], ssize[i], tssize[i]);
      }
      else {
        alpha = sqrt(fills[i] - (subfills[2*i]+subfills[2*i+1]))/ssize[i];
        tssize[i] = sqrt(oldavg - (subfills[2*i]+subfills[2*i+1]))/alpha;
        mprintf("< %2D %7D %7D %6D %6D\n", i, fills[i], subfills[2*i]+subfills[2*i+1], ssize[i], tssize[i]);
      }
    }
  }

      
  /* Perform a simple greedy refinement */
  RefineTopLevelSeparators(ctrl, graph, nparts, part, spart, fills, subfills, ssize, tssize);

  /* Go and compute the fill-ins of each partition */
  for (i=0; i<nparts; i++) 
    ComputePartitionFillIn(ctrl, graph, nparts, part, spart, i, &fills[i], &subfills[2*i]);

for (i=0; i<nparts; i++) {
  if (tssize[i] < ssize[i])
    mprintf("%3D %6D %6D %8D\n", i, ssize[i], tssize[i], fills[i]);
}

  newavg    = idxsum(nparts, fills, 1)/nparts;
  newmax    = fills[idxargmax(nparts, fills)];
  newcut    = ComputeCut(graph, part);
  newmaxcut = ComputeMaxCut(graph, nparts, part);

  mprintf("AverageFill: %10D, MaxFill: %10D, LoadImbalance: %.3f, EdgeCut: %5D, MaxCut: %5D\n", 
         newavg, newmax, 1.0*newmax/newavg, newcut, newmaxcut);

  mprintf("Effective Load Imbalance: %.3f, Sum Fill/Comm Ratio: %6.2f, Max Fill/Comm Ratio: %6.2f\n", 
          1.0*newmax/oldavg, 1.0*(oldmax-newmax)/(newcut-oldcut), 1.0*(oldmax-newmax)/(newmaxcut-oldmaxcut));

  gk_free((void **)&fills, &subfills, &ssize, &spart, &ppart, &vmap, &vimap, LTERM);

}



/*************************************************************************
* This function computes a fill-reducing ordering for a particular 
* partition, whose top-level separator is provided as input
**************************************************************************/
void ComputePartitionFillIn(CtrlType *ctrl, GraphType *graph, idxtype nparts, 
       idxtype *part, idxtype *spart, idxtype pid, idxtype *r_fill, 
       idxtype *r_subfill)
{
  idxtype i, j, k, pnvtxs, options[10], numflag=0;
  idxtype *vmap, *vimap, *ppart, *iperm;
  idxtype *lvmap, *lvimap, *rvmap, *rvimap, *lperm, *liperm, *rperm, *riperm;
  GraphType *pgraph, *lgraph, *rgraph;

  /* Extract the partition graph */
  vmap   = idxmalloc(graph->nvtxs, "ComputePartitionFillIn: vmap");
  vimap  = idxmalloc(graph->nvtxs, "ComputePartitionFillIn: vimap");

  pgraph = ExtractPartitionGraph(graph, part, pid, vmap, vimap);

  pnvtxs  = pgraph->nvtxs;

  iperm  = idxmalloc(pnvtxs, "ComputePartitionFillIn: iperm");
  ppart  = idxmalloc(pnvtxs, "ComputePartitionFillIn: ppart");

  lvmap  = idxmalloc(pnvtxs, "ComputePartitionFillIn: lvmap");
  lvimap = idxmalloc(pnvtxs, "ComputePartitionFillIn: lvimap");
  rvmap  = idxmalloc(pnvtxs, "ComputePartitionFillIn: rvmap");
  rvimap = idxmalloc(pnvtxs, "ComputePartitionFillIn: rvimap");
  lperm  = idxmalloc(pnvtxs, "ComputePartitionFillIn: liperm");
  liperm = idxmalloc(pnvtxs, "ComputePartitionFillIn: liperm");
  rperm  = idxmalloc(pnvtxs, "ComputePartitionFillIn: riperm");
  riperm = idxmalloc(pnvtxs, "ComputePartitionFillIn: riperm");


  /* project down the separator-partition */
  for (i=0; i<pnvtxs; i++)
    ppart[i] = spart[vimap[i]];

  /* Check to see if the separators are indeed separators */
  for (i=0; i<pnvtxs; i++) {
    for (j=pgraph->xadj[i]; j<pgraph->xadj[i+1]; j++) {
      if (ppart[i] != ppart[pgraph->adjncy[j]]) {
        if (!(ppart[i] == 2 || ppart[pgraph->adjncy[j]] == 2)) 
          mprintf("Error: %D %D %D\n", i, ppart[i], ppart[pgraph->adjncy[j]]);
      }
    }
  }
        

  lgraph = ExtractPartitionGraph(pgraph, ppart, 0, lvmap, lvimap);
  rgraph = ExtractPartitionGraph(pgraph, ppart, 1, rvmap, rvimap);

  options[0] = 1;
  options[1] = 3;
  options[2] = 1;
  options[3] = 1;
  options[4] = 0;
  options[5] = 0;
  options[6] = 0;
  options[7] = 4;
  METIS_NodeND(&lgraph->nvtxs, lgraph->xadj, lgraph->adjncy, &numflag, options, lperm, liperm);
  METIS_NodeND(&rgraph->nvtxs, rgraph->xadj, rgraph->adjncy, &numflag, options, rperm, riperm);
    
  /* Put the ordering together */
  for (i=0; i<lgraph->nvtxs; i++) 
    iperm[lvimap[i]] = liperm[i];
  for (i=0; i<rgraph->nvtxs; i++) 
    iperm[rvimap[i]] = lgraph->nvtxs + riperm[i];

  for (j=lgraph->nvtxs+rgraph->nvtxs, i=0; i<pnvtxs; i++) {
    if (ppart[i] == 2) 
      iperm[i] = j++;
  }


  *r_fill      = ComputeFillIn2(pgraph, iperm);
  r_subfill[0] = ComputeFillIn2(lgraph, liperm);
  r_subfill[1] = ComputeFillIn2(rgraph, riperm);

  mprintf("%4D %5D %5D %5D %10D %10D %10D\n", pid, pnvtxs-lgraph->nvtxs-rgraph->nvtxs, 
          lgraph->nvtxs, rgraph->nvtxs, r_subfill[0], r_subfill[1], *r_fill);

  FreeGraph(pgraph, 1);
  FreeGraph(lgraph, 1);
  FreeGraph(rgraph, 1);

  gk_free((void **)&vmap, &vimap, &iperm, &ppart, &lvmap, &lvimap, &rvmap, &rvimap,
         &lperm, &liperm, &rperm, &riperm, LTERM);

}



/*************************************************************************
* This function extracts the graph of a partition and returns it.
**************************************************************************/
GraphType *ExtractPartitionGraph(GraphType *graph, idxtype *part, idxtype pid, 
                                 idxtype *vmap, idxtype *vimap)
{
  idxtype i, j, k, nvtxs, pnvtxs, pnedges;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  idxtype *pxadj, *padjncy, *pvwgt, *padjwgt;
  GraphType *pgraph;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  vwgt   = graph->vwgt;
  adjwgt = graph->adjwgt;


  /* Determine the size of the partition graph */
  for (pnvtxs=pnedges=0, i=0; i<nvtxs; i++) {
    if (part[i] != pid)
      continue;

    vmap[i] = pnvtxs;
    vimap[pnvtxs] = i;
    pnvtxs++;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (part[adjncy[j]] == pid)
        pnedges++;
    }
  }


  /* Allocate memory for the partition graph */
  pgraph = CreateGraph();
  pgraph->nvtxs  = pnvtxs;
  pgraph->nedges = pnedges;

  pxadj   = pgraph->xadj          = idxmalloc(pnvtxs+1, "ExtractPartitionGraph: xadj");
  pvwgt   = pgraph->vwgt          = idxmalloc(pnvtxs,   "ExtractPartitionGraph: vwgt");
  padjncy = pgraph->adjncy        = idxmalloc(pnedges,  "ExtractPartitionGraph: adjncy");
  padjwgt = pgraph->adjwgt        = idxmalloc(pnedges,  "ExtractPartitionGraph: adjwgt");
  pgraph->adjwgtsum               = idxmalloc(pnvtxs,   "ExtractPartitionGraph: adjwgtsum");
  pgraph->cmap                    = idxmalloc(pnvtxs,   "ExtractPartitionGraph: cmap");


  /* Go and extract the partition graph */
  pxadj[0] = pnvtxs = pnedges = 0;
  for (i=0; i<nvtxs; i++) {
    if (part[i] != pid) 
      continue;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (part[adjncy[j]] == pid) {
        padjncy[pnedges] = vmap[adjncy[j]];
        padjwgt[pnedges] = adjwgt[j];
        pnedges++;
      }
    }
    pvwgt[pnvtxs] = vwgt[i]; 
    pxadj[++pnvtxs] = pnedges;
  }

  return pgraph;
}



/*************************************************************************
* This function implements a very simple greedy k-way refinement routine 
* for top level separators
**************************************************************************/
void RefineTopLevelSeparators(CtrlType *ctrl, GraphType *graph, idxtype nparts,  
                              idxtype *part, idxtype *spart, idxtype *fills, idxtype *subfills, 
                              idxtype *ssizes, idxtype *tssizes)
{
  idxtype i, ii, j, k, from, to, nvtxs, mincut, ndegrees, pass, avgssize, 
      avgfill, id, maxpwgt, maxspwgt, nmoves, avgdegree, nskip[5];
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  idxtype *perm, *pwgts, *pdegrees[4], *pmarker, *pids, *spwgts[3];
  float *alphas[3];

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;


  pwgts       = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: pwgts");
  spwgts[0]   = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: spwgts");
  spwgts[1]   = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: spwgts");
  spwgts[2]   = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: spwgts");
  pdegrees[0] = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: pdegrees");
  pdegrees[1] = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: pdegrees");
  pdegrees[2] = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: pdegrees");
  pdegrees[3] = idxsmalloc(nparts, 0, "RefineTopLevelSeparators: pdegrees");
  pmarker     = idxsmalloc(nparts, -1, "RefineTopLevelSeparators: pmarker ");
  pids        = idxmalloc(nparts, "RefineTopLevelSeparators: pids ");
  perm        = idxmalloc(nvtxs, "RefineTopLevelSeparators: perm");
  alphas[0]   = gk_fmalloc(nparts, "RefineTopLevelSeparators: alphas[0]");
  alphas[1]   = gk_fmalloc(nparts, "RefineTopLevelSeparators: alphas[1]");
  alphas[2]   = gk_fmalloc(nparts, "RefineTopLevelSeparators: alphas[2]");

  mincut = ComputeCut(graph, part);
  
  /* Compute partition and separator sizes */
  for (i=0; i<nvtxs; i++) {
    pwgts[part[i]] += vwgt[i];
    spwgts[spart[i]][part[i]] += vwgt[i];
  }

  /* Estimate the exponents for partition-size-fill approximation */
  for (i=0; i<nparts; i++) {
    alphas[0][i] = log(subfills[2*i])/log(spwgts[0][i]);
    alphas[1][i] = log(subfills[2*i+1])/log(spwgts[1][i]);
    alphas[2][i] = log(fills[i] - (subfills[2*i]+subfills[2*i+1]))/log(ssizes[i]);
  }

  maxpwgt = (int)((float)((UNBALANCE_FRACTION*idxsum(nparts, pwgts, 1))/nparts));
  maxspwgt = 1.10*(idxsum(nparts, pwgts, 1)-idxsum(nparts, ssizes, 1))/(2.0*nparts);

  avgssize  = idxsum(nparts, ssizes, 1)/nparts;
  avgfill   = idxsum(nparts, fills, 1)/nparts;
  avgdegree = xadj[nvtxs]/(1.0*nvtxs);


mprintf("Avgssize: %D, Avgfill: %D, AvgDegree: %D\n", avgssize, avgfill, avgdegree);

  /*--------------------------------------------------------------------------
   * Go and do the refinement
   *--------------------------------------------------------------------------*/
  for (pass=0; pass<10; pass++) {
    RandomPermute(nvtxs, perm, 1);

    nskip[0] = nskip[1] = nskip[2] = nskip[3] = nskip[4] = 0;
    for (nmoves=0, ii=0; ii<nvtxs; ii++) {
      i = perm[ii];

      from = part[i];

      /* Determine the connectivity of this vertex */
      for (ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
        to = part[adjncy[j]];
        if ((k = pmarker[to]) == -1) {
          pmarker[to] = k = ndegrees++;
          pids[k] = to;
          pdegrees[0][k] = pdegrees[1][k] = pdegrees[2][k] = pdegrees[3][k] = 0;
        }

        pdegrees[spart[adjncy[j]]][k] += adjwgt[j]; 
        pdegrees[3][k]                += adjwgt[j]; 
      }

      id = (pmarker[from] != -1 ? pdegrees[3][pmarker[from]] : 0);

      if (spart[i] == 2) { /* I'm on the separator */
        for (to=-1, j=0; j<ndegrees; j++) {
          if (pids[j] == from || pwgts[pids[j]]+vwgt[i] > maxpwgt)
            continue;
          
          if (fills[pids[j]] > amax(fills[from], avgfill)) {
            nskip[0]++;
            if (tssizes[from] < ssizes[from])
              nskip[2]++;
            continue;
          }
          if (pdegrees[0][j] > 0 && spwgts[0][pids[j]] > maxspwgt) 
            continue;
          if (pdegrees[1][j] > 0 && spwgts[1][pids[j]] > maxspwgt) 
            continue;

          if (pdegrees[0][j]*pdegrees[1][j] == 0) { /* will not increase target separator */
            to = (to == -1 || pdegrees[3][to] < pdegrees[3][j] ? j : to);
          }
        }

        if (to == -1 && tssizes[from] < ssizes[from])
          nskip[1]++;

        /* See if some of the moves should be skipped */
        if (to != -1 && pdegrees[3][to] < id && tssizes[from] > ssizes[from]) 
          to = -1;

        if (to != -1 && pdegrees[3][to] < id && RandomInRange(3) == 0) {
          nskip[3]++;
          to = -1;
        }

        if (to != -1 && id - pdegrees[3][to] > avgdegree*(0.9+0.07*pass)) {
          nskip[4]++;
          to = -1;
        }
      }
      else { /* I'm not on the separator */
        for (to=-1, j=0; j<ndegrees; j++) {
          if (pids[j] == from || pwgts[pids[j]]+vwgt[i] > maxpwgt)
            continue;
          
          if (tssizes[pids[j]] < ssizes[pids[j]])
            continue;

          if (pdegrees[0][j] > 0 && spwgts[0][pids[j]] > maxspwgt)
            continue;
          if (pdegrees[1][j] > 0 && spwgts[1][pids[j]] > maxspwgt)
            continue;

          if (pdegrees[0][j]*pdegrees[1][j] == 0) { /* will not increase target separator */
            to = (to == -1 || pdegrees[3][to] < pdegrees[3][j] ? j : to);
          }
        }

        if (to != -1) {
          if (pdegrees[3][to] < id)
            to = -1;
          else if (pdegrees[3][to] == id && fills[from] < fills[pids[to]] /*pwgts[from] < pwgts[pids[to]]*/)
            to = -1;
        }
      }

      if (to != -1) {
/*
        mprintf("%2D. Moving %6D %3D -> %3D, %2D -> %2D, Gain: %+4D, SepFrom: %4D, SepTo: %3D, Cut: %5D, %4D %4D\n", 
                pass, i, from, pids[to], spart[i], (pdegrees[0][to] > pdegrees[1][to] ? 0 : 1),
                pdegrees[3][to] - id, ssizes[from], ssizes[pids[to]], mincut, pwgts[from], pwgts[pids[to]]);
*/
        /* Update the approximations of fill of 'from' */
        if (spart[i] == 2) 
          fills[from] -= (pow(ssizes[from], alphas[2][from]) - pow(ssizes[from]-vwgt[i], alphas[2][from]));
        else
          fills[from] -= (pow(spwgts[spart[i]][from], alphas[spart[i]][from]) - pow(spwgts[spart[i]][from]-vwgt[i], alphas[spart[i]][from]));


        if (spart[i] == 2) 
          ssizes[from] -= vwgt[i];

        pwgts[from]            -= vwgt[i];
        spwgts[spart[i]][from] -= vwgt[i];


        part[i]  = pids[to];
        spart[i] = (pdegrees[0][to] > pdegrees[1][to] ? 0 : 1);

        /* Update the approximations of fill of 'to' */
        fills[to] += (pow(spwgts[spart[i]][to]+vwgt[i], alphas[spart[i]][to]) - pow(spwgts[spart[i]][to], alphas[spart[i]][to]));

        pwgts[pids[to]]            += vwgt[i];
        spwgts[spart[i]][pids[to]] += vwgt[i];

        mincut = mincut - (pdegrees[3][to] - id);
        nmoves++;
      }


      /* reset pmarker */
      for (j=xadj[i]; j<xadj[i+1]; j++) 
        pmarker[part[adjncy[j]]] = -1;
    }

    mprintf("%4D. Nmoves: %4D, Cut: %6D [%4D %4D %4D %4D %4D]\n", pass, nmoves, mincut,
            nskip[0], nskip[1], nskip[2], nskip[3], nskip[4]);

    if (nmoves == 0)
      break;
  }

for (i=0; i<nparts; i++)
  mprintf("%4D %6D %8D %8D\n", i, pwgts[i], ssizes[i], fills[i]);

mprintf("%D %D %D / %D %D %D / %D %D\n", pwgts[10], spwgts[0][10], spwgts[1][10], pwgts[15], spwgts[0][15], spwgts[1][15], maxpwgt, maxspwgt);

  gk_free((void **)&pwgts, &pdegrees[0], &pdegrees[1], &pdegrees[2], &pdegrees[3], 
         &pmarker, &pids, &perm, &spwgts[0], &spwgts[1], &spwgts[2], &alphas[0], 
         &alphas[1], &alphas[2], LTERM);
}
