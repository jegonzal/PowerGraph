/*
 * Copyright 2003, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This file contains the top level routines for the contact-friendly partitioning
 * algorithm CMETIS.
 *
 * Started 4/3/03
 * George
 *
 */

#include <metislib.h>

#define DEPSILON 1.0e-012

#define _PRINTSTAT

/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void *METIS_PartGraphForContact(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, 
                double *xyzcoords, idxtype *sflag, idxtype *numflag, idxtype *nparts, 
                idxtype *options, idxtype *edgecut, idxtype *part) 
{
  idxtype i, j, ii, dim, ncon, wgtflag, mcnumflag, nnodes, nlnodes, nclean, naclean, ndirty, maxdepth, rwgtflag, rnumflag;
  idxtype *mcvwgt, *dtpart, *marker, *leafpart;
  idxtype *adjwgt;
  float rubvec[2], lbvec[2];
  GraphType graph, *cgraph;
  ContactInfoType *cinfo;
  DKeyValueType *xyzcand[3];

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  /*---------------------------------------------------------------------
   * Allocate memory for the contact info type
   *---------------------------------------------------------------------*/
  cinfo = (ContactInfoType *)gk_malloc(sizeof(ContactInfoType), "METIS_PartGraphForContact: cinfo");
  cinfo->leafptr  = idxsmalloc(*nvtxs+1, 0, "METIS_PartGraphForContact: leafptr");
  cinfo->leafind  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafind");
  cinfo->leafwgt  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafwgt");
  cinfo->part     = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: part");
  leafpart = cinfo->leafpart = idxmalloc(*nvtxs, "METIS_PartGraphForContact: leafpart");
  cinfo->dtree    = (DTreeNodeType *)gk_malloc(sizeof(DTreeNodeType)*(*nvtxs), "METIS_PartGraphForContact: cinfo->dtree");
  cinfo->nvtxs    = *nvtxs;

  /*---------------------------------------------------------------------
   * Compute the initial k-way partitioning 
   *---------------------------------------------------------------------*/
  mcvwgt = idxsmalloc(2*(*nvtxs), 0, "METIS_PartGraphForContact: mcvwgt");
  for (i=0; i<*nvtxs; i++) {
    mcvwgt[2*i+0] = 1;
    mcvwgt[2*i+1] = (sflag[i] == 0 ? 0 : 1);
  }

  adjwgt = idxmalloc(xadj[*nvtxs], "METIS_PartGraphForContact: adjwgt");
  for (i=0; i<*nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      adjwgt[j] = (sflag[i] && sflag[adjncy[j]] ? 5 : 1);
  }

  rubvec[0] = 1.03;
  rubvec[1] = 1.05;
  ncon = 2;
  mcnumflag = 0;
  wgtflag   = 1;

  METIS_mCPartGraphKway(nvtxs, &ncon, xadj, adjncy, mcvwgt, adjwgt, &wgtflag, &mcnumflag,
                        nparts, rubvec, options, edgecut, part);

  /* The following is just for stat reporting purposes */
  SetUpGraph(&graph, OP_KMETIS, *nvtxs, 2, xadj, adjncy, mcvwgt, NULL, 0);
  graph.vwgt = mcvwgt;
  ComputePartitionBalance(&graph, *nparts, part, lbvec);
  mprintf("  %D-way Edge-Cut: %7D, Balance: %5.2f %5.2f\n", *nparts, ComputeCut(&graph, part), lbvec[0], lbvec[1]);


  /*---------------------------------------------------------------------
   * Induce the decission tree
   *---------------------------------------------------------------------*/
  dtpart = idxmalloc(*nvtxs, "METIS_PartGraphForContact: dtpart");
  marker = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: marker");

  for (dim=0; dim<3; dim++) {
    xyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*(*nvtxs), "METIS_PartGraphForContact: xyzcand[dim]");
    for (i=0; i<*nvtxs; i++) {
      xyzcand[dim][i].key = xyzcoords[3*i+dim];
      xyzcand[dim][i].val = i;
    }
    idkeysort(*nvtxs, xyzcand[dim]);
  }


  nnodes = nlnodes = nclean = naclean = ndirty = maxdepth = 0;
  InduceDecissionTree(*nvtxs, xyzcand, sflag, *nparts, part,
                      *nvtxs/(20*(*nparts)), *nvtxs/(20*(*nparts)*(*nparts)), 0.90,
                      &nnodes, &nlnodes, cinfo->dtree, leafpart, dtpart,
                      &nclean, &naclean, &ndirty, &maxdepth, marker);

  mprintf("NNodes: %5D, NLNodes: %5D, NClean: %5D, NAClean: %5D, NDirty: %5D, MaxDepth: %3D\n", nnodes, nlnodes, nclean, naclean, ndirty, maxdepth);


  /*---------------------------------------------------------------------
   * Create the tree-induced coarse graph and refine it
   *---------------------------------------------------------------------*/
  cgraph = CreatePartitionGraphForContact(*nvtxs, xadj, adjncy, mcvwgt, adjwgt, nlnodes, leafpart);

  for (i=0; i<*nvtxs; i++)
    part[leafpart[i]] = dtpart[i];

  ComputePartitionBalance(cgraph, *nparts, part, lbvec);
  mprintf("  %D-way Edge-Cut: %7D, Balance: %5.2f %5.2f\n", *nparts, ComputeCut(cgraph, part), lbvec[0], lbvec[1]);


  rwgtflag = 3;
  rnumflag = 0;
  METIS_mCRefineGraphKway(&(cgraph->nvtxs), &ncon, cgraph->xadj, cgraph->adjncy, cgraph->vwgt, 
                          cgraph->adjwgt, &rwgtflag, &rnumflag, nparts, rubvec, options, edgecut, 
                          part);

  ComputePartitionBalance(cgraph, *nparts, part, lbvec);
  mprintf("  %D-way Edge-Cut: %7D, Balance: %5.2f %5.2f\n", *nparts, ComputeCut(cgraph, part), lbvec[0], lbvec[1]);


  /*---------------------------------------------------------------------
   * Use that to compute the partition of the original graph
   *---------------------------------------------------------------------*/
  idxcopy(cgraph->nvtxs, part, dtpart);
  for (i=0; i<*nvtxs; i++)
    part[i] = dtpart[leafpart[i]];

  ComputePartitionBalance(&graph, *nparts, part, lbvec);
  idxset(*nvtxs, 1, graph.vwgt);
  mprintf("  %D-way Edge-Cut: %7D, Volume: %7D, Balance: %5.2f %5.2f\n", *nparts, 
           ComputeCut(&graph, part), ComputeVolume(&graph, part), lbvec[0], lbvec[1]);


  /*---------------------------------------------------------------------
   * Induce the final decission tree
   *---------------------------------------------------------------------*/
  nnodes = nlnodes = nclean = naclean = ndirty = maxdepth = 0;
  InduceDecissionTree(*nvtxs, xyzcand, sflag, *nparts, part,
                      *nvtxs/((40)*(*nparts)), 1, 1.00,
                      &nnodes, &nlnodes, cinfo->dtree, leafpart, dtpart, 
                      &nclean, &naclean, &ndirty, &maxdepth, marker);

  mprintf("NNodes: %5D, NLNodes: %5D, NClean: %5D, NAClean: %5D, NDirty: %5D, MaxDepth: %3D\n", nnodes, nlnodes, nclean, naclean, ndirty, maxdepth);

  
  /*---------------------------------------------------------------------
   * Populate the remaining fields of the cinfo data structure
   *---------------------------------------------------------------------*/
  cinfo->nnodes = nnodes;
  cinfo->nleafs = nlnodes;
  idxcopy(*nvtxs, part, cinfo->part);

  BuildDTLeafContents(cinfo, sflag);

  CheckDTree(*nvtxs, xyzcoords, part, cinfo);

  gk_free((void **)&mcvwgt, &dtpart, &xyzcand[0], &xyzcand[1], &xyzcand[2], &marker, &adjwgt, LTERM);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);

  return (void *)cinfo;
}


/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void METIS_UpdateContactInfo(void *raw_cinfo, idxtype *nvtxs, double *xyzcoords, idxtype *sflag)
{
  idxtype i, root, nchanges;
  ContactInfoType *cinfo;
  DTreeNodeType *dtree;

  cinfo = (ContactInfoType *)raw_cinfo;
  dtree = cinfo->dtree;

  if (cinfo->nvtxs != *nvtxs)
    errexit("The provided number of vertices do not match the initial information: %d %d\n", *nvtxs, cinfo->nvtxs);

  /* Go and reset the nsvtxs fields of all the nodes */
  for (i=0; i<cinfo->nnodes; i++) {
    dtree[i].nvtxs = 0;
    dtree[i].nsvtxs = 0;
  }

  /* Go and traverse each surface node and see where it gets assigned (We only focus on surface nodes. Right???) */
  for (nchanges=0, i=0; i<*nvtxs; i++) {
    if (1 || sflag[i]) { 
      for (root=0; dtree[root].leafid == -1;) 
        root = (xyzcoords[3*i+dtree[root].dim] <= dtree[root].value ? dtree[root].left : dtree[root].right);

      if (cinfo->leafpart[i] != dtree[root].leafid && sflag[i])
        nchanges++;

      cinfo->leafpart[i] = dtree[root].leafid;
      dtree[root].nvtxs++;
      if (sflag[i])
        dtree[root].nsvtxs++;
    }
  }

  mprintf("NChanges: %D\n", nchanges);

  BuildDTLeafContents(cinfo, sflag);

return;
for (i=0; i<cinfo->nnodes; i++) {
  if (dtree[i].leafid != -1)
    mprintf("%4D %4D %4D %4D\n", dtree[i].nvtxs, dtree[i].nsvtxs, dtree[i].leafid, dtree[i].partid);
}
}


/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void *METIS_SetupContact0(idxtype *nvtxs, double *xyzcoords, idxtype *sflag, 
                idxtype *nparts, idxtype *part) 
{
  idxtype i, j, ii, dim, ncontacts, wgtflag, mcnumflag, nnodes, nlnodes, nclean, naclean, ndirty, maxdepth, rwgtflag, rnumflag;
  idxtype *mcvwgt, *dtpart, *marker, *leafpart, *csflag, *cpart;
  idxtype *adjwgt;
  GraphType graph, *cgraph;
  ContactInfoType *cinfo;
  DKeyValueType *xyzcand[3];


  /*---------------------------------------------------------------------
   * Allocate memory for the contact info type
   *---------------------------------------------------------------------*/
  cinfo = (ContactInfoType *)gk_malloc(sizeof(ContactInfoType), "METIS_PartGraphForContact: cinfo");
  cinfo->leafptr  = idxsmalloc(*nvtxs+1, 0, "METIS_PartGraphForContact: leafptr");
  cinfo->leafind  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafind");
  cinfo->leafwgt  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafwgt");
  cinfo->part     = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: part");
  leafpart = cinfo->leafpart = idxmalloc(*nvtxs, "METIS_PartGraphForContact: leafpart");
  cinfo->dtree    = (DTreeNodeType *)gk_malloc(sizeof(DTreeNodeType)*(*nvtxs), "METIS_PartGraphForContact: cinfo->dtree");
  cinfo->nvtxs    = *nvtxs;


  /*---------------------------------------------------------------------
   * Induce the decission tree
   *---------------------------------------------------------------------*/
  dtpart = idxmalloc(*nvtxs, "METIS_PartGraphForContact: dtpart");
  marker = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: marker");

  for (dim=0; dim<3; dim++) 
    xyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*(*nvtxs), "METIS_PartGraphForContact: xyzcand[dim]");
  for (ncontacts=0, i=0; i<*nvtxs; i++) {
    if (sflag[i]) {
      for (dim=0; dim<3; dim++) {
        xyzcand[dim][ncontacts].key = xyzcoords[3*i+dim];
        xyzcand[dim][ncontacts].val = i;
      }
      ncontacts++;
    }
  }
  for (dim=0; dim<3; dim++) 
    idkeysort(ncontacts, xyzcand[dim]);


  nnodes = nlnodes = nclean = naclean = ndirty = maxdepth = 0;
  InduceDecissionTree(ncontacts, xyzcand, sflag, *nparts, part, ncontacts, 1, 1.00, 
                      &nnodes, &nlnodes, cinfo->dtree, leafpart, dtpart, &nclean, 
                      &naclean, &ndirty, &maxdepth, marker);

  mprintf("NNodes: %5D, NLNodes: %5D, NClean: %5D, NAClean: %5D, NDirty: %5D, MaxDepth: %3D\n", nnodes, nlnodes, nclean, naclean, ndirty, maxdepth);


  /*---------------------------------------------------------------------
   * Populate the remaining fields of the cinfo data structure
   *---------------------------------------------------------------------*/
  cinfo->nnodes = nnodes;
  cinfo->nleafs = nlnodes;
  idxcopy(*nvtxs, part, cinfo->part);

  BuildDTLeafContents(cinfo, sflag);

  CheckDTreeSurface(*nvtxs, xyzcoords, part, cinfo, sflag);

  gk_free((void **)&dtpart, &xyzcand[0], &xyzcand[1], &xyzcand[2], &marker, LTERM);
  
/*
for (i=0; i<nnodes; i++)
  mprintf("%4D %2D %4D\n", i, cinfo->dtree[i].leafid, cinfo->dtree[i].nsvtxs);
*/

  return (void *)cinfo;
}


/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void *METIS_SetupContact(idxtype *nvtxs, double *xyzcoords, idxtype *sflag, 
                idxtype *nparts, idxtype *part) 
{
  idxtype i, j, ii, dim, ncontacts, wgtflag, mcnumflag, nnodes, nlnodes, nclean, naclean, ndirty, maxdepth, rwgtflag, rnumflag;
  idxtype *mcvwgt, *dtpart, *marker, *leafpart, *csflag, *cpart;
  idxtype *adjwgt;
  GraphType graph, *cgraph;
  ContactInfoType *cinfo;
  DKeyValueType *xyzcand[3];


  /*---------------------------------------------------------------------
   * Allocate memory for the contact info type
   *---------------------------------------------------------------------*/
  cinfo = (ContactInfoType *)gk_malloc(sizeof(ContactInfoType), "METIS_PartGraphForContact: cinfo");
  cinfo->leafptr  = idxsmalloc(*nvtxs+1, 0, "METIS_PartGraphForContact: leafptr");
  cinfo->leafind  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafind");
  cinfo->leafwgt  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafwgt");
  cinfo->part     = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: part");
  leafpart = cinfo->leafpart = idxmalloc(*nvtxs, "METIS_PartGraphForContact: leafpart");
  cinfo->dtree    = (DTreeNodeType *)gk_malloc(sizeof(DTreeNodeType)*(*nvtxs), "METIS_PartGraphForContact: cinfo->dtree");
  cinfo->nvtxs    = *nvtxs;


  /*---------------------------------------------------------------------
   * Induce the decission tree
   *---------------------------------------------------------------------*/
  dtpart = idxmalloc(*nvtxs, "METIS_PartGraphForContact: dtpart");
  marker = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: marker");

  for (dim=0; dim<3; dim++) {
    xyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*(*nvtxs), "METIS_PartGraphForContact: xyzcand[dim]");
    for (i=0; i<*nvtxs; i++) {
      xyzcand[dim][i].key = xyzcoords[3*i+dim];
      xyzcand[dim][i].val = i;
    }
    idkeysort(*nvtxs, xyzcand[dim]);
  }


  nnodes = nlnodes = nclean = naclean = ndirty = maxdepth = 0;
  InduceDecissionTree(*nvtxs, xyzcand, sflag, *nparts, part,
                      *nvtxs, 1, 1.00,
                      &nnodes, &nlnodes, cinfo->dtree, leafpart, dtpart, 
                      &nclean, &naclean, &ndirty, &maxdepth, marker);

  mprintf("NNodes: %5D, NLNodes: %5D, NClean: %5D, NAClean: %5D, NDirty: %5D, MaxDepth: %3D\n", nnodes, nlnodes, nclean, naclean, ndirty, maxdepth);


  /*---------------------------------------------------------------------
   * Populate the remaining fields of the cinfo data structure
   *---------------------------------------------------------------------*/
  cinfo->nnodes = nnodes;
  cinfo->nleafs = nlnodes;
  idxcopy(*nvtxs, part, cinfo->part);

  BuildDTLeafContents(cinfo, sflag);

  CheckDTree(*nvtxs, xyzcoords, part, cinfo);

  gk_free((void **)&dtpart, &xyzcand[0], &xyzcand[1], &xyzcand[2], &marker, LTERM);
  
/*
for (i=0; i<nnodes; i++)
  mprintf("%4D %2D %4D\n", i, cinfo->dtree[i].leafid, cinfo->dtree[i].nsvtxs);
*/

  return (void *)cinfo;
}


/*************************************************************************
* This function is the entry point for detecting contacts between 
* bounding boxes and surface nodes
**************************************************************************/
void METIS_FindContacts(void *raw_cinfo, idxtype *nboxes, double *boxcoords, idxtype *nparts, 
               idxtype **r_cntptr, idxtype **r_cntind)
{
  idxtype i, ncnts, tncnts, maxtncnts;
  idxtype *cntptr, *cntind, *auxcntind, *stack, *marker;
  ContactInfoType *cinfo;

  cinfo = (ContactInfoType *)raw_cinfo;

  maxtncnts = 6*(*nboxes);
  cntptr    = idxsmalloc(*nboxes+1, 0, "METIS_FindContacts: cntptr");
  cntind    = idxmalloc(maxtncnts, "METIS_FindContacts: cntind");
  auxcntind = idxmalloc(*nparts, "METIS_FindContacts: auxcntind");
  stack     = idxmalloc(cinfo->nnodes, "METIS_FindContacts: stack");
  marker    = idxsmalloc(*nparts, 0, "METIS_FindContacts: marker");
  

  /* Go through each box and determine its contacting partitions */
  for (tncnts=0, i=0; i<*nboxes; i++) {
    ncnts = FindBoxContacts(cinfo, boxcoords+i*6, stack, auxcntind, marker);

    if (ncnts == 0)
      mprintf("CSearchError: Box has no contacts!\n");
  
    if (ncnts + tncnts >= maxtncnts) {
      maxtncnts += (tncnts+ncnts)*(*nboxes-i)/i;
      if ((cntind = (idxtype *)realloc(cntind, maxtncnts*sizeof(idxtype))) == NULL)
        errexit("Realloc failed! of %d words!\n", maxtncnts);
    }
    cntptr[i] = ncnts;
    idxcopy(ncnts, auxcntind, cntind+tncnts);
    tncnts += ncnts;
  }
  MAKECSR(i, *nboxes, cntptr); 

  *r_cntptr = cntptr;
  *r_cntind = cntind;

  gk_free((void **)&auxcntind, &stack, &marker, LTERM);

}


/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void METIS_FreeContactInfo(void *raw_cinfo)
{
  ContactInfoType *cinfo;

  cinfo = (ContactInfoType *)raw_cinfo;

  gk_free((void **)&(cinfo->leafptr), &(cinfo->leafind), &(cinfo->leafwgt), &(cinfo->part), &(cinfo->leafpart), &(cinfo->dtree), &cinfo, LTERM);
}



/******************************************************************************
* This function creates a coarse graph corresponding to the partitioning vector
*******************************************************************************/
GraphType *CreatePartitionGraphForContact(idxtype nvtxs, idxtype *xadj, idxtype *adjncy, 
                idxtype *vwgt, idxtype *adjwgt, idxtype cnvtxs, idxtype *part)
{
  idxtype i, ii, j, jj, k, cnedges;
  idxtype *cxadj, *cadjncy, *cvwgt, *cadjwgt;
  idxtype *ptr, *ind, *marker;
  GraphType *cgraph;

  ptr    = idxsmalloc(cnvtxs+1, 0, "CreatePartitionGraph: ptr");
  ind    = idxmalloc(nvtxs, "CreatePartitionGraph: ind");
  marker = idxsmalloc(cnvtxs, -1, "CreatePartitionGraph: marker");

  cgraph = CreateGraph();

  cgraph->ncon  = 2;
  cgraph->nvtxs = cnvtxs;
  cxadj   = cgraph->xadj   = idxsmalloc(cnvtxs+1, 0, "CreatePartitionGraph: cxadj");
  cadjncy = cgraph->adjncy = idxmalloc(xadj[nvtxs], "CreatePartitionGraph: cadjncy");
  cvwgt   = cgraph->vwgt   = idxmalloc(2*cnvtxs, "CreatePartitionGraph: cvwgt");
  cadjwgt = cgraph->adjwgt = idxmalloc(xadj[nvtxs], "CreatePartitionGraph: cadjwgt");


  for (i=0; i<nvtxs; i++) 
    ptr[part[i]]++;
  MAKECSR(i, cnvtxs, ptr);
  for (i=0; i<nvtxs; i++) 
    ind[ptr[part[i]]++] = i;
  SHIFTCSR(i, cnvtxs, ptr);

  
  /* Create the adjacency list of the coarse/partition graph */
  for (cxadj[0]=cnedges=0, i=0; i<cnvtxs; i++) {
    cvwgt[2*i+0] = cvwgt[2*i+1] = 0;
    for (j=ptr[i]; j<ptr[i+1]; j++) {
      ii = ind[j];
      cvwgt[2*i+0] += vwgt[2*ii+0];
      cvwgt[2*i+1] += vwgt[2*ii+1];

      for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
        if ((k = part[adjncy[jj]]) == i)
          continue;
        if (marker[k] == -1) {
          cadjncy[cnedges] = k;
          cadjwgt[cnedges] = adjwgt[jj];
          marker[k] = cnedges++;
        }
        else 
          cadjwgt[marker[k]] += adjwgt[jj];
      }
    }
    cxadj[i+1] = cnedges;

    /* reset the marker array */
    for (j=cxadj[i]; j<cnedges; j++)
      marker[cadjncy[j]] = -1;
  }

  gk_free((void **)&ptr, &ind, &marker, LTERM);

  return cgraph;
}

    

/*************************************************************************
* This function induces a DT that satisfied the given size requirements
**************************************************************************/
idxtype InduceDecissionTree(idxtype nvtxs, DKeyValueType **xyzcand, idxtype *sflag, 
          idxtype nparts, idxtype *part, idxtype maxnvtxs, idxtype minnvtxs, float minfrac, idxtype *r_nnodes, 
          idxtype *r_nlnodes, DTreeNodeType *dtree, idxtype *dtpart, idxtype *dtipart, 
          idxtype *r_nclean, idxtype *r_naclean, idxtype *r_ndirty, idxtype *r_maxdepth, idxtype *marker)
{
  idxtype i, nr, nl, j, k, mynodeID, dim, pid, bestdim, nleft, nright, 
      lmaxdepth, rmaxdepth, isclean, isleaf, cleanmarked, nsvtxs;
  idxtype *tpwgts, *pwgts[2], lnvtxs[3];
  double points[3], scores[3], sum2[2], newscore, bestscore, bestpoint;
  DKeyValueType *lxyzcand[3], *rxyzcand[3];

  *r_maxdepth = 1;

  mynodeID = (*r_nnodes)++;
  dtree[mynodeID].nvtxs  = nvtxs;
  dtree[mynodeID].nsvtxs = 0;
  dtree[mynodeID].leafid = -1;


  /* Determine overall pwgts */
  tpwgts = idxsmalloc(nparts, 0, "InduceDecissionTree: tpwgts");

  for (i=0; i<nvtxs; i++) {
    dtree[mynodeID].nsvtxs += (sflag[xyzcand[0][i].val] ? 1 : 0);
    tpwgts[part[xyzcand[0][i].val]]++;
  }

  pid = idxargmax(nparts, tpwgts);

  /*-----------------------------------------------------------------------
   * Check the exit conditions
   *-----------------------------------------------------------------------*/
  isclean = (tpwgts[pid] == nvtxs);
  isleaf  = 0;

  if (nvtxs <= maxnvtxs && isclean) { /* Determine if single class */
    //mprintf("LEAF1: %5D             Pure node!\n", nvtxs);
    *r_nclean += nvtxs;
    isleaf = 1;
  }
  else if (nvtxs < minnvtxs) { /* Determine if too small to continue */
    for (k=0, i=0; i<nparts; i++)
      k += (tpwgts[i] > 0 ? 1 : 0);
    //mprintf("LEAF3: %5D %5D       Skipping small node!\n", nvtxs, k);
    
    *r_ndirty += nvtxs*k;
    isleaf = 1;
  } else if (nvtxs < maxnvtxs && tpwgts[pid] >= (int)(minfrac*nvtxs)) { /* Determine if mostly one class */
    //mprintf("LEAF2: %5D %5D %4D  Almost pure node!\n", nvtxs, tpwgts[idxargmax(nparts, tpwgts)], idxargmax(nparts, tpwgts));
    *r_naclean += nvtxs;
    isleaf = 1;
  } else  { /* Check if all coordinates are the same */
    for (dim=0; dim<3; dim++) {
      for (i=1; i<nvtxs; i++) {
        if (fabs(xyzcand[dim][0].key - xyzcand[dim][i].key) > DEPSILON)
          break;
      }
      if (i != nvtxs)
        break;
    }

    if (dim == 3) { /* All coordinates matched! Treat it as a dirty node! */
      for (k=0, i=0; i<nparts; i++)
        k += (tpwgts[i] > 0 ? 1 : 0);
      mprintf("LEAF4: %5D %5D       Skipping same coord-nodes! (%D %D)\n", nvtxs, k, isclean, part[xyzcand[0][0].val]);
    
      *r_ndirty += nvtxs*k;
      isleaf = 1;
    }
  }

  if (isleaf) {
    for (i=0; i<nvtxs; i++) {
      dtpart[xyzcand[0][i].val]  = *r_nlnodes;
      dtipart[xyzcand[0][i].val] = pid;
    }

    dtree[mynodeID].leafid = (*r_nlnodes)++;
    dtree[mynodeID].partid = pid;
    dtree[mynodeID].left = dtree[mynodeID].right = -1;

    gk_free((void **)&tpwgts, LTERM);
    return mynodeID;
  }



  /*-----------------------------------------------------------------------
   * Find the best splitting point
   *-----------------------------------------------------------------------*/
  pwgts[0] = idxmalloc(nparts, "InduceDecissionTree: pwgts[0]");
  pwgts[1] = idxmalloc(nparts, "InduceDecissionTree: pwgts[1]");

  /* Go and scan each dimension */
  for (dim=0; dim<3; dim++) {
    /* Establish initial conditions for the scan */
    idxcopy(nparts, tpwgts, pwgts[1]);
    idxset(nparts, 0, pwgts[0]);

    sum2[0] = sum2[1] = 0.0;
    for (j=0; j<nparts; j++) 
      sum2[1] += (double)pwgts[1][j]*pwgts[1][j];

    scores[dim] = -1.0;
    nleft = 0;

    /* Scan until you find a well-separated vertex */
    for (i=0; i<nvtxs-1; i++) {
      pid = part[xyzcand[dim][i].val];
      sum2[0] += 1.0 + (double)2*pwgts[0][pid];
      sum2[1] += 1.0 - (double)2*pwgts[1][pid];
      pwgts[0][pid]++;
      pwgts[1][pid]--;
      nleft++;

      if (fabs(xyzcand[dim][i].key < xyzcand[dim][i+1].key) > DEPSILON) {
        scores[dim] = sqrt(sum2[0])+sqrt(sum2[1]);
        points[dim] = (xyzcand[dim][i].key+xyzcand[dim][i+1].key)/2.0;
        lnvtxs[dim] = nleft;
        break;
      }
    }

#ifdef PRINTSTAT
    if (i == nvtxs-1)
      mprintf("DTree: Initial Scan Along dim %D failed!\n", dim);
#endif


    /* Continue with the rest */
    for (i++; i<nvtxs-1; i++) {
      pid = part[xyzcand[dim][i].val];
      sum2[0] += 1.0 + (double)2*pwgts[0][pid];
      sum2[1] += 1.0 - (double)2*pwgts[1][pid];
      pwgts[0][pid]++;
      pwgts[1][pid]--;
      nleft++;

      newscore = sqrt(sum2[0])+sqrt(sum2[1]);
      if (isclean) {
        if (i >= nvtxs/2) {
          if (fabs(xyzcand[dim][i].key - xyzcand[dim][i+1].key) < DEPSILON) 
            continue;

          scores[dim] = xyzcand[dim][nvtxs-1].key - xyzcand[dim][0].key; /* Use the axis span as the score */
          points[dim] = (xyzcand[dim][i].key+xyzcand[dim][i+1].key)/2.0;
          lnvtxs[dim] = nleft;
          break;
        }
      }
      else {
        if (newscore > scores[dim]) {
          if (fabs(xyzcand[dim][i].key - xyzcand[dim][i+1].key) < DEPSILON) 
            continue;

          scores[dim] = newscore;
          points[dim] = (xyzcand[dim][i].key+xyzcand[dim][i+1].key)/2.0;
          lnvtxs[dim] = nleft;
        }
      }

      //mprintf("%5D %f %f %f\n", nleft, newscore, sum2[0], sum2[1]);
    }

#ifdef PRINTSTAT
    /* Print some Stats */
    if (scores[dim] >= 0) {
      mprintf("Dim: %3D, Score: %f, Point: %f [%5D %5D] [%f %f]\n", dim, scores[dim], 
              points[dim], lnvtxs[dim], nvtxs-lnvtxs[dim], xyzcand[dim][0].key, xyzcand[dim][nvtxs-1].key);
      idxcopy(nparts, tpwgts, pwgts[1]);
      idxset(nparts, 0, pwgts[0]);
      for (i=0; i<lnvtxs[dim]; i++) {
        pid = part[xyzcand[dim][i].val];
        pwgts[0][pid]++;
        pwgts[1][pid]--;
      }
      for (j=0; j<nparts; j++) 
        if (pwgts[0][j]+pwgts[1][j] > 0)
          mprintf("%5D => %5D %5D\n", j, pwgts[0][j], pwgts[1][j]);
    }
#endif
  }

  /* Determine the best overall score */
  bestdim   = 0;
  bestscore = scores[0];
  bestpoint = points[0];
  for (dim=1; dim<3; dim++) {
    if (scores[dim] > bestscore) {
      bestscore = scores[dim];
      bestpoint = points[dim];
      bestdim   = dim;
    }
  }

  if (bestscore <= 0.0)
    errexit("Major Failure... Non-seperable! %4d nodes. Needs to be fixed!\n", nvtxs);

  
  dtree[mynodeID].dim   = bestdim;
  dtree[mynodeID].value = bestpoint;

  //mprintf("BestDim: %D!\n", bestdim);

  
  /*-----------------------------------------------------------------------
   * Ok, now go and do the split 
   *-----------------------------------------------------------------------*/
  nleft  = lnvtxs[bestdim];
  nright = nvtxs - nleft;

  for (dim=0; dim<3; dim++) {
    lxyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*nleft,  "InduceDecissionTree: lxyzcand[dim]");
    rxyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*nright, "InduceDecissionTree: rxyzcand[dim]");
  }

  /* Mark the left vertices */
  for (i=0; i<nleft; i++) 
    marker[xyzcand[bestdim][i].val] = 1;

  for (dim=0; dim<3; dim++) {
    for (nl=nr=0, i=0; i<nvtxs; i++) {
      if (marker[xyzcand[dim][i].val]) 
        lxyzcand[dim][nl++] = xyzcand[dim][i];
      else
        rxyzcand[dim][nr++] = xyzcand[dim][i];
    }
  }
    
  /* Reset the marking */
  for (i=0; i<nleft; i++) 
    marker[xyzcand[bestdim][i].val] = 0;

  gk_free((void **)&tpwgts, &pwgts[0], &pwgts[1], LTERM);


  dtree[mynodeID].left  = InduceDecissionTree(nleft, lxyzcand, sflag, nparts, part, maxnvtxs, 
                                minnvtxs, minfrac, r_nnodes, r_nlnodes, dtree, dtpart, dtipart, r_nclean, 
                                r_naclean, r_ndirty, &lmaxdepth, marker);
  dtree[mynodeID].right = InduceDecissionTree(nright, rxyzcand, sflag, nparts, part, maxnvtxs, 
                                minnvtxs, minfrac, r_nnodes, r_nlnodes, dtree, dtpart, dtipart, r_nclean, 
                                r_naclean, r_ndirty, &rmaxdepth, marker);

  *r_maxdepth += amax(lmaxdepth, rmaxdepth);
 
  gk_free((void **)&lxyzcand[0], &lxyzcand[1], &lxyzcand[2], &rxyzcand[0], &rxyzcand[1], &rxyzcand[2], LTERM);
  
  return mynodeID;
}


/***********************************************************************************
* This function determines the partitions of the vertices assigned to each leaf
************************************************************************************/
void BuildDTLeafContents(ContactInfoType *cinfo, idxtype *sflag)
{
  idxtype i, j, k, nvtxs, ncontacts, nleafs, nind, tnind, tcomm;
  idxtype *part, *leafptr, *leafind, *leafwgt, *leafpart;
  KeyValueType *cand;

  nvtxs    = cinfo->nvtxs;
  nleafs   = cinfo->nleafs;
  part     = cinfo->part;
  leafpart = cinfo->leafpart;
  leafptr  = cinfo->leafptr;
  leafind  = cinfo->leafind;
  leafwgt  = cinfo->leafwgt;

  cand = (KeyValueType *)gk_malloc(sizeof(KeyValueType)*nvtxs, "BuildDTLeafContents: cand");

  for (ncontacts=0, i=0; i<nvtxs; i++) {
    if (sflag[i]) {
      cand[ncontacts].key   = leafpart[i];
      cand[ncontacts++].val = part[i];
    }
  }
  ikeyvalsort(ncontacts, cand);

/*
  for (i=0; i<ncontacts; i++)
    mprintf("%4D %5D %5D\n", i, cand[i].key, cand[i].val);
*/

  idxset(nleafs, 0, leafptr);

  leafind[0] = cand[0].val;
  leafwgt[0] = 1;
  nind = tnind = 1;
  for (i=1; i<ncontacts; i++) {
    if (cand[i].key != cand[i-1].key) {
      leafptr[cand[i-1].key] = nind;
      leafind[tnind]         = cand[i].val;
      leafwgt[tnind++]       = 1;
      nind = 1;
    }
    else if (cand[i].val != cand[i-1].val) {
      leafind[tnind]   = cand[i].val;
      leafwgt[tnind++] = 1;
      nind++;
    }
    else {
      leafwgt[tnind-1]++;
    }
  }
  leafptr[cand[i-1].key] = nind;

  MAKECSR(i, nleafs, leafptr);

  for (tcomm=0; i<nleafs; i++) {
    tcomm += (leafptr[i+1]-leafptr[i]-1)*idxsum(leafptr[i+1]-leafptr[i], leafwgt+leafptr[i], 1);

/*
    if (leafptr[i+1]-leafptr[i] > 1) {
      mprintf("%4D, ", i);
      for (j=leafptr[i]; j<leafptr[i+1]; j++)
        mprintf("[%3D %4D] ", leafind[j], leafwgt[j]);
      mprintf("\n");
    }
*/

  }
      

  mprintf("NLeafs: %D, NLeafIndices: %D, EstimComm: %D\n", nleafs, leafptr[nleafs], tcomm);

/*
for (i=0; i<nleafs; i++) {
  mprintf("Leaf: %D => ", i);
  for (j=leafptr[i]; j<leafptr[i+1]; j++)
    mprintf("[%D %D] ", leafind[j], leafwgt[j]);
  mprintf("\n");
}
*/

  gk_free((void **)&cand, LTERM);
}


/***********************************************************************************
* This function determines the leaf-nodes that are intersected by a bbox
************************************************************************************/
idxtype FindBoxContacts(ContactInfoType *cinfo, double *coords, idxtype *stack, idxtype *cntind, idxtype *marker)
{
  idxtype i, j, l, k, ncnts, root, leafid;
  idxtype *leafptr, *leafind;
  DTreeNodeType *dtree;

  leafptr = cinfo->leafptr;
  leafind = cinfo->leafind;
  dtree   = cinfo->dtree;

//mprintf("%e %e %e %e %e %e\n", coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);

  stack[0] = 0;
  k = 1;
  ncnts = 0;
  while (k > 0) {
    root = stack[--k]; 

    /* See if you hit a leaf */
    if ((leafid = dtree[root].leafid) != -1) {
//mprintf("Got to a leaf....  %D %D %D\n", leafid, dtree[root].nsvtxs, leafptr[leafid+1]-leafptr[leafid]);

      if (dtree[root].nsvtxs > 0) {
        /* Add unique processor IDs. */
        for (j=leafptr[leafid]; j<leafptr[leafid+1]; j++) {
          if (!marker[leafind[j]]) { 
            cntind[ncnts++] = leafind[j];
            marker[leafind[j]] = 1;
          }
        }
//        mprintf("[%4D %4D %2D %2D] ", root, dtree[root].nsvtxs, dtree[root].partid, ncnts);
      }

      continue;
    }

//mprintf("Making decissions... %D %e %D %D = > ", dtree[root].dim, dtree[root].value, dtree[root].left, dtree[root].right);

    if (coords[dtree[root].dim] <= dtree[root].value) { /* min <= value */
      stack[k++] = dtree[root].left;
//mprintf(" %D", stack[k-1]);
    }
    
    if (coords[3+dtree[root].dim] >= dtree[root].value) { /* max >= value */
      stack[k++] = dtree[root].right;
//mprintf(" %D", stack[k-1]);
    }

//mprintf("\n");

  }

  for (i=0; i<ncnts; i++)
    marker[cntind[i]] = 0;

//mprintf("\n");

  return ncnts;

}



/***********************************************************************************
* This function checks the DT to see if it properly "classifies" all points
************************************************************************************/
void CheckDTree(idxtype nvtxs, double *xyzcoords, idxtype *part, ContactInfoType *cinfo)
{
  idxtype i, j, k, root;
  idxtype *leafptr, *leafind;
  DTreeNodeType *dtree;

  leafptr = cinfo->leafptr;
  leafind = cinfo->leafind;
  dtree   = cinfo->dtree;

  for (i=0; i<nvtxs; i++) {
    for (root = 0; dtree[root].leafid == -1;) 
      root = (xyzcoords[3*i+dtree[root].dim] <= dtree[root].value ? dtree[root].left : dtree[root].right);

    if (cinfo->leafpart[i] != dtree[root].leafid)
      mprintf("DTError! %4D %4D %4D %4D %4D\n", i, cinfo->leafpart[i], dtree[root].leafid, part[i], leafind[leafptr[dtree[root].leafid]]);
  }

}


/***********************************************************************************
* This function checks the DT to see if it properly "classifies" all points
************************************************************************************/
void CheckDTreeSurface(idxtype nvtxs, double *xyzcoords, idxtype *part, ContactInfoType *cinfo, idxtype *sflag)
{
  idxtype i, j, k, root;
  idxtype *leafptr, *leafind;
  DTreeNodeType *dtree;

  leafptr = cinfo->leafptr;
  leafind = cinfo->leafind;
  dtree   = cinfo->dtree;

  for (i=0; i<nvtxs; i++) {
    if (!sflag[i])
      continue;

    for (root = 0; dtree[root].leafid == -1;) 
      root = (xyzcoords[3*i+dtree[root].dim] <= dtree[root].value ? dtree[root].left : dtree[root].right);

    if (cinfo->leafpart[i] != dtree[root].leafid)
      mprintf("SDTError! %4D %4D %4D %4D %4D\n", i, cinfo->leafpart[i], dtree[root].leafid, part[i], leafind[leafptr[dtree[root].leafid]]);
  }

}







/*************************************************************************
* This function is the entry point for KMETIS
**************************************************************************/
void *METIS_PartSurfForContactRCB(idxtype *nvtxs, double *xyzcoords, idxtype *sflag, 
                idxtype *nparts, idxtype *part, idxtype *bestdims) 
{
  idxtype i, j, nsurf, dim, ncon, nnodes, nlnodes;
  idxtype *marker, *spart;
  ContactInfoType *cinfo;
  DKeyValueType *xyzcand[3];
  double *myxyzcoords;


  /*---------------------------------------------------------------------
   * Allocate memory for the contact info type
   *---------------------------------------------------------------------*/
  cinfo = (ContactInfoType *)gk_malloc(sizeof(ContactInfoType), "METIS_PartGraphForContact: cinfo");
  cinfo->leafptr  = idxsmalloc(*nvtxs+1, 0, "METIS_PartGraphForContact: leafptr");
  cinfo->leafind  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafind");
  cinfo->leafwgt  = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: leafwgt");
  cinfo->part     = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: part");
  cinfo->leafpart = idxmalloc(*nvtxs, "METIS_PartGraphForContact: leafpart");
  cinfo->dtree    = (DTreeNodeType *)gk_malloc(sizeof(DTreeNodeType)*(*nvtxs), "METIS_PartGraphForContact: cinfo->dtree");

  /*---------------------------------------------------------------------
   * Induce the decission tree
   *---------------------------------------------------------------------*/
  myxyzcoords = gk_dmalloc(3*(*nvtxs), "METIS_PartSurfForContactRCB: myxyzcoords");
  marker = idxsmalloc(*nvtxs, 0, "METIS_PartGraphForContact: marker");

  for (dim=0; dim<3; dim++) {
    xyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*(*nvtxs), "METIS_PartGraphForContact: xyzcand[dim]");

    for (nsurf=0, i=0; i<*nvtxs; i++) {
      if (sflag[i]) {
        myxyzcoords[3*nsurf+dim] = xyzcoords[3*i+dim];
        xyzcand[dim][nsurf].key  = xyzcoords[3*i+dim];
        xyzcand[dim][nsurf].val  = nsurf++;
      }
    }
    idkeysort(nsurf, xyzcand[dim]);
  }

  spart = idxsmalloc(nsurf, 0, "METIS_PartGraphForContact: spart");

  nnodes = nlnodes = 0;
  InduceRCBTree(nsurf, xyzcand, 0, *nparts, &nnodes, &nlnodes, cinfo->dtree, cinfo->leafpart, spart, marker, bestdims);

  mprintf("NNodes: %5D, NLNodes: %5D\n", nnodes, nlnodes);

  /* Project the partition back to the original space */
  for (nsurf=0, i=0; i<*nvtxs; i++) 
    part[i] = (sflag[i] ? spart[nsurf++] : -1);


  /*---------------------------------------------------------------------
   * Populate the remaining fields of the cinfo data structure
   *---------------------------------------------------------------------*/
  cinfo->nvtxs  = nsurf;
  cinfo->nnodes = nnodes;
  cinfo->nleafs = nlnodes;
  idxcopy(nsurf, spart, cinfo->part);

  idxset(nsurf, 1, marker);

  BuildDTLeafContents(cinfo, marker);

  CheckDTree(nsurf, myxyzcoords, spart, cinfo);

  gk_free((void **)&xyzcand[0], &xyzcand[1], &xyzcand[2], &myxyzcoords, &marker, &spart, LTERM);
  
  for (i=0; i<cinfo->nnodes; i++)
    bestdims[i] = cinfo->dtree[i].dim;

  return (void *)cinfo;
}



/*************************************************************************
* This function induces a DT that satisfied the given size requirements
**************************************************************************/
idxtype InduceRCBTree(idxtype nvtxs, DKeyValueType **xyzcand, idxtype firstPID, idxtype nparts, 
          idxtype *r_nnodes, idxtype *r_nlnodes, DTreeNodeType *dtree, idxtype *leafpart,
          idxtype *part, idxtype *marker, idxtype *oldBestDims)
{
  idxtype i, nr, nl, j, k, mynodeID, dim, bestdim, lnvtxs, rnvtxs, lnparts, rnparts;
  DKeyValueType *lxyzcand[3], *rxyzcand[3];
  double bestpoint;


  mynodeID = (*r_nnodes)++;
  dtree[mynodeID].nvtxs  = nvtxs;
  dtree[mynodeID].nsvtxs = nvtxs;
  dtree[mynodeID].leafid = -1;


  /*-----------------------------------------------------------------------
   * Check the exit conditions
   *-----------------------------------------------------------------------*/
  if (nparts == 1) {
    //mprintf("Pid:%D, Size:%D\n", firstPID, nvtxs);
    for (i=0; i<nvtxs; i++) {
      leafpart[xyzcand[0][i].val] = *r_nlnodes;
      part[xyzcand[0][i].val]     = firstPID;
    }

    dtree[mynodeID].leafid = (*r_nlnodes)++;
    dtree[mynodeID].partid = firstPID;
    dtree[mynodeID].left = dtree[mynodeID].right = -1;

    return mynodeID;
  }


  /*-----------------------------------------------------------------------
   * Find the best splitting point
   *-----------------------------------------------------------------------*/
  lnparts = nparts/2;
  rnparts = nparts - lnparts;
  lnvtxs = nvtxs*lnparts/nparts;
  rnvtxs = -1;

  if ((bestdim = oldBestDims[mynodeID]) != -1) {
    if (fabs(xyzcand[bestdim][lnvtxs].key - xyzcand[bestdim][lnvtxs+1].key) < DEPSILON)
      lnvtxs = 0.99*lnvtxs;

    for (; lnvtxs<nvtxs; lnvtxs++) {
      if (fabs(xyzcand[bestdim][lnvtxs].key - xyzcand[bestdim][lnvtxs+1].key) > DEPSILON) 
        break;
    }
    lnvtxs++;
    rnvtxs = nvtxs - lnvtxs;
  }

  if (rnvtxs <= 0) {
    if (bestdim != -1)
      mprintf("Finding a dimension for %D points...\n", nvtxs);

    lnparts = nparts/2;
    rnparts = nparts - lnparts;
    lnvtxs = nvtxs*lnparts/nparts;

    for (bestdim=0, i=1; i<3; i++) 
      bestdim = (xyzcand[i][nvtxs-1].key - xyzcand[i][0].key > xyzcand[bestdim][nvtxs-1].key - xyzcand[bestdim][0].key ? i : bestdim);

    for (; lnvtxs<nvtxs; lnvtxs++) {
      if (fabs(xyzcand[bestdim][lnvtxs].key - xyzcand[bestdim][lnvtxs+1].key) > DEPSILON) 
        break;
    }
    lnvtxs++;
    rnvtxs = nvtxs - lnvtxs;
  }


  dtree[mynodeID].dim   = bestdim;
  dtree[mynodeID].value = bestpoint = (xyzcand[bestdim][lnvtxs-1].key+xyzcand[bestdim][lnvtxs].key)/2;

  /* Print some Stats */
  //mprintf("Dim: %3D, Point: %f [%5D %5D] [%3D %3D]\n", bestdim, bestpoint, lnvtxs, rnvtxs, lnparts, rnparts);

  /*-----------------------------------------------------------------------
   * Ok, now go and do the split 
   *-----------------------------------------------------------------------*/
  for (dim=0; dim<3; dim++) {
    lxyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*lnvtxs,  "InduceDecissionTree: lxyzcand[dim]");
    rxyzcand[dim] = (DKeyValueType *)gk_malloc(sizeof(DKeyValueType)*rnvtxs, "InduceDecissionTree: rxyzcand[dim]");
  }

  /* Mark the left vertices */
  for (i=0; i<lnvtxs; i++) 
    marker[xyzcand[bestdim][i].val] = 1;

  for (dim=0; dim<3; dim++) {
    for (nl=nr=0, i=0; i<nvtxs; i++) {
      if (marker[xyzcand[dim][i].val]) 
        lxyzcand[dim][nl++] = xyzcand[dim][i];
      else
        rxyzcand[dim][nr++] = xyzcand[dim][i];
    }
  }
    
  /* Reset the marking */
  for (i=0; i<lnvtxs; i++) 
    marker[xyzcand[bestdim][i].val] = 0;


  dtree[mynodeID].left  = InduceRCBTree(lnvtxs, lxyzcand, firstPID, lnparts, r_nnodes, r_nlnodes, dtree, 
                                leafpart, part, marker, oldBestDims);
  dtree[mynodeID].right = InduceRCBTree(rnvtxs, rxyzcand, firstPID+lnparts, rnparts, r_nnodes, r_nlnodes, 
                                dtree, leafpart, part, marker, oldBestDims);

  gk_free((void **)&lxyzcand[0], &lxyzcand[1], &lxyzcand[2], &rxyzcand[0], &rxyzcand[1], &rxyzcand[2], LTERM);
  
  return mynodeID;
}

