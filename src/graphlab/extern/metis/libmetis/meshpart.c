/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * meshpart.c
 *
 * This file contains routines for partitioning finite element meshes.
 *
 * Started 9/29/97
 * George
 *
 * $Id: meshpart.c,v 1.2 2002/08/10 06:29:31 karypis Exp $
 *
 */

#include <metislib.h>


/*************************************************************************
* This function partitions a finite element mesh by partitioning its nodal
* graph using KMETIS and then assigning elements in a load balanced fashion.
**************************************************************************/
void METIS_PartMeshNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, 
idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)

{
  idxtype i, j, k, me;
  idxtype *xadj, *adjncy, *pwgts;
  idxtype options[10], pnumflag=0, wgtflag=0;
  idxtype nnbrs, nbrind[200], nbrwgt[200], maxpwgt;
  idxtype esize, esizes[] = {-1, 3, 4, 8, 4, 2};

  esize = esizes[*etype];

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esize, elmnts);

  xadj = idxmalloc(*nn+1, "METIS_MESHPARTNODAL: xadj");
  adjncy = idxmalloc(20*(*nn), "METIS_MESHPARTNODAL: adjncy");

  METIS_MeshToNodal(ne, nn, elmnts, etype, &pnumflag, xadj, adjncy);

  adjncy = realloc(adjncy, xadj[*nn]*sizeof(idxtype));

  options[0] = 0;
  METIS_PartGraphKway(nn, xadj, adjncy, NULL, NULL, &wgtflag, &pnumflag, nparts, options, edgecut, npart);

  /* OK, now compute an element partition based on the nodal partition npart */
  idxset(*ne, -1, epart);
  pwgts = idxsmalloc(*nparts, 0, "METIS_MESHPARTNODAL: pwgts");
  for (i=0; i<*ne; i++) {
    me = npart[elmnts[i*esize]];
    for (j=1; j<esize; j++) {
      if (npart[elmnts[i*esize+j]] != me)
        break;
    }
    if (j == esize) {
      epart[i] = me;
      pwgts[me]++;
    }
  }

  maxpwgt = 1.03*(*ne)/(*nparts);
  for (i=0; i<*ne; i++) {
    if (epart[i] == -1) { /* Assign the boundary element */
      nnbrs = 0;
      for (j=0; j<esize; j++) {
        me = npart[elmnts[i*esize+j]];
        for (k=0; k<nnbrs; k++) {
          if (nbrind[k] == me) {
            nbrwgt[k]++;
            break;
          }
        }
        if (k == nnbrs) {
          nbrind[nnbrs] = me;
          nbrwgt[nnbrs++] = 1;
        }
      }
      /* Try to assign it first to the domain with most things in common */
      j = idxargmax(nnbrs, nbrwgt);
      if (pwgts[nbrind[j]] < maxpwgt) {
        epart[i] = nbrind[j];
      }
      else {
        /* If that fails, assign it to a light domain */
        for (j=0; j<nnbrs; j++) {
          if (pwgts[nbrind[j]] < maxpwgt) {
            epart[i] = nbrind[j];
            break;
          }
        }
        if (j == nnbrs) 
          epart[i] = nbrind[idxargmax(nnbrs, nbrwgt)];
      }
      pwgts[epart[i]]++;
    }
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering2((*ne)*esize, elmnts, *ne, *nn, epart, npart);

  gk_free((void **)&xadj, &adjncy, &pwgts, LTERM);

}




/*************************************************************************
* This function partitions a Mixed element mesh by partitioning its nodal
* graph using KMETIS and then assigning elements in a load balanced fashion.
**************************************************************************/
void METIS_PartMixedMeshNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, 
idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)

{
  idxtype i, j, k, me, tot;
  idxtype *xadj, *adjncy, *pwgts, *hash;
  idxtype options[10], pnumflag=0, wgtflag=0;
  idxtype nnbrs, nbrind[200], nbrwgt[200], maxpwgt;
  idxtype esize, sizes[] = {-1, 3, 4, 8, 4, 2};

  hash = idxsmalloc((*ne), 0, "METIS_MIXEDMESHPARTNODAL: hash");
   tot=0;
   for (i=0;i<(*ne);i++){
      hash[i]=tot;
      tot+=sizes[etype[i]];
      }

  if (*numflag == 1)
    ChangeMesh2CNumbering(tot,  elmnts);

  xadj = idxmalloc(*nn+1, "METIS_MIXEDMESHPARTNODAL: xadj");
  adjncy = idxmalloc(20*(*nn), "METIS_MIXEDMESHPARTNODAL: adjncy");

  METIS_MixedMeshToNodal(ne, nn, elmnts, etype, &pnumflag, xadj, adjncy);

  adjncy = realloc(adjncy, xadj[*nn]*sizeof(idxtype));

  options[0] = 0;
  METIS_PartGraphKway(nn, xadj, adjncy, NULL, NULL, &wgtflag, &pnumflag, nparts, options, edgecut, npart);

  /* OK, now compute an element partition based on the nodal partition npart */
  idxset(*ne, -1, epart);
  pwgts = idxsmalloc(*nparts, 0, "METIS_MIXEDMESHPARTNODAL: pwgts");
  for (i=0; i<*ne; i++) {
    me = npart[elmnts[hash[i]]];
    for (j=1; j<sizes[etype[i]]; j++) {
      if (npart[elmnts[hash[i]+j]] != me)
        break;
    }
    if (j == sizes[etype[i]]) {
      epart[i] = me;
      pwgts[me]++;
    }
  }

  maxpwgt = 1.03*(*ne)/(*nparts);
  for (i=0; i<*ne; i++) {
    if (epart[i] == -1) { /* Assign the boundary element */
      nnbrs = 0;
      for (j=0; j<sizes[etype[i]]; j++) {
        me = npart[elmnts[hash[i]+j]];
        for (k=0; k<nnbrs; k++) {
          if (nbrind[k] == me) {
            nbrwgt[k]++;
            break;
          }
        }
        if (k == nnbrs) {
          nbrind[nnbrs] = me;
          nbrwgt[nnbrs++] = 1;
        }
      }
      /* Try to assign it first to the domain with most things in common */
      j = idxargmax(nnbrs, nbrwgt);
      if (pwgts[nbrind[j]] < maxpwgt) {
        epart[i] = nbrind[j];
      }
      else {
        /* If that fails, assign it to a light domain */
        for (j=0; j<nnbrs; j++) {
          if (pwgts[nbrind[j]] < maxpwgt) {
  epart[i] = nbrind[j];
            break;
          }
        }
        if (j == nnbrs)
          epart[i] = nbrind[idxargmax(nnbrs, nbrwgt)];
      }
      pwgts[epart[i]]++;
    }
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering2(tot, elmnts, *ne, *nn, epart, npart);

  gk_free((void **)&xadj, &adjncy, &pwgts, &hash, LTERM);

}



/*************************************************************************
* This function partitions a finite element mesh by partitioning its dual
* graph using KMETIS and then assigning nodes in a load balanced fashion.
**************************************************************************/
void METIS_PartMeshDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, 
idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart,
idxtype wgtflag, idxtype * vwgt) 
{
  idxtype i, j, k, me;
  idxtype *xadj, *adjncy, *pwgts, *nptr, *nind, *elms, cnt;
  idxtype options[10], pnumflag=0;
  idxtype nnbrs, nbrind[200], nbrwgt[200], maxpwgt;
  idxtype esize, esizes[] = {-1, 3, 4, 8, 4, 2};

  esize = esizes[*etype];

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esize, elmnts);

  xadj = idxmalloc(*ne+1, "METIS_MESHPARTNODAL: xadj");
  elms = idxsmalloc(*ne+1, 0, "main: elms");


  cnt=METIS_MeshToDualCount(ne, nn, elmnts, elms, etype, &pnumflag);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MeshToDual(ne, nn, elmnts, elms, etype, &pnumflag, xadj, adjncy);


  options[0] = 0;
  METIS_PartGraphKway(ne, xadj, adjncy, vwgt, NULL, &wgtflag, &pnumflag, nparts, options, edgecut, epart);

  /* Construct the node-element list */
  nptr = idxsmalloc(*nn+1, 0, "METIS_MESHPARTDUAL: nptr");
  for (j=esize*(*ne), i=0; i<j; i++) 
    nptr[elmnts[i]]++;
  MAKECSR(i, *nn, nptr);

  nind = idxmalloc(nptr[*nn], "METIS_MESHPARTDUAL: nind");
  for (k=i=0; i<(*ne); i++) {
    for (j=0; j<esize; j++, k++) 
      nind[nptr[elmnts[k]]++] = i;
  }
  for (i=(*nn); i>0; i--)
    nptr[i] = nptr[i-1];
  nptr[0] = 0;


  /* OK, now compute a nodal partition based on the element partition npart */
  idxset(*nn, -1, npart);
  pwgts = idxsmalloc(*nparts, 0, "METIS_MESHPARTDUAL: pwgts");
  for (i=0; i<*nn; i++) {
    me = epart[nind[nptr[i]]];
    for (j=nptr[i]+1; j<nptr[i+1]; j++) {
      if (epart[nind[j]] != me)
        break;
    }
    if (j == nptr[i+1]) {
      npart[i] = me;
      pwgts[me]++;
    }
  }

  maxpwgt = 1.03*(*nn)/(*nparts);
  for (i=0; i<*nn; i++) {
    if (npart[i] == -1) { /* Assign the boundary element */
      nnbrs = 0;
      for (j=nptr[i]; j<nptr[i+1]; j++) {
        me = epart[nind[j]];
        for (k=0; k<nnbrs; k++) {
          if (nbrind[k] == me) {
            nbrwgt[k]++;
            break;
          }
        }
        if (k == nnbrs) {
          nbrind[nnbrs] = me;
          nbrwgt[nnbrs++] = 1;
        }
      }
      /* Try to assign it first to the domain with most things in common */
      j = idxargmax(nnbrs, nbrwgt);
      if (pwgts[nbrind[j]] < maxpwgt) {
        npart[i] = nbrind[j];
      }
      else {
        /* If that fails, assign it to a light domain */
        npart[i] = nbrind[0];
        for (j=0; j<nnbrs; j++) {
          if (pwgts[nbrind[j]] < maxpwgt) {
            npart[i] = nbrind[j];
            break;
          }
        }
      }
      pwgts[npart[i]]++;
    }
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering2((*ne)*esize, elmnts, *ne, *nn, epart, npart);

  gk_free((void **)&xadj, &adjncy, &pwgts, &nptr, &nind, LTERM);

}



/*************************************************************************
* This function partitions a finite element mesh by partitioning its dual
* graph using KMETIS and then assigning nodes in a load balanced fashion.
**************************************************************************/
void METIS_PartMixedMeshDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, 
idxtype *numflag,  idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart,
idxtype *conmat, idxtype custom, idxtype wgtflag, idxtype *vwgt)
{
  idxtype i, j, k, l, me, tot, cnt;
  idxtype *xadj, *adjncy, *pwgts, *nptr, *nind, *hash, *elms;
  idxtype options[10], pnumflag=0;
  idxtype nnbrs, nbrind[200], nbrwgt[200], maxpwgt;
  idxtype esize, sizes[] = {-1, 3, 4, 8, 4, 2};

  hash = idxsmalloc((*ne), 0, "METIS_MIXEDMESHPARTNODAL: hash");
   tot=0;
   for (i=0;i<(*ne);i++){
      hash[i]=tot;
      tot+=sizes[etype[i]];
      }
 

  if (*numflag == 1)
    ChangeMesh2CNumbering(tot, elmnts);

  xadj = idxmalloc(*ne+1, "METIS_MESHPARTNODAL: xadj");
  elms = idxsmalloc(*ne+1, 0, ": elms"); 

  if (custom==1){
  cnt=METIS_MixedMeshToDualCount(ne, nn, elmnts, elms, etype, &pnumflag, conmat, 1);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MixedMeshToDual(ne, nn, elmnts, elms, etype, &pnumflag, xadj, adjncy, conmat, 1);
  }
  else{
  cnt=METIS_MixedMeshToDualCount(ne, nn, elmnts, elms, etype, &pnumflag, conmat, 0);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MixedMeshToDual(ne, nn, elmnts, elms, etype, &pnumflag, xadj, adjncy, conmat, 0);
  }

  options[0] = 0;
  METIS_PartGraphKway(ne, xadj, adjncy, vwgt, NULL, &wgtflag, &pnumflag, nparts, options, edgecut, epart);

  /* Construct the node-element list */
  nptr = idxsmalloc(*nn+1, 0, "METIS_MESHPARTDUAL: nptr");
 
   for(i=0;i<(*ne);i++){
   l=hash[i];
   for (j=(sizes[etype[i]]), k=0; k<j; k++)
     nptr[elmnts[l+k]]++;

   }
 
  MAKECSR(i, *nn, nptr);

  nind = idxmalloc(nptr[*nn], "METIS_MESHPARTDUAL: nind");
  for (k=i=0; i<(*ne); i++) {
    for (j=0; j<sizes[etype[i]]; j++, k++)
      nind[nptr[elmnts[k]]++] = i;
  }
  for (i=(*nn); i>0; i--)
    nptr[i] = nptr[i-1];
  nptr[0] = 0;


  /* OK, now compute a nodal partition based on the element partition npart */
  idxset(*nn, -1, npart);
  pwgts = idxsmalloc(*nparts, 0, "METIS_MESHPARTDUAL: pwgts");
  for (i=0; i<*nn; i++) {
    me = epart[nind[nptr[i]]];
    for (j=nptr[i]+1; j<nptr[i+1]; j++) {
      if (epart[nind[j]] != me)
        break;
    }
    if (j == nptr[i+1]) {
      npart[i] = me;
      pwgts[me]++;
    }
  }

  maxpwgt = 1.03*(*nn)/(*nparts);
  for (i=0; i<*nn; i++) {
    if (npart[i] == -1) { /* Assign the boundary element */
      nnbrs = 0;
      for (j=nptr[i]; j<nptr[i+1]; j++) {
        me = epart[nind[j]];
        for (k=0; k<nnbrs; k++) {
          if (nbrind[k] == me) {
            nbrwgt[k]++;
            break;
          }
        }
        if (k == nnbrs) {
          nbrind[nnbrs] = me;
          nbrwgt[nnbrs++] = 1;
        }
      }
  /* Try to assign it first to the domain with most things in common */
      j = idxargmax(nnbrs, nbrwgt);
      if (pwgts[nbrind[j]] < maxpwgt) {
        npart[i] = nbrind[j];
      }
      else {
        /* If that fails, assign it to a light domain */
        npart[i] = nbrind[0];
        for (j=0; j<nnbrs; j++) {
          if (pwgts[nbrind[j]] < maxpwgt) {
            npart[i] = nbrind[j];
            break;
          }
        }
      }
      pwgts[npart[i]]++;
    }
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering2((*ne)*esize, elmnts, *ne, *nn, epart, npart);

  gk_free((void **)&xadj, &adjncy, &pwgts, &nptr, &nind, LTERM);

}


