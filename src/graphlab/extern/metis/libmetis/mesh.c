/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh.c
 *
 * This file contains routines for converting 3D and 4D finite element
 * meshes into dual or nodal graphs
 *
 * Started 8/18/97
 * George
 *
 * $Id: mesh.c,v 1.2 2002/08/10 06:29:31 karypis Exp $
 *
 */

#include <metislib.h>

/******************************************************************************
* This function counts the neighbours of each element
******************************************************************************/
idxtype METIS_MeshToDualCount(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *elms, 
idxtype *etype, idxtype *numflag)
{
  idxtype cnt, esizes[] = {-1, 3, 4, 8, 4, 2};

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esizes[*etype], elmnts);

  cnt=GENDUALMETIS_COUNT(*ne, *nn, *etype, elmnts, elms);

  if (*numflag == 1)
    ChangeMesh2FNumbering3((*ne)*esizes[*etype], elmnts);
  
  return cnt;
}



/*****************************************************************************
* This function creates a graph corresponding to the dual of a finite element
* mesh. At this point the supported elements are triangles, tetrahedrons, and
* bricks.
******************************************************************************/
void METIS_MeshToDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *elms, 
idxtype *etype, idxtype *numflag,  idxtype *dxadj, idxtype *dadjncy)
{
  idxtype esizes[] = {-1, 3, 4, 8, 4, 2};

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esizes[*etype], elmnts);

  GENDUALMETIS(*ne, *nn, *etype, elmnts, elms, dxadj, dadjncy);

  if (*numflag == 1)
    ChangeMesh2FNumbering((*ne)*esizes[*etype], elmnts, *ne, dxadj, dadjncy);
}





/******************************************************************************
* This function counts the neighbours of each element 
******************************************************************************/
idxtype METIS_MixedMeshToDualCount(idxtype *ne, idxtype *nn, idxtype *elmnts, 
idxtype * elms, idxtype *etype, idxtype *numflag, idxtype *conmat, idxtype custom)
{

   idxtype i, j, jj, k, kk, kkk, l, m, tot,  n, nedges, mask, cnt=0;
   idxtype *nptr, *nind;
   idxtype *mark, ind[200], wgt[200];
   idxtype sizes[] = {-1, 3, 4, 8, 4, 2},
       mgcnums[] = {-1, 2, 3, 4, 2};

   idxtype *hash;

   idxtype mgcnum[6][6] = {-1, -1, -1, -1, -1, -1,
                       -1, 2, 3, 3, 2 , 2,
                       -1, 3, 3, 3, 3, 2,
                       -1, 3, 3, 4, 4, 2,
                       -1, 2, 3, 4, 2, 2,
                       -1, 2, 2, 2, 2, 1} ;


   if (custom==1) /* External magic numbers supplied */
      for (i=0,k=0;i<6;i++)
       for (j=0;j<6;j++) 
           mgcnum[i][j]=conmat[k++];


   hash = idxsmalloc((*ne)+1, 0, "MXNODALMETIS: hash");
   
   tot=0;
   for (i=0;i<(*ne);i++){
      hash[i]=tot;
      tot+=sizes[etype[i]];
      }

   if (*numflag == 1)
    ChangeMesh2CNumbering(tot, elmnts);




   mask = (1<<11)-1;
   mark = idxsmalloc(mask+1, -1, "GENDUALMETIS: mark");


   /* Construct the node-element list first */
   nptr = idxsmalloc((*nn)+1, 0, "MXDUALMETIS: nptr");

   for(i=0;i<(*ne);i++){
   l=hash[i];
   for (j=(sizes[etype[i]]), k=0; k<j; k++)
     nptr[elmnts[l+k]]++;

   }


   MAKECSR(i, *nn, nptr);

   nind = idxmalloc(nptr[*nn], "MXDUALMETIS: nind");
   for (k=i=0; i<(*ne); i++) {
     for (j=0; j<sizes[etype[i]]; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }

  for (i=(*nn); i>0; i--)
     nptr[i] = nptr[i-1];
     nptr[0] = 0;



   for (i=0; i<(*ne); i++) {
     for (m=j=0; j<sizes[etype[i]]; j++) {
       n = elmnts[hash[i]+j];
       for (k=nptr[n+1]-1; k>=nptr[n]; k--) {
         if ((kk = nind[k]) <= i)
           break;

         kkk = kk&mask;
         if ((l = mark[kkk]) == -1) {
           ind[m] = kk;
           wgt[m] = 1;
           mark[kkk] = m++;
         }
         else if (ind[l] == kk) {
           wgt[l]++;
         }
         else {
           for (jj=0; jj<m; jj++) {
             if (ind[jj] == kk) {
               wgt[jj]++;
               break;
             }
           }
           if (jj == m) {
             ind[m] = kk;
             wgt[m++] = 1;
           }
         }
       }
     }
  for (j=0; j<m; j++) {
       if (wgt[j] >= mgcnum[etype[i]][etype[ind[j]]]){ 
             elms[i]++;
             elms[ind[j]]++; 
             cnt+=2;
       }   
       mark[ind[j]&mask] = -1;
     }
   }


   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);
   gk_free((void **)&hash, LTERM);

 if (*numflag == 1)
    ChangeMesh2FNumbering3(tot, elmnts);
 
 return cnt; 

}






/*****************************************************************************
* This function creates a graph corresponding to the dual of a mixed element
* mesh. At this point the supported elements are triangles, tetrahedrons, 
* bricks and lines.
******************************************************************************/
void METIS_MixedMeshToDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *elms, 
idxtype *etype, idxtype *numflag,idxtype *dxadj, idxtype *dadjncy,idxtype *conmat,
idxtype custom)
{

   idxtype i, j, jj, k, kk, kkk, l, m, tot,  n, nedges, mask;
   idxtype *nptr, *nind;
   idxtype *mark, ind[200], wgt[200];
   idxtype sizes[] = {-1, 3, 4, 8, 4, 2},
       mgcnums[] = {-1, 2, 3, 4, 2};

   idxtype *hash,*mhash; 

   idxtype mgcnum[6][6] = {-1, -1, -1, -1, -1, -1,
                       -1, 2, 3, 3, 2 , 2,
                       -1, 3, 3, 3, 3, 2,
                       -1, 3, 3, 4, 4, 2,
                       -1, 2, 3, 4, 2, 2,
                       -1, 2, 2, 2, 2, 1} ; 


   if (custom==1)  /*External magic numbers supplied */
      for (i=0,k=0;i<6;i++)
       for (j=0;j<6;j++)
           mgcnum[i][j]=conmat[k++];


   hash = idxsmalloc((*ne), 0, "MXNODALMETIS: hash");
   mhash = idxsmalloc((*ne), 0, "MXNODALMETIS: hash");
   tot=0;
   for (i=0;i<(*ne);i++){
      hash[i]=tot;
      tot+=sizes[etype[i]];
      }

   if (*numflag == 1)
    ChangeMesh2CNumbering(tot, elmnts);


	   

   mask = (1<<11)-1;
   mark = idxsmalloc(mask+1, -1, "GENDUALMETIS: mark");


   /* Construct the node-element list first */
   nptr = idxsmalloc((*nn)+1, 0, "MXDUALMETIS: nptr");

   for(i=0;i<(*ne);i++){
   l=hash[i];
   for (j=(sizes[etype[i]]), k=0; k<j; k++)
     nptr[elmnts[l+k]]++;

   }


   MAKECSR(i, *nn, nptr);

   nind = idxmalloc(nptr[*nn], "MXDUALMETIS: nind");
   for (k=i=0; i<(*ne); i++) {
     for (j=0; j<sizes[etype[i]]; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }

  for (i=(*nn); i>0; i--)
     nptr[i] = nptr[i-1];
     nptr[0] = 0;
   

   dxadj[0]=0;
 
   for (i=1; i<(*ne); i++)
         mhash[i]=dxadj[i]=dxadj[i-1]+elms[i-1];  
         

   for (i=0; i<(*ne); i++) {
     for (m=j=0; j<sizes[etype[i]]; j++) {
       n = elmnts[hash[i]+j];
       for (k=nptr[n+1]-1; k>=nptr[n]; k--) {
         if ((kk = nind[k]) <= i)
           break;

         kkk = kk&mask;
         if ((l = mark[kkk]) == -1) {
           ind[m] = kk;
           wgt[m] = 1;
           mark[kkk] = m++;
         }
         else if (ind[l] == kk) {
           wgt[l]++;
         }
         else {
           for (jj=0; jj<m; jj++) {
             if (ind[jj] == kk) {
               wgt[jj]++;
               break;
             }
           }
           if (jj == m) {
             ind[m] = kk;
             wgt[m++] = 1;
           }
         }
       }
     }
  for (j=0; j<m; j++) {
       if (wgt[j] >= mgcnum[etype[i]][etype[ind[j]]]) {
         k = ind[j];
         dadjncy[dxadj[i]++] = k;
         dadjncy[dxadj[k]++] = i;
       }
       mark[ind[j]&mask] = -1;
     }
   }

   /* Go and consolidate the dxadj and dadjncy */
 for (j=i=0; i<*ne; i++) {
     for (k=mhash[i]; k<dxadj[i]; k++,j++)
       dadjncy[j] = dadjncy[k];
     dxadj[i] = j;
   }
   for (i=*ne; i>0; i--)
     dxadj[i] = dxadj[i-1];
   dxadj[0] = 0;

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);
   gk_free((void **)&hash, LTERM);

 if (*numflag == 1)
    ChangeMesh2FNumbering(tot, elmnts, *nn, dxadj, dadjncy);


}






/*****************************************************************************
* This function creates a graph corresponding to the finite element mesh. 
* At this point the supported elements are triangles, tetrahedrons.
******************************************************************************/
void METIS_MeshToNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, 
idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  idxtype esizes[] = {-1, 3, 4, 8, 4, 2};

  if (*numflag == 1)
    ChangeMesh2CNumbering((*ne)*esizes[*etype], elmnts);

  switch (*etype) {
    case 1:
      TRINODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
    case 2:
      TETNODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
    case 3:
      HEXNODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
    case 4:
      QUADNODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
    case 5:
      LINENODALMETIS(*ne, *nn, elmnts, dxadj, dadjncy);
      break;
  }

  if (*numflag == 1)
    ChangeMesh2FNumbering((*ne)*esizes[*etype], elmnts, *nn, dxadj, dadjncy);
}


/*****************************************************************************
* This function creates a graph corresponding to the finite mixed element mesh. 
* At this point the supported elements are triangles, tetrahedrons, hexahedras,
* quadilaterals and lines.
******************************************************************************/
void METIS_MixedMeshToNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, 
idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  idxtype sizes[]  = {-1, 3,4,8,4,2};
  idxtype i, j, jj, k, kk, kkk, l, m, n, nedges;
  idxtype *nptr, *nind;
  idxtype *mark;
  idxtype *hash;
  idxtype tableH[8][3] = {1, 3, 4,
                      0, 2, 5,
                      1, 3, 6,
                      0, 2, 7,
                      0, 5, 7,
                      1, 4, 6,
                      2, 5, 7,
                      3, 4, 6};

  idxtype tableQ[4][2] = {1, 3,
                      0, 2,
                      1, 3,
                      0, 2};




  hash = idxsmalloc((*ne), 0, "MXNODALMETIS: hash");
  m=0; 
  for (i=0;i<(*ne);i++){
      hash[i]=m;
      m+=sizes[etype[i]];
      }
 
  if (*numflag == 1)
    ChangeMesh2CNumbering(m, elmnts);


   
   /* Construct the node-element list first */
   nptr = idxsmalloc((*nn)+1, 0, "MXNODALMETIS: nptr");
  
   for(i=0;i<(*ne);i++){
   l=hash[i];
   for (j=(sizes[etype[i]]), k=0; k<j; k++)
     nptr[elmnts[l+k]]++;
   
   }


   MAKECSR(i, *nn, nptr);

   nind = idxmalloc(nptr[*nn], "MXNODALMETIS: nind");
   for (k=i=0; i<(*ne); i++) {
     for (j=0; j<sizes[etype[i]]; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
  
  for (i=*nn; i>0; i--)
     nptr[i] = nptr[i-1];
     nptr[0] = 0;


   mark = idxsmalloc(*nn, -1, "MXNODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<(*nn); i++) {
     mark[i] = i;
     dxadj[i+1]=dxadj[i];
     for (j=nptr[i]; j<nptr[i+1]; j++) {

     l=hash[nind[j]];

     switch(etype[nind[j]])
     {
     case 1:

        for (jj=l, k=0; k<3; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
        }
        break; 
    
     case 2:

        for (jj=l, k=0; k<4; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
       }
       break;

     
     case 3:
      
       jj=l;
       for (k=0; k<8; k++) {
         if (elmnts[jj+k] == i)
           break;
       }
       ASSERT(k != 8);

       /* You found the index, now go and put the 3 neighbors */
       kk = elmnts[jj+tableH[k][0]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+tableH[k][1]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+tableH[k][2]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       break;

     case 4:

       jj=l;
       for (k=0; k<4; k++) {
         if (elmnts[jj+k] == i)
           break;
       }
       ASSERT(k != 4);

       /* You found the index, now go and put the 2 neighbors */
       kk = elmnts[jj+tableQ[k][0]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+tableQ[k][1]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       break; 

     
     case 5:
       
        for (jj=l, k=0; k<2; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
       }
       break;

 
   
   }
   dxadj[i+1] = nedges;

   }

}
   
     

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);




  if (*numflag == 1)
    ChangeMesh2FNumbering(m, elmnts, *nn, dxadj, dadjncy);

}


/*****************************************************************************
* This function counts dual neighbours in finite element mesh
******************************************************************************/
idxtype GENDUALMETIS_COUNT(idxtype nelmnts, idxtype nvtxs, idxtype etype, idxtype *elmnts, idxtype *elms)
{
 idxtype i, j, jj, k, kk, kkk, l, m, n, nedges, mask, cnt=0;
   idxtype *nptr, *nind;
   idxtype *mark, ind[200], wgt[200];
   idxtype esize, esizes[] = {-1, 3, 4, 8, 4, 2},
       mgcnum, mgcnums[] = {-1, 2, 3, 4, 2, 1};

   mask = (1<<11)-1;
   mark = idxsmalloc(mask+1, -1, "GENDUALMETIS: mark");

   /* Get the element size and magic number for the particular element */
   esize = esizes[etype];
   mgcnum = mgcnums[etype];

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "GENDUALMETIS: nptr");
   for (j=esize*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "GENDUALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<esize; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   for (i=0; i<nelmnts; i++) {
     for (m=j=0; j<esize; j++) {
       n = elmnts[esize*i+j];
       for (k=nptr[n+1]-1; k>=nptr[n]; k--) {
         if ((kk = nind[k]) <= i)
           break;

         kkk = kk&mask;
         if ((l = mark[kkk]) == -1) {
           ind[m] = kk;
           wgt[m] = 1;
           mark[kkk] = m++;
         }
         else if (ind[l] == kk) {
           wgt[l]++;
         }
         else {
           for (jj=0; jj<m; jj++) {
             if (ind[jj] == kk) {
               wgt[jj]++;
               break;
             }
           }
           if (jj == m) {
             ind[m] = kk;
             wgt[m++] = 1;
           }
    }
       }
     }
     for (j=0; j<m; j++) {
       if (wgt[j] == mgcnum) {
         k = ind[j];
         elms[i]++;
         elms[k]++;
         cnt+=2;
       }
       mark[ind[j]&mask] = -1;
     }
   }


   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);
   
   return cnt;

}


/*****************************************************************************
* This function creates the dual of a finite element mesh
******************************************************************************/
void GENDUALMETIS(idxtype nelmnts, idxtype nvtxs, idxtype etype, idxtype *elmnts, idxtype *elms, idxtype *dxadj, idxtype *dadjncy)
{
 idxtype i, j, jj, k, kk, kkk, l, m, n, nedges, mask;
   idxtype *nptr, *nind, *mhash;
   idxtype *mark, ind[200], wgt[200];
   idxtype esize, esizes[] = {-1, 3, 4, 8, 4, 2},
       mgcnum, mgcnums[] = {-1, 2, 3, 4, 2, 1};

   mask = (1<<11)-1;
   mark = idxsmalloc(mask+1, -1, "GENDUALMETIS: mark");

   /* Get the element size and magic number for the particular element */
   esize = esizes[etype];
   mgcnum = mgcnums[etype];

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "GENDUALMETIS: nptr");
   for (j=esize*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "GENDUALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<esize; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;

   mhash = idxsmalloc(nelmnts, 0, "MXNODALMETIS: hash");
   
   dxadj[0]=0;
   for (i=1; i<nelmnts; i++)
         mhash[i]=dxadj[i]=dxadj[i-1]+elms[i-1];


   for (i=0; i<nelmnts; i++) {
     for (m=j=0; j<esize; j++) {
       n = elmnts[esize*i+j];
       for (k=nptr[n+1]-1; k>=nptr[n]; k--) {
         if ((kk = nind[k]) <= i)
           break;

         kkk = kk&mask;
         if ((l = mark[kkk]) == -1) {
           ind[m] = kk;
           wgt[m] = 1;
           mark[kkk] = m++;
         }
         else if (ind[l] == kk) {
           wgt[l]++;
         }
         else {
           for (jj=0; jj<m; jj++) {
             if (ind[jj] == kk) {
               wgt[jj]++;
               break;
             }
           }
           if (jj == m) {
             ind[m] = kk;
             wgt[m++] = 1;
           }
         }
       }
     }
     for (j=0; j<m; j++) {
       if (wgt[j] == mgcnum) {
         k = ind[j];
         dadjncy[dxadj[i]++] = k;
         dadjncy[dxadj[k]++] = i;
       }
       mark[ind[j]&mask] = -1;
     }
   }

   /* Go and consolidate the dxadj and dadjncy */
 for (j=i=0; i<nelmnts; i++) {
     for (k=mhash[i]; k<dxadj[i]; k++, j++)
       dadjncy[j] = dadjncy[k];
     dxadj[i] = j;
   }
   for (i=nelmnts; i>0; i--)
     dxadj[i] = dxadj[i-1];
   dxadj[0] = 0;

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);
   gk_free((void **)&mhash, LTERM);


}




/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void TRINODALMETIS(idxtype nelmnts, idxtype nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   idxtype i, j, jj, k, kk, kkk, l, m, n, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "TRINODALMETIS: nptr");
   for (j=3*nelmnts, i=0; i<j; i++) 
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "TRINODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<3; j++, k++) 
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "TRINODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       for (jj=3*nind[j], k=0; k<3; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
       }
     }
     dxadj[i+1] = nedges;
   }

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);

}


/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void TETNODALMETIS(idxtype nelmnts, idxtype nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   idxtype i, j, jj, k, kk, kkk, l, m, n, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "TETNODALMETIS: nptr");
   for (j=4*nelmnts, i=0; i<j; i++) 
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "TETNODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<4; j++, k++) 
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "TETNODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       for (jj=4*nind[j], k=0; k<4; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
       }
     }
     dxadj[i+1] = nedges;
   }

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);

}


/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void HEXNODALMETIS(idxtype nelmnts, idxtype nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   idxtype i, j, jj, k, kk, kkk, l, m, n, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;
   idxtype table[8][3] = {1, 3, 4,
                      0, 2, 5,
                      1, 3, 6,
                      0, 2, 7,
                      0, 5, 7,
                      1, 4, 6,
                      2, 5, 7,
                      3, 4, 6};

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "HEXNODALMETIS: nptr");
   for (j=8*nelmnts, i=0; i<j; i++) 
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "HEXNODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<8; j++, k++) 
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "HEXNODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       jj=8*nind[j];
       for (k=0; k<8; k++) {
         if (elmnts[jj+k] == i)
           break;
       }
       ASSERT(k != 8);

       /* You found the index, now go and put the 3 neighbors */
       kk = elmnts[jj+table[k][0]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+table[k][1]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+table[k][2]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
     }
     dxadj[i+1] = nedges;
   }

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);

}


/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void QUADNODALMETIS(idxtype nelmnts, idxtype nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy)
{
   idxtype i, j, jj, k, kk, kkk, l, m, n, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;
   idxtype table[4][2] = {1, 3, 
                      0, 2,
                      1, 3, 
                      0, 2}; 

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "QUADNODALMETIS: nptr");
   for (j=4*nelmnts, i=0; i<j; i++) 
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "QUADNODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<4; j++, k++) 
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "QUADNODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       jj=4*nind[j];
       for (k=0; k<4; k++) {
         if (elmnts[jj+k] == i)
           break;
       }
       ASSERT(k != 4);

       /* You found the index, now go and put the 2 neighbors */
       kk = elmnts[jj+table[k][0]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
       kk = elmnts[jj+table[k][1]];
       if (mark[kk] != i) {
         mark[kk] = i;
         dadjncy[nedges++] = kk;
       }
     }
     dxadj[i+1] = nedges;
   }

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);

}


/*****************************************************************************
* This function creates the nodal graph of a finite element mesh
******************************************************************************/
void LINENODALMETIS(idxtype nelmnts, idxtype nvtxs, idxtype *elmnts, idxtype *dxadj, idxtype *dadjncy
)
{
   idxtype i, j, jj, k, kk, kkk, l, m, n, nedges;
   idxtype *nptr, *nind;
   idxtype *mark;

   /* Construct the node-element list first */
   nptr = idxsmalloc(nvtxs+1, 0, "TRINODALMETIS: nptr");
   for (j=2*nelmnts, i=0; i<j; i++)
     nptr[elmnts[i]]++;
   MAKECSR(i, nvtxs, nptr);

   nind = idxmalloc(nptr[nvtxs], "TRINODALMETIS: nind");
   for (k=i=0; i<nelmnts; i++) {
     for (j=0; j<2; j++, k++)
       nind[nptr[elmnts[k]]++] = i;
   }
   for (i=nvtxs; i>0; i--)
     nptr[i] = nptr[i-1];
   nptr[0] = 0;


   mark = idxsmalloc(nvtxs, -1, "TRINODALMETIS: mark");

   nedges = dxadj[0] = 0;
   for (i=0; i<nvtxs; i++) {
     mark[i] = i;
     for (j=nptr[i]; j<nptr[i+1]; j++) {
       for (jj=2*nind[j], k=0; k<2; k++, jj++) {
         kk = elmnts[jj];
         if (mark[kk] != i) {
           mark[kk] = i;
           dadjncy[nedges++] = kk;
         }
       }
     }
     dxadj[i+1] = nedges;
   }

   gk_free((void **)&mark, LTERM);
   gk_free((void **)&nptr, LTERM);
   gk_free((void **)&nind, LTERM);

}

