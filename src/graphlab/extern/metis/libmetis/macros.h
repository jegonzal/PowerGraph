/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * macros.h
 *
 * This file contains macros used in multilevel
 *
 * Started 9/25/94
 * George
 *
 * $Id: macros.h,v 1.6 2003/04/30 12:42:05 karypis Exp $
 *
 */


/*************************************************************************
* The following macro returns a random number in the specified range
**************************************************************************/
#define AND(a, b) ((a) < 0 ? ((-(a))&(b)) : ((a)&(b)))
#define OR(a, b) ((a) < 0 ? -((-(a))|(b)) : ((a)|(b)))
#define XOR(a, b) ((a) < 0 ? -((-(a))^(b)) : ((a)^(b)))

#define idxcopy(n, a, b) (idxtype *)memcpy((void *)(b), (void *)(a), sizeof(idxtype)*(n)) 

#define HASHFCT(key, size) ((key)%(size))


/*************************************************************************
* Datatype related macros
**************************************************************************/
#define idxtype_abs(x) ((x) > 0 ? (x) : -(x))





/*************************************************************************
* These macros insert and remove nodes from the boundary list
**************************************************************************/
#define BNDInsert(nbnd, bndind, bndptr, vtx) \
   do { \
     ASSERT(bndptr[vtx] == -1); \
     bndind[nbnd] = vtx; \
     bndptr[vtx] = nbnd++;\
   } while(0) 

#define BNDDelete(nbnd, bndind, bndptr, vtx) \
   do { \
     ASSERT(bndptr[vtx] != -1); \
     bndind[bndptr[vtx]] = bndind[--nbnd]; \
     bndptr[bndind[nbnd]] = bndptr[vtx]; \
     bndptr[vtx] = -1; \
   } while(0) 


