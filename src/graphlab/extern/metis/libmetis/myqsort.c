/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * myqsort.c
 * 
 * This file contains a fast idxtype increasing qsort algorithm.
 * Addopted from TeX
 * 
 * Started 10/18/96
 * George
 * 
 * $Id: myqsort.c,v 1.2 2002/08/10 06:29:33 karypis Exp $
 */

#include <metislib.h>			/* only for type declarations */




/*************************************************************************
* Entry point of idxtype increasing sort
**************************************************************************/
void iidxsort(size_t n, idxtype *base)
{
#define idxtype_lt(a, b) ((*a) < (*b))
  GKQSORT(idxtype, base, n, idxtype_lt);
}


/*************************************************************************
* Entry point of KeyVal increasing sort, ONLY key part
**************************************************************************/
void ikeysort(size_t n, KeyValueType *base)
{
#define keyvalue_lt(a, b) ((a)->key < (b)->key)
  GKQSORT(KeyValueType, base, n, keyvalue_lt);
}



/*************************************************************************
* Entry point of KeyVal increasing sort, BOTH key and val part
**************************************************************************/
void ikeyvalsort(size_t n, KeyValueType *base)
{
#define keyvalueboth_lt(a, b) ((a)->key < (b)->key || ((a)->key == (b)->key && (a)->val < (b)->val))
  GKQSORT(KeyValueType, base, n, keyvalueboth_lt);
}


/*************************************************************************
* Entry point of DKeyValueType increasing sort based on keys
**************************************************************************/
void idkeysort(size_t n, DKeyValueType *base)
{
#define dkeyvalue_lt(a, b) ((a)->key < (b)->key)
  GKQSORT(DKeyValueType, base, n, dkeyvalue_lt);
}
