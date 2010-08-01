/*!
\file  dfkvkselect.c
\brief Sorts only the largest k values
 
\date   Started 7/14/00
\author George
\version\verbatim $Id: dfkvkselect.c 1421 2007-04-06 14:37:41Z karypis $\endverbatim
*/


#include <GKlib.h>

/* Byte-wise swap two items of size SIZE. */
#define QSSWAP(a, b, stmp) do { stmp = (a); (a) = (b); (b) = stmp; } while (0)


/******************************************************************************
* This function puts the 'topk' largest values in the beginning of the array
*******************************************************************************/
int gk_dfkvkselect(int n, int topk, gk_fkv_t *cand)
{
  int i, j, lo, hi, mid;
  gk_fkv_t pivot, stmp;

  if (n <= topk)
    return n; /* return if the array has fewer elements than we want */

  for (lo=0, hi=n-1; hi-lo > 2;) {
    mid = lo + ((hi-lo) >> 1);

    /* select the pivot */
    if (cand[mid].key > cand[lo].key)
      QSSWAP(cand[mid], cand[lo], stmp);
    if (cand[hi].key > cand[mid].key)
      QSSWAP(cand[mid], cand[hi], stmp);
    else 
      goto jump_over;
    if (cand[mid].key > cand[lo].key)
      QSSWAP(cand[mid], cand[lo], stmp);

jump_over:;
    pivot = cand[mid];

    /* Here's the famous ``collapse the walls'' section of quicksort. */
    for (i=lo+1, j=hi-1; i<=j;) {
      for (; cand[i].key > pivot.key; i++);
      for (; pivot.key > cand[j].key; j--);

      if (i < j) {
        QSSWAP (cand[i], cand[j], stmp);
        i++; 
        j--;
      }
      else if (i == j) {
        i++;
        j--;
      }
    } 

    if (i > topk) 
      hi = i;
    else if (i < topk)
      lo = i;
    else
      break;
  }

  if (hi-lo == 2) {
    if (cand[lo].key < cand[lo+1].key)
      QSSWAP(cand[lo], cand[lo+1], stmp);
  }

/*
  if (cand[lo].key < cand[hi].key)
    printf("Hmm Error: %d %d %d %f %f\n", i, lo, hi, cand[lo].key, cand[hi].key);


  for (i=topk; i<n; i++) {
    for (j=0; j<topk; j++)
      if (cand[i].key > cand[j].key)
        printf("Hmm Error: %d %d %f %f %d %d\n", i, j, cand[i].key, cand[j].key, lo, hi);
  }
*/    

  return topk;
}
