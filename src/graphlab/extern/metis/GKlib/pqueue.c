/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pqueue.c
 *
 * This file contains functions for manipulating the bucket list
 * representation of the gains associated with each vertex in a graph.
 * These functions are used by the refinement algorithms
 *
 * Started 9/2/94
 * George
 *
 * $Id: pqueue.c 1414 2007-04-05 02:52:46Z karypis $
 *
 */

#include <GKlib.h>


/*************************************************************************
* This function initializes the data structures of the priority queue
**************************************************************************/
void gk_PQueueInit(gk_PQueue_t *queue, int maxnodes)
{
  queue->nnodes = 0;
  queue->maxnodes = maxnodes;

  queue->heap    = gk_fkvmalloc(maxnodes, "gk_PQueueInit: heap");
  queue->locator = gk_ismalloc(maxnodes, -1, "gk_PQueueInit: locator");
}


/*************************************************************************
* This function resets the buckets
**************************************************************************/
void gk_PQueueReset(gk_PQueue_t *queue)
{
  queue->nnodes = 0;

  gk_iset(queue->maxnodes, -1, queue->locator);
}


/*************************************************************************
* This function frees the buckets
**************************************************************************/
void gk_PQueueFree(gk_PQueue_t *queue)
{
  gk_free((void *)&queue->heap, &queue->locator, LTERM);
  queue->maxnodes = 0;
}


/*************************************************************************
* This function returns the number of nodes in the queue
**************************************************************************/
int gk_PQueueGetSize(gk_PQueue_t *queue)
{
  return queue->nnodes;
}


/*************************************************************************
* This function adds an item in the priority queue
**************************************************************************/
int gk_PQueueInsert(gk_PQueue_t *queue, int node, float key)
{
  int i, j;
  int *locator;
  gk_fkv_t *heap;

  ASSERT(gk_CheckHeap(queue));

  heap = queue->heap;
  locator = queue->locator;

  ASSERT(locator[node] == -1);

  i = queue->nnodes++;
  while (i > 0) {
    j = (i-1)/2;
    if (heap[j].key < key) {
      heap[i] = heap[j];
      locator[heap[i].val] = i;
      i = j;
    }
    else 
      break;
  }
  ASSERT(i >= 0);
  heap[i].key = key;
  heap[i].val = node;
  locator[node] = i;

  ASSERT(gk_CheckHeap(queue));

  return 0;
}


/*************************************************************************
* This function deletes an item from the priority queue
**************************************************************************/
int gk_PQueueDelete(gk_PQueue_t *queue, int node)
{
  int i, j;
  float newkey, oldkey;
  int *locator;
  gk_fkv_t *heap;

  heap = queue->heap;
  locator = queue->locator;

  ASSERT(locator[node] != -1);
  ASSERT(heap[locator[node]].val == node);

  ASSERT(gk_CheckHeap(queue));

  i = locator[node];
  locator[node] = -1;

  if (--queue->nnodes > 0 && heap[queue->nnodes].val != node) {
    node = heap[queue->nnodes].val;
    newkey = heap[queue->nnodes].key;
    oldkey = heap[i].key;

    if (oldkey < newkey) { /* Filter-up */
      while (i > 0) {
        j = (i-1)>>1;
        if (heap[j].key < newkey) {
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else 
          break;
      }
    }
    else { /* Filter down */
      while ((j=2*i+1) < queue->nnodes) {
        if (heap[j].key > newkey) {
          if (j+1 < queue->nnodes && heap[j+1].key > heap[j].key)
            j = j+1;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else if (j+1 < queue->nnodes && heap[j+1].key > newkey) {
          j = j+1;
          heap[i] = heap[j];
          locator[heap[i].val] = i;
          i = j;
        }
        else
          break;
      }
    }

    heap[i].key = newkey;
    heap[i].val = node;
    locator[node] = i;
  }

  ASSERT(gk_CheckHeap(queue));

  return 0;
}



/*************************************************************************
* This function deletes a node from a partition and reinserts it with
* an updated key value
**************************************************************************/
int gk_PQueueUpdate(gk_PQueue_t *queue, int node, float newkey)
{
  int i, j;
  int *locator;
  gk_fkv_t *heap;
  float oldkey;

  heap = queue->heap;
  locator = queue->locator;

  oldkey = heap[locator[node]].key;

  ASSERT(locator[node] != -1);
  ASSERT(heap[locator[node]].val == node);
  ASSERT(gk_CheckHeap(queue));

  i = locator[node];

  if (oldkey < newkey) { /* Filter-up */
    while (i > 0) {
      j = (i-1)>>1;
      if (heap[j].key < newkey) {
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else 
        break;
    }
  }
  else { /* Filter down */
    while ((j=2*i+1) < queue->nnodes) {
      if (heap[j].key > newkey) {
        if (j+1 < queue->nnodes && heap[j+1].key > heap[j].key)
          j = j+1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else if (j+1 < queue->nnodes && heap[j+1].key > newkey) {
        j = j+1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else
        break;
    }
  }

  heap[i].key = newkey;
  heap[i].val = node;
  locator[node] = i;

  ASSERT(gk_CheckHeap(queue));

  return 0;
}



/*************************************************************************
* This function returns the item with the largest key-value and removes
* it from the priority queue
**************************************************************************/
int gk_PQueueGetMax(gk_PQueue_t *queue)
{
  int vtx, i, j, node;
  int *locator;
  gk_fkv_t *heap;
  float key;

  if (queue->nnodes == 0)
    return -1;

  queue->nnodes--;

  heap = queue->heap;
  locator = queue->locator;

  vtx = heap[0].val;
  locator[vtx] = -1;

  if ((i = queue->nnodes) > 0) {
    key = heap[i].key;
    node = heap[i].val;
    i = 0;
    while ((j=2*i+1) < queue->nnodes) {
      if (heap[j].key > key) {
        if (j+1 < queue->nnodes && heap[j+1].key > heap[j].key)
          j = j+1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else if (j+1 < queue->nnodes && heap[j+1].key > key) {
        j = j+1;
        heap[i] = heap[j];
        locator[heap[i].val] = i;
        i = j;
      }
      else
        break;
    }

    heap[i].key = key;
    heap[i].val = node;
    locator[node] = i;
  }

  ASSERT(gk_CheckHeap(queue));
  return vtx;
}
      

/*************************************************************************
* This function returns the vertex with the largest gain from a partition
**************************************************************************/
int gk_PQueueSeeMaxVal(gk_PQueue_t *queue)
{
  return (queue->nnodes == 0 ? -1 : queue->heap[0].val);
}

/*************************************************************************
* This function returns the vertex with the largest gain from a partition
**************************************************************************/
float gk_PQueueSeeMaxKey(gk_PQueue_t *queue)
{
  return (queue->nnodes == 0 ? FLT_MAX : queue->heap[0].key);
}
      

/*************************************************************************
* This function returns the length of the queue
**************************************************************************/
int gk_PQueueLength(gk_PQueue_t *queue)
{
  return queue->nnodes;
}

/*************************************************************************
* This function returns the vertex with the largest gain from a partition
**************************************************************************/
int gk_PQueueSeeConstraintMax(gk_PQueue_t *queue, float maxwgt, double *wgts)
{
  int i;
  if (queue->nnodes == 0)
    return -1;

  if (maxwgt <= 1000)
    return gk_PQueueSeeMaxVal(queue);

  for (i=0; i<queue->nnodes; i++) {
    if (queue->heap[i].key > 0) {
      if (wgts[queue->heap[i].val] <= maxwgt)
        return queue->heap[i].val;
    }
    else {
      if (queue->heap[i/2].key <= 0)
        break;
    }
  }

  return queue->heap[0].val;

}


/*************************************************************************
* This function returns the key of a specific node
**************************************************************************/
float gk_PQueueSeeKey(gk_PQueue_t *queue, int node)
{
  int *locator;
  gk_fkv_t *heap;

  heap = queue->heap;
  locator = queue->locator;

  return heap[locator[node]].key;
}



/*************************************************************************
* This functions checks the consistency of the heap
**************************************************************************/
int gk_CheckHeap(gk_PQueue_t *queue)
{
  int i, j, nnodes;
  int *locator;
  gk_fkv_t *heap;

  heap = queue->heap;
  locator = queue->locator;
  nnodes = queue->nnodes;

  if (nnodes == 0)
    return 1;

  ASSERT(locator[heap[0].val] == 0);
  for (i=1; i<nnodes; i++) {
    ASSERTP(locator[heap[i].val] == i, ("%d %d %d %d\n", nnodes, i, heap[i].val, locator[heap[i].val])); 
    ASSERTP(heap[i].key <= heap[(i-1)/2].key, ("%d %d %d %f %f\n", i, (i-1)/2, nnodes, heap[i].key, heap[(i-1)/2].key));
  }
  for (i=1; i<nnodes; i++)
    ASSERT(heap[i].key <= heap[0].key);

  for (j=i=0; i<queue->maxnodes; i++) {
    if (locator[i] != -1)
      j++;
  }
  ASSERTP(j == nnodes, ("%d %d\n", j, nnodes));

  return 1;
}
