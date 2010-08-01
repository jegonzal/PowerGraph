/*!
\file  memory.c
\brief This file contains various allocation routines 

The allocation routines included are for 1D and 2D arrays of the 
most datatypes that GKlib support. Many of these routines are 
defined with the help of the macros in gk_memory.h. These macros 
can be used to define other memory allocation routines.

\date   Started 4/3/2007
\author George
\version\verbatim $Id: memory.c 1432 2007-04-07 20:06:19Z karypis $ \endverbatim
*/


#include <GKlib.h>

#ifdef GKMSPACE
/* This is the mspace for all the gk_malloc() calls. This is a thread local allocation */
#ifdef __linux__
static __thread mspace gk_mspace = 0;
#else
static  mspace gk_mspace = 0;
#endif

/* This function is mostly for debugging */
void gk_printmspaceaddr() { printf("mspace: %p\n", (void *)gk_mspace); }
#endif


/*************************************************************************/
/*! Routines for allocating an array of a particular data type */
/**************************************************************************/
GK_XMALLOC(gk_cmalloc, char)
GK_XMALLOC(gk_imalloc, int)
GK_XMALLOC(gk_fmalloc, float)
GK_XMALLOC(gk_dmalloc, double)
GK_XMALLOC(gk_idxmalloc, gk_idx_t)

GK_XMALLOC(gk_ckvmalloc, gk_ckv_t)
GK_XMALLOC(gk_ikvmalloc, gk_ikv_t)
GK_XMALLOC(gk_fkvmalloc, gk_fkv_t)
GK_XMALLOC(gk_dkvmalloc, gk_dkv_t)
GK_XMALLOC(gk_skvmalloc, gk_skv_t)
GK_XMALLOC(gk_idxkvmalloc, gk_idxkv_t)


/*************************************************************************/
/*! Routines for reallocating an array of a particular data type */
/**************************************************************************/
GK_XREALLOC(gk_crealloc, char)
GK_XREALLOC(gk_irealloc, int)
GK_XREALLOC(gk_frealloc, float)
GK_XREALLOC(gk_drealloc, double)
GK_XREALLOC(gk_idxrealloc, gk_idx_t)

GK_XREALLOC(gk_ckvrealloc, gk_ckv_t)
GK_XREALLOC(gk_ikvrealloc, gk_ikv_t)
GK_XREALLOC(gk_fkvrealloc, gk_fkv_t)
GK_XREALLOC(gk_dkvrealloc, gk_dkv_t)
GK_XREALLOC(gk_skvrealloc, gk_skv_t)
GK_XREALLOC(gk_idxkvrealloc, gk_idxkv_t)


/*************************************************************************/
/*! Routines for allocating and setting array of a particular data type */
/**************************************************************************/
GK_XSMALLOC(gk_csmalloc, char, gk_cset)
GK_XSMALLOC(gk_ismalloc, int, gk_iset)
GK_XSMALLOC(gk_fsmalloc, float, gk_fset)
GK_XSMALLOC(gk_dsmalloc, double, gk_dset)
GK_XSMALLOC(gk_idxsmalloc, gk_idx_t, gk_idxset)


/*************************************************************************/
/*! Routines for allocating/setting 2D arrays of a particular data type */
/**************************************************************************/
GK_ALLOCMATRIX(gk_cAllocMatrix, char, gk_csmalloc)
GK_ALLOCMATRIX(gk_iAllocMatrix, int, gk_ismalloc)
GK_ALLOCMATRIX(gk_fAllocMatrix, float, gk_fsmalloc)
GK_ALLOCMATRIX(gk_dAllocMatrix, double, gk_dsmalloc)
GK_ALLOCMATRIX(gk_idxAllocMatrix, gk_idx_t, gk_idxsmalloc)


/*************************************************************************/
/*! Routines for freeing 2D arrays of a particular data type */
/**************************************************************************/
GK_FREEMATRIX(gk_cFreeMatrix, char)
GK_FREEMATRIX(gk_iFreeMatrix, int)
GK_FREEMATRIX(gk_fFreeMatrix, float)
GK_FREEMATRIX(gk_dFreeMatrix, double)
GK_FREEMATRIX(gk_idxFreeMatrix, gk_idx_t)


/*************************************************************************/
/*! Routines for filling a 2D array with a constant value */
/**************************************************************************/
GK_SETMATRIX(gk_cSetMatrix, char)
GK_SETMATRIX(gk_iSetMatrix, int)
GK_SETMATRIX(gk_fSetMatrix, float)
GK_SETMATRIX(gk_dSetMatrix, double)
GK_SETMATRIX(gk_idxSetMatrix, gk_idx_t)




/*************************************************************************
* This function allocates a two-dimensional matrix 
**************************************************************************/
void gk_AllocMatrix(void ***r_matrix, size_t elmlen, size_t ndim1, size_t ndim2)
{
  gk_idx_t i;
  void **matrix;

  matrix = (void **)gk_malloc(ndim1*sizeof(void *), "GKAllocMatrix: matrix");
  for (i=0; i<ndim1; i++) 
    matrix[i] = (void *)gk_malloc(ndim2*elmlen, "GKAllocMatrix: matrix[i]");

  *r_matrix = matrix;
}


/*************************************************************************
* This function frees a two-dimensional matrix 
**************************************************************************/
void gk_FreeMatrix(void ***r_matrix, size_t ndim1, size_t ndim2)
{
  gk_idx_t i;
  void **matrix;

  matrix = *r_matrix;

  for (i=0; i<ndim1; i++) 
    gk_free((void **)&matrix[i], LTERM);

  gk_free((void **)matrix, LTERM);

  *r_matrix = NULL;

}



/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *gk_malloc(size_t nbytes, char *msg)
{
  return malloc(nbytes);
  void *ptr;

  if (nbytes == 0)
    return NULL;

#ifdef GKMSPACE
  if (gk_mspace == 0)
    gk_mspace = create_mspace(0, 0);

  if (gk_mspace == NULL)
    gk_errexit(SIGMEM, "***Memory allocation failed for creating gk_mspace.");

  ptr = (void *)mspace_malloc(gk_mspace, nbytes);
#else
  ptr = (void *)dlmalloc(nbytes);
#endif

  if (ptr == NULL) {
#ifdef WIN32
    printf("   Maximum memory used:              %10Iu bytes\n", gk_GetMaxMemoryUsed());
    printf("   Current memory used:              %10Iu bytes\n", gk_GetCurMemoryUsed());
#endif
#ifdef SUNOS
    printf("   Maximum memory used:              %10lu bytes\n", gk_GetMaxMemoryUsed());
    printf("   Current memory used:              %10lu bytes\n", gk_GetCurMemoryUsed());
#endif
#ifdef LINUX
    printf("   Maximum memory used:              %10zu bytes\n", gk_GetMaxMemoryUsed());
    printf("   Current memory used:              %10zu bytes\n", gk_GetCurMemoryUsed());
#endif

    gk_errexit(SIGMEM, "***Memory allocation failed for %s. Requested size: %zd bytes", msg, nbytes);
  }

  return ptr;
}


/*************************************************************************
* This function is my wrapper around realloc
**************************************************************************/
void *gk_realloc(void *oldptr, size_t nbytes, char *msg)
{
  return realloc(oldptr, nbytes);
  void *ptr;

  if (nbytes == 0) {
    gk_free((void **)&oldptr, LTERM);
    return NULL;
  }

#ifdef GKMSPACE
  if (gk_mspace == 0)
    gk_mspace = create_mspace(0, 0);

  if (gk_mspace == NULL)
    gk_errexit(SIGMEM, "***Memory allocation failed for creating gk_mspace.");

  ptr = (void *)mspace_realloc(gk_mspace, oldptr, nbytes);
#else
  ptr = (void *)dlrealloc(oldptr, nbytes);
#endif

  if (ptr == NULL) {
#ifdef WIN32
    printf("   Maximum memory used:              %10Iu bytes\n", gk_GetMaxMemoryUsed());
    printf("   Current memory used:              %10Iu bytes\n", gk_GetCurMemoryUsed());
#endif
#ifdef SUNOS
    printf("   Maximum memory used:              %10lu bytes\n", gk_GetMaxMemoryUsed());
    printf("   Current memory used:              %10lu bytes\n", gk_GetCurMemoryUsed());
#endif
#ifdef LINUX
    printf("   Maximum memory used:              %10zu bytes\n", gk_GetMaxMemoryUsed());
    printf("   Current memory used:              %10zu bytes\n", gk_GetCurMemoryUsed());
#endif

    gk_errexit(SIGMEM, "***Memory re-allocation failed for %s. Requested size: %zd bytes", msg, nbytes);
  }

  return ptr;
}


/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void gk_free(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
    free(*ptr1);
  *ptr1 = NULL;

  va_start(plist, ptr1);

  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
      free(*ptr);
    *ptr = NULL;
  }

  va_end(plist);
}            


/*************************************************************************
* This function cleans up the memory that has been allocated thus far.
* This work only if code has been compiled with GKMSPACE 
**************************************************************************/
void gk_malloc_cleanup()
{
#ifdef GKMSPACE
  if (gk_mspace != 0)
    destroy_mspace(gk_mspace);
  gk_mspace = 0;
#endif
}


/*************************************************************************
* This function returns the current ammount of dynamically allocated
* memory that is used by the system
**************************************************************************/
size_t gk_GetCurMemoryUsed()
{
  struct mallinfo meminfo;
  size_t cused=0;

#ifdef GKMSPACE
  if (gk_mspace != 0) {
    meminfo = mspace_mallinfo(gk_mspace);
    cused = meminfo.uordblks;
  }
#else
  meminfo = dlmallinfo();
  cused = meminfo.uordblks;
#endif

  return cused;
}


/*************************************************************************
* This function returns the maximum ammount of dynamically allocated 
* memory that was used by the system
**************************************************************************/
size_t gk_GetMaxMemoryUsed()
{
  struct mallinfo meminfo;
  size_t mused=0;

#ifdef GKMSPACE
  if (gk_mspace != 0) {
    meminfo = mspace_mallinfo(gk_mspace);
    mused = meminfo.usmblks;
  }
#else
  meminfo = dlmallinfo();
  mused = meminfo.usmblks;
#endif

  return mused;
}


