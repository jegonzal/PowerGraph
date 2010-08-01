/*!
\file gk_proto.h
\brief This file contains function prototypes

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_proto.h 1432 2007-04-07 20:06:19Z karypis $ \endverbatim
*/

#ifndef _GK_PROTO_H_
#define _GK_PROTO_H_


/*-------------------------------------------------------------
 * blas.c 
 *-------------------------------------------------------------*/
GK_SET_PROTO(gk_cset, char)
GK_SET_PROTO(gk_iset, int)
GK_SET_PROTO(gk_fset, float)
GK_SET_PROTO(gk_dset, double)
GK_SET_PROTO(gk_idxset, gk_idx_t)
GK_INCSET_PROTO(gk_cincset, char)
GK_INCSET_PROTO(gk_iincset, int)
GK_INCSET_PROTO(gk_fincset, float)
GK_INCSET_PROTO(gk_dincset, double)
GK_INCSET_PROTO(gk_idxincset, gk_idx_t)
GK_ARGMAX_PROTO(gk_cargmax, char)
GK_ARGMAX_PROTO(gk_iargmax, int)
GK_ARGMAX_PROTO(gk_fargmax, float)
GK_ARGMAX_PROTO(gk_dargmax, double)
GK_ARGMAX_PROTO(gk_idxargmax, gk_idx_t)
GK_ARGMIN_PROTO(gk_cargmin, char)
GK_ARGMIN_PROTO(gk_iargmin, int)
GK_ARGMIN_PROTO(gk_fargmin, float)
GK_ARGMIN_PROTO(gk_dargmin, double)
GK_ARGMIN_PROTO(gk_idxargmin, gk_idx_t)
GK_ARGMAX_N_PROTO(gk_cargmax_n, char)
GK_ARGMAX_N_PROTO(gk_iargmax_n, int)
GK_ARGMAX_N_PROTO(gk_fargmax_n, float)
GK_ARGMAX_N_PROTO(gk_dargmax_n, double)
GK_ARGMAX_N_PROTO(gk_idxargmax_n, gk_idx_t)
GK_SUM_PROTO(gk_csum, char, intmax_t)
GK_SUM_PROTO(gk_isum, int, intmax_t)
GK_SUM_PROTO(gk_fsum, float, float)
GK_SUM_PROTO(gk_dsum, double, double)
GK_SUM_PROTO(gk_idxsum, gk_idx_t, intmax_t)
GK_SCALE_PROTO(gk_cscale, char)
GK_SCALE_PROTO(gk_iscale, int)
GK_SCALE_PROTO(gk_fscale, float)
GK_SCALE_PROTO(gk_dscale, double)
GK_SCALE_PROTO(gk_idxscale, gk_idx_t)
GK_NORM2_PROTO(gk_cnorm2, char, intmax_t)
GK_NORM2_PROTO(gk_inorm2, int, intmax_t)
GK_NORM2_PROTO(gk_fnorm2, float, float)
GK_NORM2_PROTO(gk_dnorm2, double, double)
GK_NORM2_PROTO(gk_idxnorm2, gk_idx_t, intmax_t)
GK_DOT_PROTO(gk_cdot, char, intmax_t)
GK_DOT_PROTO(gk_idot, int, intmax_t)
GK_DOT_PROTO(gk_fdot, float, float)
GK_DOT_PROTO(gk_ddot, double, double)
GK_DOT_PROTO(gk_idxdot, gk_idx_t, intmax_t)
GK_AXPY_PROTO(gk_caxpy, char)
GK_AXPY_PROTO(gk_iaxpy, int)
GK_AXPY_PROTO(gk_faxpy, float)
GK_AXPY_PROTO(gk_daxpy, double)
GK_AXPY_PROTO(gk_idxaxpy, gk_idx_t)




/*-------------------------------------------------------------
 * io.c
 *-------------------------------------------------------------*/
FILE *gk_fopen(char *, char *, const char *);
void gk_fclose(FILE *);
gk_loop_t gk_getline(char **lineptr, size_t *n, FILE *stream);



/*-------------------------------------------------------------
 * fs.c
 *-------------------------------------------------------------*/
int gk_fexists(char *);
int gk_dexists(char *);
intmax_t gk_getfsize(char *);
void gk_getfilestats(char *fname, int *r_nlines, int *r_ntokens, int *r_nbytes);
char *gk_getbasename(char *path);
char *gk_getextname(char *path);
char *gk_getfilename(char *path);
char *gk_getpathname(char *path);
int gk_mkpath(char *);
int gk_rmpath(char *);



/*-------------------------------------------------------------
 * memory.c
 *-------------------------------------------------------------*/
GK_XMALLOC_PROTO(gk_cmalloc, char)
GK_XMALLOC_PROTO(gk_imalloc, int)
GK_XMALLOC_PROTO(gk_fmalloc, float)
GK_XMALLOC_PROTO(gk_dmalloc, double) 
GK_XMALLOC_PROTO(gk_idxmalloc, gk_idx_t) 
GK_XMALLOC_PROTO(gk_ckvmalloc, gk_ckv_t)
GK_XMALLOC_PROTO(gk_ikvmalloc, gk_ikv_t)
GK_XMALLOC_PROTO(gk_fkvmalloc, gk_fkv_t)
GK_XMALLOC_PROTO(gk_dkvmalloc, gk_dkv_t)
GK_XMALLOC_PROTO(gk_skvmalloc, gk_skv_t)
GK_XMALLOC_PROTO(gk_idxkvmalloc, gk_idxkv_t)

GK_XREALLOC_PROTO(gk_crealloc, char)
GK_XREALLOC_PROTO(gk_irealloc, int)
GK_XREALLOC_PROTO(gk_frealloc, float)
GK_XREALLOC_PROTO(gk_drealloc, double) 
GK_XREALLOC_PROTO(gk_idxrealloc, gk_idx_t) 
GK_XREALLOC_PROTO(gk_ckvrealloc, gk_ckv_t)
GK_XREALLOC_PROTO(gk_ikvrealloc, gk_ikv_t)
GK_XREALLOC_PROTO(gk_fkvrealloc, gk_fkv_t)
GK_XREALLOC_PROTO(gk_dkvrealloc, gk_dkv_t)
GK_XREALLOC_PROTO(gk_skvrealloc, gk_skv_t)
GK_XREALLOC_PROTO(gk_idxkvrealloc, gk_idxkv_t)

GK_XSMALLOC_PROTO(gk_csmalloc, char)
GK_XSMALLOC_PROTO(gk_ismalloc, int)
GK_XSMALLOC_PROTO(gk_fsmalloc, float)
GK_XSMALLOC_PROTO(gk_dsmalloc, double)
GK_XSMALLOC_PROTO(gk_idxsmalloc, gk_idx_t)
GK_ALLOCMATRIX_PROTO(gk_cAllocMatrix, char)
GK_ALLOCMATRIX_PROTO(gk_iAllocMatrix, int)
GK_ALLOCMATRIX_PROTO(gk_fAllocMatrix, float)
GK_ALLOCMATRIX_PROTO(gk_dAllocMatrix, double)
GK_ALLOCMATRIX_PROTO(gk_idxAllocMatrix, gk_idx_t)
GK_FREEMATRIX_PROTO(gk_cFreeMatrix, char)
GK_FREEMATRIX_PROTO(gk_iFreeMatrix, int)
GK_FREEMATRIX_PROTO(gk_fFreeMatrix, float)
GK_FREEMATRIX_PROTO(gk_dFreeMatrix, double)
GK_FREEMATRIX_PROTO(gk_idxFreeMatrix, gk_idx_t)
GK_SETMATRIX_PROTO(gk_cSetMatrix, char)
GK_SETMATRIX_PROTO(gk_iSetMatrix, int)
GK_SETMATRIX_PROTO(gk_fSetMatrix, float)
GK_SETMATRIX_PROTO(gk_dSetMatrix, double)
GK_SETMATRIX_PROTO(gk_idxSetMatrix, gk_idx_t)



void            gk_AllocMatrix(void ***, size_t, size_t , size_t);
void            gk_FreeMatrix(void ***, size_t, size_t);
void           *gk_malloc(size_t, char *);
void           *gk_realloc(void *, size_t, char *);
void            gk_free(void **ptr1,...);
void            gk_malloc_cleanup();
size_t          gk_GetCurMemoryUsed();
size_t          gk_GetMaxMemoryUsed();



/*-------------------------------------------------------------
 * seq.c
 *-------------------------------------------------------------*/
gk_seq_t *gk_seq_ReadGKMODPSSM(char *file_name);
gk_i2cc2i_t *gk_i2cc2i_create_common(char *alphabet);
void gk_seq_init(gk_seq_t *seq);




/*-------------------------------------------------------------
 * pdb.c
 *-------------------------------------------------------------*/
char gk_threetoone(char *res);
void gk_freepdbf(pdbf *p);
pdbf *gk_readpdbfile(char *fname);
void gk_writefullatom(pdbf *p, char *fname);
void gk_writebackbone(pdbf *p, char *fname);
void gk_writealphacarbons(pdbf *p, char *fname);
void gk_showcorruption(pdbf *p);


/*-------------------------------------------------------------
 * error.c
 *-------------------------------------------------------------*/
void errexit(char *,...);
void gk_errexit(int signum, char *,...);
void gk_SetSignalHandlers();
void gk_UnsetSignalHandlers();
void gk_NonLocalExit_Handler(int signum);
char *gk_strerror(int errnum);
void PrintBackTrace();


/*-------------------------------------------------------------
 * util.c
 *-------------------------------------------------------------*/
void  gk_RandomPermute(size_t, int *, int);
void  gk_array2csr(size_t n, size_t range, int *array, int *ptr, int *ind);
int   gk_log2(int);
int   gk_ispow2(int);
float gk_flog2(float);


/*-------------------------------------------------------------
 * time.c
 *-------------------------------------------------------------*/
gk_wclock_t gk_WClockSeconds(void);
double gk_CPUSeconds(void);

/*-------------------------------------------------------------
 * string.c
 *-------------------------------------------------------------*/
char   *gk_strchr_replace(char *str, char *fromlist, char *tolist);
int     gk_strstr_replace(char *str, char *pattern, char *replacement, char *options, char **new_str);
char   *gk_strtprune(char *, char *);
char   *gk_strhprune(char *, char *);
char   *gk_strtoupper(char *); 
char   *gk_strtolower(char *); 
char   *gk_strdup(char *orgstr);
int     gk_strcasecmp(char *s1, char *s2);
char   *gk_time2str(time_t time);
time_t  gk_str2time(char *str);
int     gk_GetStringID(gk_StringMap_t *strmap, char *key);



/*-------------------------------------------------------------
 * sort.c 
 *-------------------------------------------------------------*/
void gk_icsort(int, char *);
void gk_dcsort(int, char *);
void gk_iisort(int, int *);
void gk_disort(int, int *);
void gk_ifsort(int, float *);
void gk_dfsort(int, float *);
void gk_idsort(int, double *);
void gk_ddsort(int, double *);
void gk_iidxsort(int, gk_idx_t *);
void gk_didxsort(int, gk_idx_t *);
void gk_ickvsort(int, gk_ckv_t *);
void gk_dckvsort(int, gk_ckv_t *);
void gk_iikvsort(int, gk_ikv_t *);
void gk_dikvsort(int, gk_ikv_t *);
void gk_ifkvsort(int, gk_fkv_t *);
void gk_dfkvsort(int, gk_fkv_t *);
void gk_idkvsort(int, gk_dkv_t *);
void gk_ddkvsort(int, gk_dkv_t *);
void gk_iskvsort(int, gk_skv_t *);
void gk_dskvsort(int, gk_skv_t *);
void gk_iidxkvsort(int, gk_idxkv_t *);
void gk_didxkvsort(int, gk_idxkv_t *);

/*-------------------------------------------------------------
 * Selection routines
 *-------------------------------------------------------------*/
int  gk_dfkvkselect(int, int, gk_fkv_t *);


/*-------------------------------------------------------------
 * Priority queue 
 *-------------------------------------------------------------*/
void  gk_PQueueInit(gk_PQueue_t *, int);
void  gk_PQueueReset(gk_PQueue_t *);
void  gk_PQueueFree(gk_PQueue_t *);
int   gk_PQueueGetSize(gk_PQueue_t *);
int   gk_PQueueInsert(gk_PQueue_t *, int, float);
int   gk_PQueueDelete(gk_PQueue_t *, int);
int   gk_PQueueUpdate(gk_PQueue_t *, int, float);
int   gk_PQueueGetMax(gk_PQueue_t *);
int   gk_PQueueSeeMaxVal(gk_PQueue_t *queue);
float gk_PQueueSeeMaxKey(gk_PQueue_t *queue);
int   gk_PQueueLength(gk_PQueue_t *);
int   gk_PQueueSeeConstraintMax(gk_PQueue_t *queue, float maxwgt, double *wgts);
float gk_PQueueSeeKey(gk_PQueue_t *queue, int node);
int   gk_CheckHeap(gk_PQueue_t *);

/*-------------------------------------------------------------
 * HTable routines
 *-------------------------------------------------------------*/
gk_HTable_t *HTable_Create(int nelements);
void         HTable_Reset(gk_HTable_t *htable);
void         HTable_Resize(gk_HTable_t *htable, int nelements);
void         HTable_Insert(gk_HTable_t *htable, int key, int val);
void         HTable_Delete(gk_HTable_t *htable, int key);
int          HTable_Search(gk_HTable_t *htable, int key);
int          HTable_GetNext(gk_HTable_t *htable, int key, int *val, int type);
int          HTable_SearchAndDelete(gk_HTable_t *htable, int key);
void         HTable_Destroy(gk_HTable_t *htable);
int          HTable_HFunction(int nelements, int key);
 

/*-------------------------------------------------------------
 * Tokenizer routines
 *-------------------------------------------------------------*/
void gk_strtokenize(char *line, char *delim, gk_Tokens_t *tokens);
void gk_freetokenslist(gk_Tokens_t *tokens);

/*-------------------------------------------------------------
 * Encoder/Decoder
 *-------------------------------------------------------------*/
void encodeblock(unsigned char *in, unsigned char *out);
void decodeblock(unsigned char *in, unsigned char *out);
void GKEncodeBase64(int nbytes, unsigned char *inbuffer, unsigned char *outbuffer);
void GKDecodeBase64(int nbytes, unsigned char *inbuffer, unsigned char *outbuffer);

/*-------------------------------------------------------------
 * OpenMP fake functions
 *-------------------------------------------------------------*/
#if !defined(__OPENMP__)
void omp_set_num_threads(int num_threads);
int omp_get_num_threads(void);
int omp_get_max_threads(void);
int omp_get_thread_num(void);
int omp_get_num_procs(void);
int omp_in_parallel(void);
void omp_set_dynamic(int num_threads);
int omp_get_dynamic(void);
void omp_set_nested(int nested);
int omp_get_nested(void);
#endif /* __OPENMP__ */

#endif

