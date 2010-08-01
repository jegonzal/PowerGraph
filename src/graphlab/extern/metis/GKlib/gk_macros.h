/*!
\file gk_macros.h
\brief This file contains various macros

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_macros.h 1432 2007-04-07 20:06:19Z karypis $ \endverbatim
*/

#ifndef _GK_MACROS_H_
#define _GK_MACROS_H_

/*-------------------------------------------------------------
 * Usefull commands 
 *-------------------------------------------------------------*/
#define amax(a, b) ((a) >= (b) ? (a) : (b))
#define amin(a, b) ((a) >= (b) ? (b) : (a))
#define amax3(a, b, c) ((a) >= (b) && (a) >= (c) ? (a) : ((b) >= (a) && (b) >= (c) ? (b) : (c)))
#define SWAP(a, b, tmp) do {(tmp) = (a); (a) = (b); (b) = (tmp);} while(0) 
#define INC_DEC(a, b, val) do {(a) += (val); (b) -= (val);} while(0)
#define gk_ccopy(n, a, b) (int *)memmove((void *)(b), (void *)(a), sizeof(char)*(n)) 
#define gk_icopy(n, a, b) (int *)memmove((void *)(b), (void *)(a), sizeof(int)*(n)) 
#define gk_fcopy(n, a, b) (float *)memmove((void *)(b), (void *)(a), sizeof(float)*(n))
#define gk_dcopy(n, a, b) (double *)memmove((void *)(b), (void *)(a), sizeof(double)*(n))
#define gk_idxcopy(n, a, b) (double *)memmove((void *)(b), (void *)(a), sizeof(gk_idx_t)*(n))
#define sign(a, b) ((a >= 0 ? b : -b))

#define ONEOVERRANDMAX (1.0/(RAND_MAX+1.0))
#define RandomInRange(u) ((int) (ONEOVERRANDMAX*(u)*rand()))


/*-------------------------------------------------------------
 * Timing macros
 *-------------------------------------------------------------*/
#define gk_clearcputimer(tmr) (tmr = 0.0)
#define gk_startcputimer(tmr) (tmr -= gk_CPUSeconds())
#define gk_stopcputimer(tmr)  (tmr += gk_CPUSeconds())
#define gk_getcputimer(tmr)   (tmr)

#define gk_clearwctimer(tmr) (tmr = 0.0)
#define gk_startwctimer(tmr) (tmr -= gk_WClockSeconds())
#define gk_stopwctimer(tmr)  (tmr += gk_WClockSeconds())
#define gk_getwctimer(tmr)   (tmr)

/*-------------------------------------------------------------
 * dbglvl handling macros
 *-------------------------------------------------------------*/
#define IFSET(a, flag, cmd) if ((a)&(flag)) (cmd);


/*-------------------------------------------------------------
 * gracefull library exit macro
 *-------------------------------------------------------------*/
#define GKSETJMP() (setjmp(gk_return_to_entry))
 

/*-------------------------------------------------------------
 * Debuging memory leaks
 *-------------------------------------------------------------*/
#ifdef DMALLOC
#   define MALLOC_CHECK(ptr)                                          \
    if (malloc_verify((ptr)) == DMALLOC_VERIFY_ERROR) {  \
        printf("***MALLOC_CHECK failed on line %d of file %s: " #ptr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define MALLOC_CHECK(ptr) ;
#endif 


/*-------------------------------------------------------------
 * CSR conversion macros
 *-------------------------------------------------------------*/
#define MAKECSR(i, n, a) \
   do { \
     for (i=1; i<n; i++) a[i] += a[i-1]; \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0) 

#define SHIFTCSR(i, n, a) \
   do { \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0) 




/*-------------------------------------------------------------
 * Program Assertions
 *-------------------------------------------------------------*/
#ifdef DEBUG
#   define ASSERT(expr)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        abort();                                                \
    }
#else
#   define ASSERT(expr) ;
#endif 

#ifdef DEBUG
#   define ASSERTP(expr,msg)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        printf msg ; \
        printf("\n"); \
        abort();                                                \
    }
#else
#   define ASSERTP(expr,msg) ;
#endif 

#endif
