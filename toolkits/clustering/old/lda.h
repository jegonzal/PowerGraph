#ifndef _LDA_H
#define _LDA_H



/**
 * This code implements Blei's LDA algorithm as described on:
 * Blei, David M.; Ng, Andrew Y.; Jordan, Michael I (January 2003). Lafferty, John. ed. "Latent Dirichlet allocation". Journal of Machine Learning Research 3 (4â€“5): pp. 993â€“1022. doi:10.1162/jmlr.2003.3.4-5.993.
 * The code is based on the code: http://chasen.org/~daiti-m/dist/lda/
 * by Daichi Mochihashi.
 *
 * */
double digamma(double x);
double trigamma(double x);

void lda_em_update_function(gl_types::iscope & scope,
			    gl_types::icallback & scheduler);

#define psi(x)  digamma(x)
#define ppsi(x) trigamma(x)
#define MAX_RECURSION_LIMIT  20
#define MAX_NEWTON_ITERATION 20

   
//intermediate memory space allocated for each computing thread
//to preform computation
struct scratch_buffer{
  double * ap;
  double * nt;
  double * pnt;
  double **q;

  scratch_buffer(){
    ap = NULL;
    nt = NULL;
    pnt = NULL;
    q = NULL;
  }
  ~scratch_buffer();
};



#endif

