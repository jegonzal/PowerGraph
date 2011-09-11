#ifndef _LDA_H
#define _LDA_H

double digamma(double x);
double trigamma(double x);

void lda_em_update_function(gl_types::iscope & scope,
			    gl_types::icallback & scheduler);

#define psi(x)  digamma(x)
#define ppsi(x) trigamma(x)
#define MAX_RECURSION_LIMIT  20
#define MAX_NEWTON_ITERATION 20



#endif

