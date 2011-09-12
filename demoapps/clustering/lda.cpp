#include "clustering.h"
#include "lda.h"
#include "../gabp/advanced_config.h"

extern problem_setup ps; 
extern advanced_config ac;
//void accum_gammas (double **gammas, double *_gamma, int n, int K);
void accum_betas (double **betas, int K, vertex_data & data);
double ** dmatrix (int rows, int cols);
void free_dmatrix (double **matrix, int rows);
void newton_alpha (double *alpha, double **gammas, int M, int K, int level);
void normalize_matrix_col (double **dst, double **src, int rows, int cols);
void normalize_matrix_row (double **dst, double **src, int rows, int cols);
double lda_lik (double **beta, double **gammas, int m);
int converged (double *u, double *v, int n, double threshold);
void fill_output_lda();

double *alpha;
double **beta;
double  **gammas, **betas;

scratch_buffer * scratch;

scratch_buffer::~scratch_buffer(){
    free(ap);
    free(nt);
    free(pnt);
    free_dmatrix(q, ps.N);
}

char * rtime (double t)
{
	int hour, min, sec;
	static char buf[BUFSIZ];

	hour = (int)floor((int)t / 60 / 60);
	min  = (int)floor(((int)t % (60 * 60)) / 60);
	sec  = (int)floor((int)t % 60);
	sprintf(buf, "%2d:%02d:%02d", hour, min, sec);
	
	return buf;
}


int doublecmp (double *x, double *y)
{
        return (*x == *y) ? 0 : ((*x < *y) ? 1 : -1);
}

void init_scratch_buffers(){

      	if ((scratch = (scratch_buffer *)calloc(ac.ncpus, sizeof(scratch_buffer))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate ap.\n");
		return;
	}

        for (int i=0; i< ac.ncpus; i++){
	if ((scratch[i].ap = (double *)calloc(ac.K, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate ap.\n");
		return;
	}
	if ((scratch[i].nt = (double *)calloc(ac.K, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate nt.\n");
		return;
	}
	if ((scratch[i].pnt = (double*)calloc(ac.K, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate pnt.\n");
		return;
	}
	if ((scratch[i].q = dmatrix(ps.N, ac.K)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate q.\n");
		return;
	}
 
        }
}


void
lda_learn (double *alpha, double **beta)
{
	double lik, plik = 0;
	double z;
	int i, j, n;
	int elapsed;
	
	/*
	 *  randomize a seed
	 *
	 */
	srand(time(NULL));
	
	/*
	 *  initialize parameters
	 *
	 */

        n=ps.M;

	for (i = 0; i < ac.K; i++)
		alpha[i] = itpp::randu();
	for (i = 0, z = 0; i < ac.K; i++)
		z += alpha[i];
	for (i = 0; i < ac.K; i++)
		alpha[i] = alpha[i] / z;
	qsort(alpha, ac.K, sizeof(double), // sort alpha initially
	      (int (*)(const void *, const void *))doublecmp);

	for (i = 0; i < ps.N; i++)
		for (j = 0; j < ac.K; j++)
			beta[i][j] = (double) 1 / ps.N;

	/*
	 *  initialize posteriors
	 *
	 */
	if ((gammas = dmatrix(ps.M, ac.K)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate gammas.\n");
		return;
	}
	if ((betas = dmatrix(ps.N, ac.K)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate betas.\n");
		return;
	}


	init_scratch_buffers();
	/*
	 *  initialize buffers
	 *
	 */
	/*if ((q = dmatrix(ps.N, ac.K)) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate q.\n");
		return;
	}
	if ((_gamma = (double *)calloc(ac.K, sizeof(double))) == NULL)
	{
		fprintf(stderr, "lda_learn:: cannot allocate _gamma.\n");
		return;
	}*/
	/*if ((ap = (double *)calloc(ac.K, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate ap.\n");
		return;
	}
	if ((nt = (double *)calloc(ac.K, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate nt.\n");
		return;
	}
	if ((pnt = (double*)calloc(ac.K, sizeof(double))) == NULL) {
		fprintf(stderr, "lda_learn:: cannot allocate pnt.\n");
		return;
	}*/

	printf("Number of documents          = %d\n", ps.M);
	printf("Number of words              = %d\n", ps.N);
	printf("Number of latent classes     = %d\n", ac.K);
	printf("Number of outer EM iteration = %d\n", ac.iter);
	printf("Number of inner EM iteration = %d\n", ac.em_max_inner_iter);
	printf("Convergence threshold        = %g\n", ac.threshold);

	/*
	 *  learn main
	 *
	 */
	ps.gt.start();
        int & t = ps.iiter;
	for (t = 0; t < ac.iter; t++)
	{
		printf("%g) iteration %d/%d.. string E-Step\n", ps.gt.current_time(), t + 1, ac.iter);
		/*
		 *  VB-E step
		 *
		 */

                ps.glcore->start();
	        ps.glcore->add_task_to_all(lda_em_update_function, 1);

		printf("%g) iteration %d/%d.. accumulating beta\n", ps.gt.current_time(), t + 1, ac.iter);
		/* iterate for data */
		for (i = 0; i < ps.M; i++)
		{
			//vbem(dp, gamma, q, nt, pnt, ap,
			//     alpha, (const double **)beta, dp->len, ac.K, ac.em_max_inner_iter);
			//accum_gammas(gammas, _gamma, i, ac.K);
			accum_betas(betas, ac.K, ps.g->vertex_data(i));
		}
		printf("%g) iteration %d/%d.. starting M-step\n", ps.gt.current_time(), t + 1, ac.iter);
		/*
		 *  VB-M step
		 *
		 */
		/* Newton-Raphson for alpha */
		newton_alpha(alpha, gammas, n, ac.K, 0);
		/* MLE for beta */
		normalize_matrix_col(beta, betas, ps.N, ac.K);
		/* clean buffer */
		for (i = 0; i < ps.N; i++)
			for (j = 0; j < ac.K; j++)
				betas[i][j] = 0;
		/*
		 *  converge?
		 *
		 */
		lik = lda_lik(beta, gammas, ps.M);
		elapsed = ps.gt.current_time();
		logstream(LOG_INFO) << "likelihood " << lik << std::endl;
		if ((t > 1) && (fabs((lik - plik)/lik) < ac.threshold)) {
			if (t < 5) {
				free_dmatrix(gammas, n);
				free_dmatrix(betas, ps.N);
				//free(_gamma);
				printf("\nearly convergence. restarting..\n");
				lda_learn (alpha, beta);
				return;
			} else {
				printf("\nconverged. [%s]\n", rtime(elapsed));
				break;
			}
		}
		plik = lik;
		/* 
		 * ETA
		 *
		 */
		printf("%g) ETA:%s (%d sec/step)\r",
		       ps.gt.current_time(),
		       rtime(elapsed * ((double) ac.iter / (t + 1) - 1)),
		       (int)((double) elapsed / (t + 1) + 0.5));
	}
 

        fill_output_lda();
       
	free_dmatrix(gammas, n);
	free_dmatrix(betas, ps.N);
        free(scratch);	
	return;
}


void fill_output_lda(){

     ps.output_clusters = zeros(ps.K, ps.N);
        //for (int i=0; i<ps.K; i++)
        //  ps.output_clusters.set_row(i, ps.clusts.cluster_vec[i].location);
     //TODO
     //
     ps.output_assignements = zeros(ps.M, ac.K);
     for (int i=0; i< ps.M; i++){ 
        ps.output_assignements.set_row(i,vec(gammas[i], ac.K));
     }	
}

/*void
accum_gammas (double **gammas, double *_gamma, int n, int K)
{*/
	/* gammas(n,:) = gamma for Newton-Raphson of alpha */
/*	int k;
	for (k = 0; k < K; k++)
		gammas[n][k] = _gamma[k];
	return;
}*/

void
accum_betas (double **betas, int K, vertex_data & data)
{
        graphlab::timer t; t.start();
	int i, k;
	int n = data.datapoint.nnz();

	for (i = 0; i < n; i++){
		int id = data.datapoint.get_nz_index(i);
                int cnt = (int)data.datapoint.get_nz_data(i);
		for (k = 0; k < ac.K; k++)
			betas[id][k] += data.distances[i*ac.K+k] * cnt;
        }
	ps.counter[LDA_ACCUM_BETA] += t.current_time();
}


double **
dmatrix (int rows, int cols)
{
        double **matrix;
        int i;

        matrix = (double **)calloc(rows, sizeof(double *));
        if (matrix == NULL)
                return NULL;
        for (i = 0; i < rows; i++) {
                matrix[i] = (double *)calloc(cols, sizeof(double));
                if (matrix[i] == NULL)
                        return NULL;
        }

        return matrix;
}

void
free_dmatrix (double **matrix, int rows)
{
        int i;
        for (i = 0; i < rows; i++)
                free(matrix[i]);
        free(matrix);
}



double
lda_lik (double **beta, double **gammas, int m)
{
	double **egammas;
	double z, lik;
	int i, j, k;
	int n;
	lik = 0;
	graphlab::timer t; t.start();
	
	
	if ((egammas = dmatrix(m, ac.K)) == NULL) {
		fprintf(stderr, "lda_likelihood:: cannot allocate egammas.\n");
		exit(1);
	}
	normalize_matrix_row(egammas, gammas, m, ac.K);
	
	for (i = 0; i < ps.M; i++)
	{
	        vertex_data &data = ps.g->vertex_data(i);
		n = data.datapoint.nnz();
		for (j = 0; j < n; j++) {
			int pos = data.datapoint.get_nz_index(j);
		        int cnt = data.datapoint.get_nz_data(j);
			for (k = 0, z = 0; k < ac.K; k++)
				z += beta[pos][k] * egammas[i][k];
			lik += cnt * log(z);
		}
	}

	free_dmatrix(egammas, m);
        ps.counter[LDA_LIKELIHOOD] += t.current_time();
	return lik;

}


void lda_em_update_function(gl_types::iscope & scope,
	                    gl_types::icallback & scheduler)
     /* const document *d, double *gamma, double **aq,
      double *nt, double *pnt, double *ap,
      const double *alpha, const double **beta,
      int L, int K, int ac.iter*/
{
	int j, k, l;
	double z;
	
	int thread_id = graphlab::thread::thread_id();
        scratch_buffer & myscratch = scratch[thread_id];

        vertex_data &data = scope.vertex_data();
        int L = data.datapoint.nnz();

        if (data.distances.size() == 0)
          data.distances = zeros(data.datapoint.nnz() * ac.K);

        double L_over_K = (double)L / ac.K;
	for (k = 0; k < ac.K; k++)
		myscratch.nt[k] = L_over_K;
	
	for (j = 0; j < ac.em_max_inner_iter; j++)
	{
		/* vb-estep */
		for (k = 0; k < ac.K; k++)
			myscratch.ap[k] = exp(psi(alpha[k] + myscratch.nt[k]));

		/* accumulate q */
		for (l = 0; l < L; l++){
	                int pos = data.datapoint.get_nz_index(l);
			for (k = 0; k < ac.K; k++)
				myscratch.q[l][k] = beta[pos][k] * myscratch.ap[k];
                }
		/* normalize q */
		for (l = 0; l < L; l++) {
			z = 0;
			for (k = 0; k < ac.K; k++)
				z += myscratch.q[l][k];
			for (k = 0; k < ac.K; k++)
				myscratch.q[l][k] /= z;
		}
		/* vb-mstep */
		for (k = 0; k < ac.K; k++) {
			z = 0;
			for (l = 0; l < L; l++){
				int cnt = (int)data.datapoint.get_nz_data(l);
				z += myscratch.q[l][k] * cnt;
                        }
			myscratch.nt[k] = z;
		}
		/* converge? */
		if ((j > 0) && converged(myscratch.nt, myscratch.pnt, ac.K, 1.0e-2))
			break;
		//for (k = 0; k < ac.K; k++)
		//	pnt[k] = myscratch.nt[k];
	        memcpy(myscratch.pnt, myscratch.nt, ac.K*sizeof(double));
        }

	for (k = 0; k < ac.K; k++)
	//	_gamma[k] = alpha[k] + nt[k];
        	gammas[scope.vertex()][k] = alpha[k] + myscratch.nt[k]; 


	for (int i = 0; i < L; i++){
		for (k = 0; k < ac.K; k++)
		    data.distances[i*ac.K+k] = myscratch.q[i][k];
			/*betas[id][k] += q[i][k] * cnt;*/
        }

	return;
}

int
converged (double *u, double *v, int n, double threshold)
{
        /* return 1 if |a - b|/|a| < threshold */
        double us = 0;
        double ds = 0;
        double d;
        int i;
        
        for (i = 0; i < n; i++)
                us += u[i] * u[i];

        for (i = 0; i < n; i++) {
                d = u[i] - v[i];
                ds += d * d;
        }

        if (sqrt(ds / us) < threshold)
                return 1;
        else
                return 0;

}

void
normalize_matrix_col (double **dst, double **src, int rows, int cols)
{
        /* column-wise normalize from src -> dst */
        double z;
        int i, j;

        for (j = 0; j < cols; j++) {
                for (i = 0, z = 0; i < rows; i++)
                        z += src[i][j];
                for (i = 0; i < rows; i++)
                        dst[i][j] = src[i][j] / z;
        }
}

void
normalize_matrix_row (double **dst, double **src, int rows, int cols)
{
        /* row-wise normalize from src -> dst */
        int i, j;
        double z;

        for (i = 0; i < rows; i++) {
                for (j = 0, z = 0; j < cols; j++)
                        z += src[i][j];
                for (j = 0; j < cols; j++)
                        dst[i][j] = src[i][j] / z;
        }
}

void
lda_main ()
{
        /* allocate parameters */
        if ((alpha = (double *)calloc(ps.K, sizeof(double))) == NULL) {
                fprintf(stderr, "lda:: cannot allocate alpha.\n");
                exit(1);
        }
        if ((beta = dmatrix(ps.N, ac.K)) == NULL) {
                fprintf(stderr, "lda:: cannot allocate beta.\n");
                exit(1);
        }
        
        lda_learn (alpha, beta);

        free_dmatrix(beta, ps.N);
        free(alpha);
}


void
newton_alpha (double *alpha, double **gammas, int M, int K, int level)
{
	int i, j, t;
	double *g, *h, *pg, *palpha;
	double z, sh, hgz;
	double psg, spg, gs;
	double alpha0, palpha0;
        graphlab::timer tt; tt.start();

	/* allocate arrays */
	if ((g = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate g.\n");
		return;
	}
	if ((h = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate h.\n");
		return;
	}
	if ((pg = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate pg.\n");
		return;
	}
	if ((palpha = (double *)calloc(K, sizeof(double))) == NULL) {
		fprintf(stderr, "newton:: cannot allocate palpha.\n");
		return;
	}

	/* initialize */
	if (level == 0)
	{
		double factor = M*K;
		for (i = 0; i < K; i++) {
			for (j = 0, z = 0; j < M; j++)
				z += gammas[j][i];
			alpha[i] = z / factor;
		}
	} else {
                double factor = M*K*pow(10,level);
		for (i = 0; i < K; i++) {
			for (j = 0, z = 0; j < M; j++)
				z += gammas[j][i];
			alpha[i] = z / factor;
		}
	}

	psg = 0;
	for (i = 0; i < M; i++) {
		for (j = 0, gs = 0; j < K; j++)
			gs += gammas[i][j];
		psg += psi(gs);
	}
	for (i = 0; i < K; i++) {
		for (j = 0, spg = 0; j < M; j++)
			spg += psi(gammas[j][i]);
		pg[i] = spg - psg;
	}

	/* main iteration */
	for (t = 0; t < MAX_NEWTON_ITERATION; t++)
	{
		for (i = 0, alpha0 = 0; i < K; i++)
			alpha0 += alpha[i];
		palpha0 = psi(alpha0);
		
		for (i = 0; i < K; i++)
			g[i] = M * (palpha0 - psi(alpha[i])) + pg[i];
		for (i = 0; i < K; i++)
			h[i] = - 1 / ppsi(alpha[i]);
		for (i = 0, sh = 0; i < K; i++)
			sh += h[i];
		
		for (i = 0, hgz = 0; i < K; i++)
			hgz += g[i] * h[i];
		hgz /= (1 / ppsi(alpha0) + sh);

		for (i = 0; i < K; i++)
			alpha[i] = alpha[i] - h[i] * (g[i] - hgz) / M;
		
		for (i = 0; i < K; i++)
			if (alpha[i] < 0) {
				if (level >= MAX_RECURSION_LIMIT) {
					fprintf(stderr, "newton:: maximum recursion limit reached.\n");
					exit(1);
				} else {
					free(g);
					free(h);
					free(pg);
					free(palpha);
					ps.counter[LDA_NEWTON_METHOD] += tt.current_time();
					return newton_alpha(alpha, gammas, M, K, 1 + level);
				}
			}

		if ((t > 0) && converged(alpha, palpha, K, 1.0e-4)) {
			free(g);
			free(h);
			free(pg);
			free(palpha);
			ps.counter[LDA_NEWTON_METHOD] += tt.current_time();
			return;
		} else
			for (i = 0; i < K; i++)
				palpha[i] = alpha[i];
		
	}
	fprintf(stderr, "newton:: maximum iteration reached. t = %d\n", t);
	ps.counter[LDA_NEWTON_METHOD] += tt.current_time();
	return;

}
