/*
   Copyright [2011] [Aapo Kyrola, Joseph Bradley, Danny Bickson, Carlos Guestrin / Carnegie Mellon University]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF Alassoprob->ny KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
// Optimization problem
//      \arg \min_x ||Ax-y|| + \lambda |x|_1
//

//#include "common.h"
#include "linear.h"
#include "advanced_config.h"

// Problem definition
//shotgun_data * lassoprob;
graph_type_shotgun* g;
extern advanced_config config;
extern problem_setup ps;

// Major optimization: always keep updated vector 'Ax'
 
void initialize_feature(vertex_data_shotgun & vdata) {
    itpp::sparse_vec& col = vdata.features;//lassoprob->A_cols[feat_idx];
    //feature& feat = lassoprob->feature_consts[feat_idx];
    //lassoprob->x[feat_idx] = 0.0;

    // Precompute covariance of a feature
    //feat.covar = 0;
    // Precompute (Ay)_i
   vdata.Ax = 0.0; 
   vdata.y = 0.0;
   for(int i=0; i<col.nnz(); i++) {
        vdata.Ax += powf(col.get_nz_data(i),2);
        vdata.y += col.get_nz_data(i)*g->vertex_data(col.get_nz_index(i)).y;
    }
    vdata.Ax *= 2;
    vdata.y *= 2;
    

}

void initialize() {
    //lassoprob->feature_consts.reserve(lassoprob->nx);
    //lassoprob->x.resize(lassoprob->nx);
    //lassoprob->Ax.resize(lassoprob->ny);
    //#pragma omp for
    for(int i=ps.m; i<ps.last_node; i++) {
        initialize_feature(g->vertex_data(i));
    }
}

double soft_thresholdO(double _lambda, double shootDiff) {
    return (shootDiff > _lambda)* (_lambda - shootDiff) + 
	               (shootDiff < -_lambda) * (-_lambda- shootDiff) ;
}

double soft_threshold(double _lambda, double shootDiff) {
  if (shootDiff > _lambda) return _lambda - shootDiff;
  if (shootDiff < -_lambda) return -_lambda - shootDiff ;
  else return 0;
}


void shotgun_lasso_update_function(gl_types_shotgun::iscope & scope,
	                           gl_types_shotgun::icallback& scheduler){
//double shoot(int x_i, double lambda) {
    //feature& feat = lassoprob->feature_consts[x_i];
    //double oldvalue = lassoprob->x[x_i];
    vertex_data_shotgun & vdata = scope.vertex_data();
    double oldvalue = vdata.x;    

    // Compute dotproduct A'_i*(Ax)
    double AtAxj = 0.0;
    itpp::sparse_vec& col = vdata.features;//lassoprob->A_cols[x_i];
    int len=col.nnz();
    for(int i=0; i<len; i++) {
        AtAxj += col.get_nz_data(i) * g->vertex_data(col.get_nz_index(i)).Ax;
    }
    
    double S_j = 2 * AtAxj - vdata.Ax * oldvalue - vdata.y;//feat.covar * oldvalue - feat.Ay_i;
    double newvalue = soft_threshold(config.shotgun_lambda,S_j)/vdata.Ax; //feat.covar;
    double delta = (newvalue - oldvalue);
    
    // Update ax
    if (delta != 0.0) {
        for(int i=0; i<len; i++) {
	       //lassoprob->Ax.add(col.idxs[i], col.values[i] * delta);
	       g->vertex_data(col.get_nz_index(i)).Ax += col.get_nz_data(i)*delta;
        }
        
        //lassoprob->x[x_i] = newvalue;
        vdata.x = newvalue;
    }
    //return std::abs(delta);
    vdata.expAx = std::fabs(delta);
}

// Find such lambda that if used for regularization,
// optimum would have all weights zero.
double compute_max_lambda() {
    double maxlambda = 0.0;
    for(int i=ps.m; i<ps.last_node; i++) {
        maxlambda = std::max(maxlambda, std::fabs(g->vertex_data(i).y)); //lassoprob->feature_consts[i].Ay_i));
    }
    return maxlambda;
}

double get_term_threshold(int k, int K, double delta_threshold) {
  // Stricter termination threshold for last step in the optimization.
  return (k == 0 ? delta_threshold  : (delta_threshold + k*(delta_threshold*50)/K));
}

 double compute_objective(double _lambda, double & l0x, double * l1x = NULL, double * l2err = NULL ) {

    double least_sqr = 0;
     
    //for(int i=0; i<lassoprob->ny; i++) {
     for (int i=0; i<(int)ps.m; i++){   
    // least_sqr += (lassoprob->Ax[i]-lassoprob->y[i])*(lassoprob->Ax[i]-lassoprob->y[i]);					
      vertex_data_shotgun & vdata = g->vertex_data(i);
      least_sqr += powf(vdata.Ax - vdata.y,2);
    }
    
     // Penalty check for feature 0
    double penalty = 0.0;
    //for(int i=0; i<lassoprob->nx; i++) {
    for (int i=ps.m; i< ps.last_node; i++){ 
       vertex_data_shotgun & vdata = g->vertex_data(i);
       //penalty += std::abs(lassoprob->x[i]);
       penalty += std::fabs(vdata.x);
       l0x += (vdata.x == 0.0);
    }
    if (l1x != NULL) *l1x = penalty;
    if (l2err != NULL) *l2err = least_sqr;
    return penalty * _lambda + least_sqr;
}



//void main_optimization_loop(double lambda, int regpathlength, double threshold, int maxiter, int verbose) {
void main_optimization_loop(gl_types_shotgun::core &glcore){
    // Empirically found heuristic for deciding how malassoprob->ny steps
    // to take on the regularization path
    int regularization_path_length = (config.shotgun_reg_path_len <= 0 ? 1+(int)(ps.n/2000) : config.shotgun_reg_path_len);
    printf("regularization_path_length = %d\n",regularization_path_length);

    ps.gt.start();   
 
    double lambda_max = compute_max_lambda();
    double lambda_min = config.shotgun_lambda;
    double alpha = pow(lambda_max/lambda_min, 1.0/(1.0*regularization_path_length));
    int regularization_path_step = regularization_path_length;

    double delta_threshold = config.threshold;
    long long int num_of_shoots = 0;
    int counter = 0;
    int iterations = 0;
    //double *delta = new double[ps.n];

    std::vector<graphlab::vertex_id_t> indices(ps.n);
    for (graphlab::vertex_id_t j=ps.m; j<(graphlab::vertex_id_t)ps.last_node; j++) 
	indices[j-ps.m] = j;


    do {
        ++iterations;
        //if (iterations >= maxiter && maxiter > 0) {
        if (iterations>= config.iter && config.iter > 0){ 
           printf("Exceeded max iterations: %d", config.iter);
            return;
        }
        double maxChange=0;
        config.shotgun_lambda = lambda_min * pow(alpha, regularization_path_step);
        
        // Parallel loop
        //#pragma omp parallel for  
       	glcore.add_tasks(indices, shotgun_lasso_update_function,1);
        glcore.start(); 

        //for(int i=ps.m; i<ps.last_node; i++) 
        //    delta[i] = shoot(i, lambda);
        maxChange = 0;
        // TODO: use OpenMP reductions (although - not much to gain)
        //for(int i=0; i<lassoprob->nx; i++)
        for (int i=ps.m ; i< ps.last_node; i++){  
          double delta = g->vertex_data(i).expAx;
          maxChange = (maxChange < delta ? delta : maxChange);
        }
        //num_of_shoots += lassoprob->nx;
        num_of_shoots += ps.n;
        counter++;

        // Convergence check.
        // We use a simple trick to converge faster for the intermediate sub-optimization problems
        // on the regularization path. This is important because we do not care about accuracy for
        // the intermediate problems, just want a good warm start for next round.
        bool converged = (maxChange <= get_term_threshold(regularization_path_step,regularization_path_length,delta_threshold));
        if (converged || counter>std::min(100, (100-regularization_path_step)*2)) {
            counter = 0;
            regularization_path_step--; 
        }
        if (config.display_cost){
          double l1x = 0, l2err = 0, l0x = 0;
          double obj = compute_objective(config.shotgun_lambda, l0x, &l1x, &l2err);
          printf("%g) Objective: %g L1 : %g L2err: %g L0: %d\n", ps.gt.current_time(), obj, l1x, l2err, (int)l0x);
          ps.last_cost = obj;
        }
    } while (regularization_path_step >= 0);
    //delete[] delta;
    printf("Num of shoots = %lld\n", num_of_shoots);
}


//double solveLasso(shotgun_data  * probdef, double lambda, int regpathlength, double threshold, int maxiter, int verbose) {
void solveLasso(gl_types_shotgun::core &  glcore){
    g = &glcore.graph();
    initialize();
    //main_optimization_loop(lambda, regpathlength, threshold, maxiter, verbose); 
    main_optimization_loop(glcore);
    return;
}


