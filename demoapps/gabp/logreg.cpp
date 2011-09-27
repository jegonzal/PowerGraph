
/*
   Copyright [2011] [Aapo Kyrola, Joseph Bradley, Danny Bickson, Carlos Guestrin / Carnegie Mellon University]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF Aps.m KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
//
// Implementation of Shotgun Logreg - a parallel Logistic Lasso solver - for OpenMP.
// Optimization problem
//      \arg \max_x \sum log(1+exp(-yAx) - \shotgun_lambda |x|_1
//
// Based on Coordinate Descent Newton algorithm described in
// Yuan, Chang et al. : 
//      A Comparison of Optimization Methods and Software for Large-scale L1-regularized Linear Classification
//  
// \author Aapo Kyrola akyrola@cs.cmu.edu
// adapted to GraphLab by Danny Bickson, CMU
// 

 
//#include "common.h"
#include "linear.h"
#ifdef HAS_ITPP
#include "cas_array.h"
#include "advanced_config.h"
#include "itpp/itbase.h"
#include <assert.h>

// Parameters
extern problem_setup ps;
extern advanced_config config;
extern graph_type_shotgun *g;

cas_array<double> Gmax;
double Gmax_old;
double Gmax_init;
bool shuffle = true;

double logloss(double x) {
     if (x > (-10) && x < 10) {
        return log(1 + exp(x));
     } else if (x<= (-10)) {
        return 0.0;
     } else {
        return x;
     }
 }

double compute_llhood() {
     double llhood = 0;
     for(size_t i= 0; i<ps.m; i++) {
        vertex_data_shotgun & vdata = g->vertex_data(i);
        itpp::sparse_vec& row = vdata.features;
        double Ax=0;
        for(int j=0; j<row.nnz(); j++) {
            assert(row.get_nz_index(j) < (int)ps.n);
            Ax += g->vertex_data(ps.m+row.get_nz_index(j)).x*row.get_nz_data(j);
        }
        llhood -= logloss(-vdata.y*Ax);
     }
     return llhood;
 }

 double compute_objective_logreg(double lambda, double * l1x = NULL, double * loglikelihood = NULL, int * _l0 = NULL,
                double * testobj = NULL) {

    double penalty = 0.0;
    int l0 = 0;
    for(size_t j=ps.m; j<(size_t)ps.last_node; j++) {
        penalty += std::fabs(g->vertex_data(j).x);
        l0 += (g->vertex_data(j).x != 0);
    }
    double llhood = compute_llhood();
    double llhood_test = 0;// //TODO compute_llhood(logregprob->x, logregprob->ytest, logregprob->Atest_rows);
    if (l1x != NULL) *l1x = penalty;
    if (loglikelihood != NULL) *loglikelihood = llhood;
    if (_l0 != NULL) *_l0 = l0;
    if (testobj != NULL) *testobj = -penalty*lambda + llhood_test; 
    return -penalty*lambda + llhood;
}



inline double sign(double a) {
    return (a<0 ? -1 : (a == 0 ? 0 : 1));
}
inline void swap(graphlab::vertex_id_t &t1, graphlab::vertex_id_t &t2) { int tmp=t2; t2=t1; t1=tmp; }
  
   
  

// Computes L_j'(0) and L_j''
// See: Yuan, Chang et al. : 
//  A Comparison of Optimization Methods and Software for Large-scale L1-regularized Linear Classification
//  equation (25)
//inline void logreg_cdn_derivandH(double shotgun_lambda, int x_i, double &G, double& H) {
void logreg_cdn_derivandH(vertex_data_shotgun& vdata, double &G, double &H){
    G = H = 0;
    itpp::sparse_vec& col = vdata.features;
    assert(col.nnz()>0);
    for(int i=0; i<col.nnz(); i++) {
       int rowi = col.get_nz_index(i);
       assert(rowi>=0 && rowi < (int)ps.m); 
       double val = col.get_nz_data(i);
       double exp_wTxind = g->vertex_data(rowi).expAx;
	double tmp1 = val/(1+exp_wTxind);
	double tmp2 =  tmp1;
	double tmp3 = tmp2*exp_wTxind;
	G += tmp2;
	H += tmp1*tmp3;
    }
    G = -G + vdata.xjneg;
    G /= config.shotgun_lambda;
    H /= config.shotgun_lambda;
    if (H<1e-5) {
        H = 1e-5;
    }
 }
 


 // L_j(x + diff*e_j)-L_j(x) (see function g_xi())
 // (eq. 18)
//inline double logreg_cdn_Ldiff(double shotgun_lambda, int x_i, double diff) {
double logreg_cdn_Ldiff(vertex_data_shotgun & vdata, double diff){
    double sum = 0.0;
    itpp::sparse_vec& col = vdata.features;
    for(int i=0; i<col.nnz(); i++) {
       int rowi = col.get_nz_index(i);
       assert(rowi >= 0 && rowi < (int)ps.m);
       double dc = diff * col.get_nz_data(i);
       double expdiff = exp(dc);
       double expAXdiff = g->vertex_data(rowi).expAx * expdiff;
       assert(expAXdiff + expdiff != 0);
       double ds = log(expAXdiff + 1.0) - log(expAXdiff+expdiff);
       sum +=  ds;
    } 
    if (std::isnan(sum)){
        fprintf(stderr, "Got numerical error: please verify that in your dataset there are no columns of matrix A with all zeros.\n");
        exit(1);
    }
    return 1.0/config.shotgun_lambda * (diff*vdata.xjneg + sum);
}
 
 // Equation 17. One-variable function that equals the change in loss function
 // when x_i is change by z
//double g_xi(double z, int x_i, double shotgun_lambda) {
double g_xi(double z, vertex_data_shotgun & vdata){
    double xv = vdata.x;
    return std::fabs(xv+z) - std::fabs(xv) + logreg_cdn_Ldiff(vdata, z);
}
 

// Compute A'A in parallel. Note: could also compute in demand?
void initialize_all() {
    //logregprob->x.resize(ps.n);
    //logregprob->Ax.resize(ps.m);
    //logregprob->expAx.resize(ps.m);
    Gmax.resize(1);

    //active = (bool *) calloc(ps.n,sizeof(bool));
    //xjneg = (double *) calloc(ps.n,sizeof(double));
    //ps.cdn_pos_y = 0;
    for(int i=0; i<(int)ps.m; i++)  {
        //logregprob->expAx[i] = 1.0; // since(exp(0) = 1)
        ps.cdn_pos_y += (g->vertex_data(i).y == 1);
    }
    ps.cdn_neg_y = ps.m-ps.cdn_pos_y;
    
    assert(ps.cdn_pos_y > 0 && ps.cdn_neg_y > 0);
    //for(int i=0; i<ps.n; i++) 
    //	active[i] = true;
    
    // use the trick that log(1+exp(x))=log(1+exp(-x))+x
    //#pragma omp parallel for
    for(int i=ps.m; i<(int)(ps.m+ps.n); i++) {
        vertex_data_shotgun &vdata = g->vertex_data(i);
        itpp::sparse_vec& col = vdata.features;
        for(int j=0; j<col.nnz(); j++) { 
             assert(col.get_nz_index(j) < ps.m);
            if (g->vertex_data(col.get_nz_index(j)).y == -1)
                vdata.xjneg += col.get_nz_data(j);
        }
    }
}

void recompute_expAx() {
    graphlab::timer t; t.start();

    //#pragma omp parallel for
    for(int i=0; i<(int)ps.m; i++) {
        double Ax=0;
        vertex_data_shotgun & vdata = g->vertex_data(i);
        itpp::sparse_vec &row = vdata.features;
        for(int j=0; j<row.nnz(); j++) {
            assert(row.get_nz_index(j) < (int)ps.n);
            Ax += g->vertex_data(ps.m+row.get_nz_index(j)).x*row.get_nz_data(j);
        }
        vdata.expAx = exp(Ax);
    }
    ps.counter[RECOMPUTE_EXP_AX_LOGREG] += t.current_time();
}
 
void recompute_partial_expAx(vertex_data_shotgun & vdata, double newval, double oldval) {
    graphlab::timer t; t.start();

    //#pragma omp parallel for
    //for(int i=0; i<(int)ps.m; i++) {
        //double Ax=0;
        //vertex_data_shotgun & vdata = g->vertex_data(i);
        //itpp::sparse_vec &row = vdata.features;
        //for(int j=0; j<row.nnz(); j++) {
        //    assert(row.get_nz_index(j) < (int)ps.n);
        //    Ax += g->vertex_data(ps.m+row.get_nz_index(j)).x*row.get_nz_data(j);
        //}
        //vdata.expAx = exp(Ax);
    //}
    itpp::sparse_vec &col = vdata.features;
    for (int j=0; j< col.nnz(); j++){
       double A_ij = col.get_nz_data(j);
       int rowi = col.get_nz_index(j);
       vertex_data_shotgun & row = g->vertex_data(rowi);
       //row.expAx = (row.expAx / exp(oldval*A_ij)) * exp(newval*A_ij);
       double fact = exp((newval - oldval)*A_ij);
       mul(&row.expAx, fact); 
    }
    ps.counter[RECOMPUTE_EXP_AX_LOGREG] += t.current_time();
}
 



void shotgun_logreg_update_function(gl_types_shotgun::iscope & scope,
				    gl_types_shotgun::icallback & scheduler){

// Yua, Chang, Hsieh and Lin: A Comparison of Optimization Methods and Software for Large-scale L1-regularized Linear Classification; p. 14
//double shoot_cdn(int x_i, double shotgun_lambda) {
    // Compute d: (equation 29), i.e the solution to the quadratic approximation of the function 
    // for weight x_i
    assert(scope.vertex() >= (graphlab::vertex_id_t)ps.m && scope.vertex() < (graphlab::vertex_id_t)ps.last_node); 
    vertex_data_shotgun &vdata = scope.vertex_data();

    if (!vdata.active || vdata.features.nnz() == 0){ 
      vdata.x = 0; 
      return;
    }
   
    double oldval = vdata.x;
    double violation = 0.0;
    double xv = vdata.x;
    
    double Ld0, Ldd0;
    //logreg_cdn_derivandH(shotgun_lambda, x_i, Ld0, Ldd0);    
    logreg_cdn_derivandH(vdata, Ld0, Ldd0);    

    double Gp = (Ld0+1);
    double Gn = (Ld0-1);

    // Shrinking (practically copied from LibLinear)
    if (xv == 0) {
        if (Gp<0) {
            violation = -Gp;   
        } else if (Gn > 0) {
            violation = Gn;
       	} else if(Gp>Gmax_old/ps.m && Gn<-Gmax_old/ps.m) {
            // Remove
            vdata.active = false;
            return;
        }   
    } else if(xv > 0)
      violation = fabs(Gp);
    else
      violation = fabs(Gn);
    
    Gmax.max(0, violation);
    if (config.debug)
    printf("node %d violation %g Ld0 %g Ldd0 %g Gp %g Gn %g xv %g \n", scope.vertex()-ps.m, violation, Ld0, Ldd0, Gp, Gn, xv);

    // Newton direction d
    double rhs = Ldd0*xv;
    double d;
    if (Gp<= rhs) {
        d = -(Gp/Ldd0);
    } else if (Gn >= rhs) {
        d = -(Gn/Ldd0); 
    } else {
        d = -xv;
    }
    
     if (std::fabs(d) < 1e-12) {
        return;
    }
    // Small optimization
    d = std::min(std::max(d,-10.0),10.0);
   
    // Backtracking line search (with Armijo-type of condition) (eq. 30)
    int iter=0;
    double gamma = 1.0; // Note: Yuan, Chang et al. use shotgun_lambda instead of gamma
    double delta = (Ld0 * d + std::fabs(xv+d) - std::fabs(xv));
    double rhs_c = config.shotgun_sigma * delta;    
    do {
        //double change_in_obj = g_xi(d, x_i, shotgun_lambda);
        double change_in_obj = g_xi(d, vdata);
        if (change_in_obj <= gamma * rhs_c) {
            // Found ok.
            vdata.x += d;
            // Update dot products (Ax)
            itpp::sparse_vec &col = vdata.features;
            //#pragma omp parallel for
            for(int i=0; i<col.nnz(); i++) {
                //logregprob->expAx.mul(col.idxs[i], exp(d * col.values[i]));
                assert(col.get_nz_index(i) < ps.m);
                g->vertex_data(col.get_nz_index(i)).expAx*= exp(d* col.get_nz_data(i));
                //TODO
            }
            return;
        }
        gamma *= 0.5;
        d *= 0.5;
    } while(++iter < config.shotgun_max_linesearch_iter); 
    recompute_partial_expAx(vdata, vdata.x, oldval);
    //vdata.x = 0;
    return;
}

 
  
void calc_and_print_objective(){ 
 double l1x=0, loglikelihood=0;
 int l0=0;
 double obj = compute_objective_logreg(config.shotgun_lambda, &l1x, &loglikelihood, &l0, NULL);
 if (l1x == 0)	
   ps.cdn_all_zero = true;
 printf("objective is: %g l1: %g loglikelihood %g l0: %d\n", obj, l1x, loglikelihood, l0); 
  ps.last_cost = obj;
}


/**
  * Main optimization loop.
  * Note: this version of code does NOT have special version for sequential
  * version. Sequential version is slightly lighter than parallel version because
  * it can use a more efficient operation to maintain the active set. In practice, this
  * has little effect. Also, it does not need to have a atomic array for maintaining Ax.
  * For the experiments in the paper, a special sequential code was used for fairness.
  */
//void compute_logreg(shotgun_data * prob, double shotgun_lambda, double config.threshold, int config.iter, int config.display_cost, bool & ps.cdn_all_zero) {
void compute_logreg(gl_types_shotgun::core & glcore){ 

    g = &glcore.graph();   
    ps.cdn_all_zero = false;
    //logregprob = prob;
    //double l1x, loglikelihood;
    ps.iiter = 0;//, t=0;
    ps.shotgun_numshoots = 0;

    std::vector<graphlab::vertex_id_t> shuffled_indices(ps.n);
    for (graphlab::vertex_id_t j=ps.m; j<(graphlab::vertex_id_t)ps.last_node; j++) 
	shuffled_indices[j-ps.m] = j;

    Gmax_old = 1e30;
    // Adjust threshold similarly as liblinear
    initialize_all();
    config.threshold =  config.threshold*std::min(ps.cdn_pos_y,ps.cdn_neg_y)/double(ps.m);
    
    while(true) {
        size_t active_size = ps.n;
        ps.shotgun_numshoots += active_size;   
        
         // Randomization
        if (shuffle)
        for(size_t j=0; j<active_size; j++) {
                graphlab::vertex_id_t i = j+rand()%(active_size-j);
                swap(shuffled_indices[i], shuffled_indices[j]);
        } 
          
 	glcore.add_tasks(shuffled_indices, shotgun_logreg_update_function,1);
        glcore.start(); 
            /* Main parallel loop */
        //#pragma omp parallel for
        //for(int s=0; s<active_size; s++) {
        //    int x_i = shuffled_indices[s];
        //    shoot_logreg(x_i, shotgun_lambda);
        //}
            
        /* Gmax handling */
        Gmax_old = Gmax[0];
        if (ps.iiter == 0) {
            Gmax_init = Gmax_old;
        }
        
        ps.iiter++;
        
        //std::cout << Gmax.get_value() <<  " " << Gmax_init << " " <<  config.threshold*Gmax_init << std::endl;
        if (ps.iiter > config.iter && config.iter>0) {
            printf("Exceeded max ps.iiter: %d\n", config.iter);
            break;
        }
        
        for(graphlab::vertex_id_t i=ps.m; i<(graphlab::vertex_id_t)ps.last_node; i++) 
          shuffled_indices[i-ps.m] = i;
        
        active_size = ps.n;
        for(size_t s=0; s<active_size; s++) {
            int j = shuffled_indices[s];
            if (!g->vertex_data(j).active) {
                active_size--;
                swap(shuffled_indices[s], shuffled_indices[active_size]);
                s--;
            }
        }
        
        if (Gmax[0] <= config.threshold*Gmax_init) {
           // std::cout << active_size << std::endl;
            if (active_size == ps.n) {
                printf("Encountered all zero solution! try to decrease shotgun_lambda\n");
 		ps.cdn_all_zero = true;
   		break;              
            } else {
                 Gmax_old = 1e30;
                 for(size_t i=ps.m; i<(size_t)ps.last_node; i++) 
                   g->vertex_data(i).active = true;
                 active_size=ps.n;
                 //recompute_partial_expAx(vdata, vdata.x, oldval);
                 //continue;
            }
        }  

     if (config.display_cost){
        calc_and_print_objective();   
     }
   }// end ps.iiter

   if (!config.display_cost){
     calc_and_print_objective();
   }


   //delete[] active;
  // delete[] xjneg;
  printf("Finished Shotgun CDN in %d ps.iiter\n", ps.iiter);
}

#endif
