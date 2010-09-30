#ifndef aapolassoAPP_HPP
#define aapolassoAPP_HPP

/*
 * Lasso Shooting implementation
 * \author akyrola
 */
#include <cmath>
#include <cstdio>
#include <string>


#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/graphlab/logger/assertions.hpp>

using namespace graphlab;
typedef graphlab::blob_types gl_types;



/**
 * Solves min_x ||Ax-y||^2 + lambda * ||x||_1 
 *     =  min_x  x'A'Ax - 2y'Ax + y'y + lambda * ||x||_1
 *     =  min_x  x'Cx - dx + lambda * ||x||_1
 *   where C = A'A and d = 2y'A are precomputed
 */

/** Vertex and edge data types **/

/*struct idx_value {
  int idx;
  double value;
  idx_value(int _idx, double _value) :idx(_idx), value(_value) {}
};*/	

struct vertex_data {
  double value;
  double d_i;
  double C_ii;

  int * C_j_indices;
  double * C_j_vals;
	
  int C_j_size;
  int C_j_counter;
  spinlock lock;
  
  // Relative penalty on the lambda for this weight
  double lambda_cofactor;
}; 


struct proc_vertex_data {
  int process_id;
  double * workarray;
  double * workarray2;
  
  std::vector<int> * myvars;
  int num_of_iters;
};

// Runnign modes
enum { COMPRESSED_SENSING = 1, ACCURATE = 0 };
int mode, convergence_mode;


// TODO!!!
double DESIRED_ACCURACY = 0.01;
int MAXZEROSWEEPS = 5;

struct proc_to_x_edge_data {
  bool my_vertex;
};

struct sparse_array_s {
  int size;
  int capacity;
  int * idxs;
  double * values;
};

//// Need to have graph and engine as global variables to access
//// from distributed handlers.
blob_graph * g;
gl_types::iengine * engine;

bool engine_stopped;

bool memorySavingMode = false;

long int compacted = 0;
long int AAsize = 0;
long int Asize = 0;

int max_endgame_iterations = 500;

unsigned int num_of_processing_nodes;

// Maintained automatically
double MAX_VALUE = 0;

atomic<int> convcounter;
bool enable_locking = false, automatically_initiated_locking = false;

#define PROC_NODE_FOR_X(x_i) (g->num_vertices()-num_of_processing_nodes+x_i%num_of_processing_nodes)

// Distributed version
#define NUMX_PERPROCESS(procid, numproc) ((nx)/numproc)
#define FIRST_X_OF_PROCESS(procid, numproc) (NUMX_PERPROCESS(procid,numproc) * procid)
#define LAST_X_OF_PROCESS(procid, numproc) ((procid == (numproc-1)) ? nx-1 : NUMX_PERPROCESS(procid,numproc) * procid + NUMX_PERPROCESS(procid,numproc) - 1)

gl_types::vertex_id_t first_proc_node;
gl_types::vertex_id_t last_proc_node;

// X-vector. Use instead of vertex to speedup...
double * XV;
int * set_to_zero;

unsigned int ny;
	
bool * running;
bool * todo;  // Keeps list of x's that are to be done - used for task stealing

std::vector<double> y;

double lambda = 0.5; // TODO: Make shared variable

// Regularization path
double lambda_max =  200, lambda_min = 0.5, alpha;
int K = 10, k;

FILE * progressFile = NULL;
size_t endgame_iterations = 0;

typedef sparse_array_s sparse_array;

std::vector<sparse_array> * Acols;	
std::vector<sparse_array> * Arows;	

graphlab::distributed_control * dc;
unsigned int procid, numprocs;

int * array_sizes;

mutex safelock;

void expand(vertex_data * var_data) {
  	int * newarray =  (int *) malloc(sizeof(int) * var_data->C_j_size * 2);
  	memcpy(newarray, var_data->C_j_indices, var_data->C_j_size* sizeof(int));
  	free(var_data->C_j_indices);
  	var_data->C_j_indices = newarray;	
 
  	double * newarrayd =  (double *) malloc(sizeof(double) * var_data->C_j_size * 2);
 	memcpy(newarrayd, var_data->C_j_vals, var_data->C_j_size* sizeof(double));
  	free(var_data->C_j_vals);
  	var_data->C_j_vals = newarrayd;
    var_data->C_j_size *= 2;
}

spinlock lock;
mutex sendlock;

void compact(vertex_data * var_data) {
  	int * newarray =  (int *) malloc(sizeof(int) * var_data->C_j_counter);
  	double * newarrayd =  (double *) malloc(sizeof(double) * var_data->C_j_counter);
	memcpy(newarray, var_data->C_j_indices, var_data->C_j_counter* sizeof(int));
  	memcpy(newarrayd, var_data->C_j_indices, var_data->C_j_counter* sizeof(double));
  	free(var_data->C_j_indices);
  	free(var_data->C_j_vals);
  	var_data->C_j_vals = newarrayd;
  	var_data->C_j_indices = newarray;	
}


void aapolasso_update_function(gl_types::iscope &scope,
                               gl_types::icallback &scheduler,
                               gl_types::ishared_data* shared_data);
                               
                               
                          
 bool terminated = false;
 
 
 double calcx_1(double * x, int nx) {
 	double s = 0;
 	for(int i=0; i<nx; i++) {
 		s += std::abs(x[i]);
 	}
 	return s;
 }
                               
                               
 double compute_objective(double * x, int nx, int * _n_nonzeros, double * __lambda, double * _least_sqr, double * _penalty) {
        double penalty = 0;
        double least_sqr = 0;
        double _lambda = lambda;
       
        
        // Ax - y
        int ny = y.size();
        for(int i=0; i<ny; i++) {
            sparse_array row = (*Arows)[i];
            double s = 0;
            for(int j=0; j<row.size; j++) {
                s += row.values[j] * x[row.idxs[j]];
            }
            least_sqr += (s-y[i])*(s-y[i]); 	   			
        }
        
        int n_nonzeros = 0;
        
        /* NOTE: penalty calc correct only if only weight for index 0 is different from 1! */
  
        // Penalty
        penalty =  std::abs(x[0]) * g->vertex_data(0).as_ptr<vertex_data>()->lambda_cofactor;
        
        for(int i=1; i<nx; i++) {
            penalty += std::abs(x[i]);
            n_nonzeros += (x[i] != 0);
        }
        
        penalty *= _lambda;
        
        *_least_sqr = least_sqr;
        *_penalty = penalty;
        *__lambda = _lambda;
        *_n_nonzeros = n_nonzeros;
        return penalty + least_sqr;
}
                               
 /***
    * OBJECTIVE COMPUTATION
    */
 class objective_computation : public runnable {      
    public:

   double * lastx; 
   double lastobjective;
   double objective_when_diverged;
   

   
   void run() {
      	lastx = NULL;
      	lastobjective = -1;
    	timer runtimer;
        runtimer.start();
        
        fprintf(progressFile, "iter,time,nx,nonzeros,obj,least_sqr,penalty,target_obj,lambda,mode\n");
   
      	printf("---- Starting objective counter ----\n");
          do { 
             
            timer t;
            t.start();
            
            int nx = Acols->size();
            double * x = (double *) malloc(sizeof(double) * nx);
            memcpy(x, XV, sizeof(double) * nx);
        
            int n_nonzeros;
            double least_sqr = 0;
 			double _lambda, penalty;
 		    double obj = compute_objective(x, nx, &n_nonzeros, &_lambda, &least_sqr, &penalty);
 	   		
 	   		printf("OBJECTIVE: %lf (%lf,%lf) LAMBDA %lf NON-ZEROS: %d (%lf sec),%lf\n", obj, least_sqr, penalty, _lambda, n_nonzeros,
 	   		            t.current_time(), runtimer.current_time()); 
 	   		
 	   		/* Write progress */
 	   		
 	   		// 1. Count iterations
 	   		size_t iterations = endgame_iterations;
 	   		for(unsigned int vid = first_proc_node; vid<last_proc_node; vid++) {
			   iterations += (g->vertex_data(vid).as_ptr<proc_vertex_data>())->num_of_iters;
			}
			fprintf(progressFile, "%ld,%lf,%d,%d,%lf,%lf,%lf,%lf,%lf,%d\n", iterations, runtimer.current_time(),
								nx, n_nonzeros, obj, least_sqr, penalty, least_sqr + penalty/_lambda * lambda_min, _lambda, mode);
			fflush(progressFile);
 	   		
 	   		if (automatically_initiated_locking) {
 	   			if (obj < 0.98 * objective_when_diverged) {
 	   				// Try disabling locking
 	   				enable_locking = false;
 	   				automatically_initiated_locking = false;
 	   				printf(" >>>>>>>> REMOVED LOCKING");
 	   			}
 	   		}
 	   		
 	   		/* Divergence check - remember previous values */
 	   		if ((lastobjective > 0) && (obj > lastobjective * 2.0)) {
 				printf("=============== DIVERGENCE DETECTED ----- INSTALL PREVIOUS VALUE ==========\n"); 
 				/*printf("========= ENABLE LOCKING =========\n");
 				enable_locking = true;
 				automatically_initiated_locking = true;
 				memcpy(XV, lastx, sizeof(double)*nx);
 				free(x);
 				sleep(2);
 				objective_when_diverged = lastobjective;
 				continue;*/
 	   		}
 	   		
 	   		if (lastx != NULL) free(lastx);
 	   		lastx = x;
 	   		
 	   		lastobjective = obj;
 	   		
 	   		// Sleep for n secs
 	   		sleep(5);
 	      } while (!terminated);
 	      if (lastx != NULL) free(lastx);
 	  }
  };
 

/**  ==================================================
  *  ============== DISTRIBUTED STUFF ==================
  *  ==================================================
  */
  
  
  
 struct xvalue {
    int x_i;
    double value;
    xvalue(int i, double v) : x_i(i), value(v) {}
 };
  
 void receive_x_value(distributed_control& dc, size_t source, void* ptr, size_t len) {
 	 printf("Received remove XV\n");
     xvalue val = * ((xvalue *) ptr);
     XV[val.x_i] = val.value;
  }
  
 void send_x_value(int x_i) {
    xvalue xval = xvalue(x_i, XV[x_i]);
    for(unsigned int i=0; i<numprocs; i++) {
        if (i != procid) {
            dc->remote_call(i, receive_x_value, (void*) (&xval), sizeof(xvalue));
        }
    }
 }
 
 
 /// Remote call to schedule all vertices
 void schedule_all(distributed_control& dc, size_t source, void* ptr, size_t len) {
    printf("Remote scheduling request received from %ld\n", source);
    size_t n = g->num_vertices();
    for(size_t vid=n-num_of_processing_nodes; vid<n; vid++) {
      engine->get_scheduler().add_task(gl_types::update_task(vid, aapolasso_update_function), 1.0);
    }
   
 }
  
 void schedule_remotes() {
      for(unsigned int i=0; i<numprocs; i++) {
        if (i != procid) {
            dc->remote_call(i, schedule_all, NULL, 0);
        }
    }
 }
 
 void remote_set_lambda(distributed_control& dc, size_t source, void *ptr, size_t len, handlerarg_t _k) {
    lambda = *((double *) ptr);     
    k = _k;
    printf("Remote lambda change %lf\n", lambda);
    ASSERT_EQ(source, 0);
 }
 
 void set_lambda(double newlambda, int k) {
    ASSERT_EQ(procid, 0);
    for(unsigned int i=0; i<numprocs; i++) {
        if (i != procid) {
            dc->remote_call(i, remote_set_lambda, (void*) (&newlambda) , sizeof(double), k);
        } 
    }
    lambda = newlambda;
  }
  
  /***
    * OBJECTIVE COMPUTATION
    */
 class value_broadcaster : public runnable {      
    public:

   double * lastx; 


   void run() {
   	   	int nx = Acols->size();
   	   	lastx = (double *) malloc(nx * sizeof(double));
 	   	memcpy(lastx, XV, sizeof(double) * nx);
      	printf("---- Starting distributed broadcaster ----\n");
          do { 
            sendlock.lock();
 	   		for(int i=0; i<nx; i++) {
 	   		    double x = XV[i];
 	   		    if (x != lastx[i]) send_x_value(i);
 	   		    lastx[i] = x;
 	   		}
 	   		sendlock.unlock();
 	   	    // Sleep for 20 microsecs
 	   		usleep(20);
 	      } while (!terminated);
 	      free(lastx);
 	  }
  };
  
  
/**  ==================================================
  *  ============== COMPUTATION ==================
  *  ==================================================
  */  
void initialize_value(int x_i, vertex_data * var_data, double * workarray, double * workarray2);  
  
void computemem() {
    size_t memd=0, memcompactd=0, memf=0, memcompactf=0;
    
    // TODO: compute memory in the Arows and Acols
    
    int nx = Acols->size();
    double * workarray = (double *) calloc(Acols->size(), sizeof(double));
    double * workarray2 = (double *) calloc(Acols->size(), sizeof(double));
    
    for(int x_i=0; x_i<nx; x_i++) {
         vertex_data * var_data = g->vertex_data(x_i).as_ptr<vertex_data>();
         initialize_value(x_i, var_data, workarray, workarray2);
         
         memd += var_data->C_j_size * (sizeof(double) + sizeof(int));
         memf += var_data->C_j_size * (sizeof(float) + sizeof(int));
         memcompactd += var_data->C_j_counter * (sizeof(double) + sizeof(int));
         memcompactf += var_data->C_j_counter * (sizeof(float) + sizeof(int));
         memd += sizeof(vertex_data);
         memf += sizeof(vertex_data) - sizeof(float) * 3;
         free(var_data->C_j_indices);
         free(var_data->C_j_vals);
         if (x_i%1000 == 0) {
            std::cout << "Memory (" << x_i << ") double: " << memd << " / " << memcompactd << " float: " << memf << " / " << memcompactf << std::endl;
         }
    }   
        std::cout << "Memory TOTAL: double: " << memd << " / " << memcompactd << " float: " << memf << " / " << memcompactf << std::endl;

}

/**
 * Computes A'A[j,:]  
 */
void initialize_value(int x_i, vertex_data * var_data, double * workarray, double * workarray2) {
  int nx = Acols->size();
  
  /* Zero out workarray */
  memset(workarray2, 0, nx * sizeof(double));
  
  /* Inject values of column x_i */
  sparse_array col_i = (*Acols)[x_i];
  
  ASSERT_NE(col_i.size, -1);
  
  var_data->C_ii = 0.0;   
  var_data->lambda_cofactor = 1.0;

  timer t;
  t.start();
  
  
  
  /* Check if is first column and all values are 1.0, then this  
     is the beta_0 term. */
  if (x_i == 0 && col_i.size == ny) {
	  var_data->lambda_cofactor = 0.0;
	  printf(">> Checking if x_0 is the offset term! \n");
  	  for(int k = 0; k<col_i.size; k++) {
		 if (col_i.values[k] != 1.0) {
		    var_data->lambda_cofactor = 1.0;
        	 printf(">> Was not %d %lf\n", k, col_i.values[k]);
		    break;
		 }
      }
      if (var_data->lambda_cofactor == 0.0) {
      	printf("=== Found beta_0 ===\n");
      }	

  }
  
  
  for(int k = 0; k<col_i.size; k++) {
    int row_idx = col_i.idxs[k];
    double xv = col_i.values[k];
    
    sparse_array rowarr = (*Arows)[row_idx];
    for(int j = 0; j<rowarr.size; j++) {
       workarray2[rowarr.idxs[j]] += xv * rowarr.values[j];
    }
  }

  // when in memorySavingMode, we can use more working memory since A'*A is not stored
  var_data->C_j_size = (array_sizes[x_i] > 0 ? array_sizes[x_i] : 2048);
  var_data->C_j_counter = 0;
  var_data->C_j_indices = (int *) malloc(sizeof(int) * var_data->C_j_size);
  var_data->C_j_vals = (double *) malloc(sizeof(double) * var_data->C_j_size);
  
  var_data->C_ii = workarray2[x_i] * 2.0;

  
  for(int j = 0; j<nx; j++) {
    double AA_ij = workarray2[j];
    
    if (AA_ij != 0.0) {
      AA_ij *= 2.0;
      var_data->C_j_indices[var_data->C_j_counter] = j;
      var_data->C_j_vals[var_data->C_j_counter] = AA_ij  / var_data->C_ii;  // Divide by variance
      var_data->C_j_counter++;

      if (var_data->C_j_counter == var_data->C_j_size) 
        expand(var_data);
      if (x_i+j == 0) printf("C(1,1)=%lf\n", AA_ij);
    }
    
  }		
  if (x_i % 1000 == 0) printf("%d %d (%d) used:%ld compacted:%ld\n", procid, x_i, col_i.size, AAsize, compacted);
  array_sizes[x_i] = var_data->C_j_counter;
  
   // Predivide
   
}




/***
 * UPDATE FUNCTION
 */
 
void acquire_locks(vertex_data * var_data, int x_i) {
	safelock.lock();
	/*
    int * indices = var_data->C_j_indices;
	int last_i = -1;
    int sz = var_data->C_j_counter;

	 for(int ii=0; ii<sz; ii++) {
  	 	 int i = *(indices++);
  	 	 if (x_i>last_i && i>x_i) {
	 	 	 g->vertex_data(x_i).as_ptr<vertex_data>()->lock.lock(); // Lock myself
  	 	 }
  	 	 g->vertex_data(i).as_ptr<vertex_data>()->lock.lock();
  	 	 last_i = i;
     }*/
}
 
 
 
void release_locks(vertex_data * var_data, int x_i) {
	safelock.unlock();
	/*	
    int * indices = var_data->C_j_indices;
	int last_i = -1;
    int sz = var_data->C_j_counter;

	 for(int ii=0; ii<sz; ii++) {
  	 	 int i = *(indices++);
  	 	 if (x_i>last_i && i>x_i) {
	 	 	 g->vertex_data(x_i).as_ptr<vertex_data>()->lock.unlock(); // Unlock myself
  	 	 }
  	 	 g->vertex_data(i).as_ptr<vertex_data>()->lock.unlock();
  	 	 last_i = i;
     }*/
}
 

// Used to determine min lambda which forces all values to 0
double calc_maxlambda_xi(vertex_data * var_data, int x_i) {
	return std::abs(var_data->d_i);
}

double calc_maxlambda() {
	timer t;
	t.start();
	int nx = Acols->size();
	double maxlambda = 0;
	for(int x_i=0; x_i<nx; x_i++) {
		maxlambda = std::max(calc_maxlambda_xi(g->vertex_data(x_i).as_ptr<vertex_data>(), x_i), 
				maxlambda);
	}
	
	printf("Max lambda = %lf (took %lf secs) \n", maxlambda, t.current_time());
	return maxlambda;
}


double calc_shootDiff(vertex_data * var_data, int x_i) {
  	int sz = var_data->C_j_counter;
  	int * indices = var_data->C_j_indices;
  	double * vals = var_data->C_j_vals;
  	register double xTCj = 0; // sse_inner_d(workarray, var_data->C_j_vals, sz);
  	register double s;
  
	for(int ii=0; ii<sz; ii++) {
  	 	 s = XV[*(indices++)] * (*(vals++)); 
  	  	 xTCj += s;
   }
  
  /* ===== 2.2 Remove C_ii*x_i and d_i ===== */
   double shootDiff = xTCj - XV[x_i] - var_data->d_i / var_data->C_ii ;   //   double shootDiff = xTCj - var_data->C_ii * XV[x_i] - var_data->d_i;
   return shootDiff;
}

double soft_threshold(double _lambda, double shootDiff) {
	  double newValue;
	  if (shootDiff > _lambda) {
		newValue = (_lambda - shootDiff); 
	  } else if (shootDiff < -_lambda) {
		newValue = (-_lambda- shootDiff);
	  } else {
		newValue = 0.0;
	  }
	  return newValue;
}

double shoot(vertex_data * var_data, int x_i, double * workarray, double * workarray2) {
  if (var_data->C_j_indices == NULL) {
  	 var_data->lock.lock();
      if (var_data->C_j_indices == NULL) {
	    initialize_value(x_i, var_data, workarray, workarray2);
	  }
   	 var_data->lock.unlock();
 
  }
  
  if (var_data->C_j_counter == 0) {
  	  return 0.0;
  }
  
  bool acquired_locks = false;
  if (enable_locking) {
  	 acquire_locks(var_data, x_i);
  	 acquired_locks = true;
  }
  
  // Note: divide by C_ii for numerical stability!!
  double _lambda = lambda/var_data->C_ii * var_data->lambda_cofactor;
  
  /* ===== 2.1 Compute dotproduct of x's and column j of C ==== */
  double shootDiff = calc_shootDiff(var_data, x_i);
  
  if ( std::isinf(shootDiff) || std::isnan(shootDiff)) {
  	  assert(false);
  	  return 1e10;
  }
  
  

  if (std::abs(shootDiff) > 1e6) {
  	printf("Warning ::: very big value %d => %lf\n", x_i, shootDiff); 
  }
  
  // Dense vector acquire locks
  double newValue = soft_threshold(_lambda, shootDiff);
    
  double possibleChange = -1;
  if (!acquired_locks && mode > 1) {
     possibleChange =  std::abs(newValue- XV[x_i]);
     if ((XV[x_i] == 0.0  && std::abs(shootDiff) > _lambda) ||(XV[x_i] != 0 && std::abs(shootDiff) <= _lambda)
     		|| (XV[x_i] != 0.0 && (possibleChange > 1.0 || possibleChange/std::abs(XV[x_i]) > 0.1))) {
  	 	
  	 	 double first_shootDiff = shootDiff;
  	 	 acquire_locks(var_data, x_i);
  	 	 acquired_locks = true;
	 	 double prev = shootDiff;
	 	 do {
	 	 	 prev = shootDiff;
		 	 shootDiff = calc_shootDiff(var_data, x_i);
		 	// printf("%d %lf - %lf\n", x_i, shootDiff, prev);
  		 } while (std::abs(prev-shootDiff) > 1e-6);
  		 
  		 // printf("%d - %d Significant change, acquire lock and compute again %lf  %lf => %lf (lamb:%lf) %lf %lf\n",
  	    // 		x_i, thread::thread_id(), XV[x_i],  first_shootDiff, shootDiff, _lambda, possibleChange/std::abs(XV[x_i]), possibleChange);
  	     		
  	      newValue = soft_threshold(_lambda, shootDiff);
  	 }
  }
  
  /* ===== 2.3 Subdifferential check and update ====== */
  //double delta =   //std::abs(newValue/XV[x_i]  - 1.0);
  //if (convergence_mode == COMPRESSED_SENSING) 
  double  delta = std::abs(newValue - XV[x_i]);
  if (XV[x_i] == 0.0) delta = std::abs(newValue-XV[x_i]);
  
  if (std::isnan(newValue) || std::isinf(newValue)) {
  	 std::cout << "NaN/inf Value encountered!" << std::endl;
  	 assert(false);
  }

  
  if (delta != 0.0) {
     XV[x_i] = newValue;
  }
  
  if (memorySavingMode) {
  	   free(var_data->C_j_indices);
  	   free(var_data->C_j_vals);
  	   var_data->C_j_indices = NULL;
  	   var_data->C_j_vals = NULL;
  	   var_data->C_j_size = 0;
  }
  
  if (acquired_locks) {
  	 release_locks(var_data,x_i);
  }
  
  if (std::isnan(delta)) delta = 1e10;
  
  return delta;
}

// Run shooting until convergence
int endgame() {
    ASSERT_EQ(numprocs, 1);
    double change=0;
    int nx = Acols->size();
    mode = 1;
    printf("Start sequential phase (nx=%d)\n", nx);

    double * workarray = (double *) calloc(Acols->size(), sizeof(double));
    double * workarray2 = (double *) calloc(Acols->size(), sizeof(double));
    // Set set_to_zero to 0
	 memset(set_to_zero, 0, sizeof(int) * Acols->size());
		      
    bool notconverged;

    int nzeros, seqiters = 0;
    double maxchange = 0.0, lastobjective = -1;
    int last_n_nonzeros = -1;
    do {
   	    maxchange = 0.0;
        seqiters++;
        change = 0;
       // timer tt; tt.start();
        nzeros = 0;
       
        for(int i=0; i<nx; i++) {
        	if (set_to_zero[i] < MAXZEROSWEEPS) {
	        	double c = shoot(g->vertex_data(i).as_ptr<vertex_data>(), i, workarray, workarray2);
    	        maxchange = std::max(c, maxchange);
        	    change += c;
            	if (c > 0 && XV[i] == 0) {
            		nzeros++;
            		set_to_zero[i]++;
            	}
            	endgame_iterations++;
            }
        }
        
        if (endgame_iterations >= max_endgame_iterations) {
            printf("Endgame iteration limit reached, ready.\n");
            break;
        }
        /* ===== CONVERGENCE CHECK ====== */
       // printf("%d %lf %lf  %lf s maxchange: %lf\n", seqiters, change, lambda, tt.current_time(), maxchange);
      notconverged = false;
      switch(convergence_mode) {
			case ACCURATE:
				notconverged = (nzeros != 0);
				notconverged = notconverged || (maxchange > (k == 0 ? DESIRED_ACCURACY  : (DESIRED_ACCURACY + (DESIRED_ACCURACY*10-DESIRED_ACCURACY)/K*k) ));
		  //	  notconverged = notconverged || (k == 0 && maxchange > 1e-4);
				 break;
			case COMPRESSED_SENSING:
				notconverged = notconverged || (maxchange > (k == 0 ? DESIRED_ACCURACY  : (DESIRED_ACCURACY + (DESIRED_ACCURACY*10-DESIRED_ACCURACY)/K*k) ));
				break;
	  }	
      // Faster convergence on high k
	  if ( notconverged && k>0 && seqiters == std::min(100, (K-k)*2)) {
	  	 printf("Enough rounds this lambda: %d\n", seqiters);
	  	 notconverged = false;
	  }  
 	    
    } while (notconverged);
    printf("Number of sequential iterations: %d\n", seqiters);
    
    free(workarray);
    free(workarray2);
    return seqiters;
}

double lastobjective = 0;
int last_n_nonzeros = -1;
bool CONVERGED = false;
bool * converged_list;

double * last_x = NULL;
double * saved_x = NULL;

void aapolasso_update_function(gl_types::iscope &scope,
                               gl_types::icallback &scheduler,
                               gl_types::ishared_data* shared_data)  {
  
  double totalChangeInXs = 0.0;
  
  /* ===== 1. PREPARATION ====== */  
  proc_vertex_data * mydata = scope.vertex_data().as_ptr<proc_vertex_data>();
  
  
  if (mydata->myvars == NULL) {
    mydata->myvars = new std::vector<int>();
    /* Collect x_i's that belong to me. (This is bit redundant, but makes 
     code less nested) */    	
    foreach(edge_id_t eid, scope.out_edge_ids()) {
      proc_to_x_edge_data * ed = scope.edge_data(eid).as_ptr<proc_to_x_edge_data>();
      if (ed->my_vertex) {
        mydata->myvars->push_back(scope.target(eid));
      }
    }
    printf("myvars->size(), %ld\n",mydata-> myvars->size());
    mydata->workarray = (double *) calloc(Acols->size(), sizeof(double));
    mydata->workarray2 = (double *) calloc(Acols->size(), sizeof(double));
  }
  
  printf("%d start\n", mydata->process_id);
  
  std::vector<int> myvars = *(mydata->myvars);

  int rounds_this_lambda = 0;
	
  while(!CONVERGED) {
	   int optimized = 0;
	  //int conflicts = 0; 
	  
	  // Shuffle
	  std::random_shuffle(myvars.begin(), myvars.end());
	  
	  double maxchange = 0.0;
	  
	  int newzeros = 0;
	  
	  /* ===== 2. SHOOTING ====== */
	  double _startlambda = lambda;
	  foreach(int x_i, myvars) {
		  vertex_data * var_data =  scope.neighbor_vertex_data(x_i).as_ptr<vertex_data>();
	
		  if (set_to_zero[x_i] < MAXZEROSWEEPS) {
			double c = shoot(var_data, x_i, mydata->workarray, mydata->workarray2);
			totalChangeInXs += c;
			if (c > 0 && XV[x_i] == 0.0) {
				newzeros++;
			    set_to_zero[x_i]++;
			} 
				//else if (c>0) set_to_zero[x_i] = 0;
			maxchange = std::max(c, maxchange);
			optimized++;
			mydata->num_of_iters++;
		  }
		  if (_startlambda != lambda) {
			maxchange = 1e10;
			break;
		  }
		  
	  }
	
	  //printf("Conflicts: %d / %d\n", conflicts, mydata->myvars->size());
	  
	  //printf("Total change :%lf / lambda %lf\n", totalChangeInXs, lambda);
	  
	  int startk = k;
	
	  if (converged_list[mydata->process_id] == false) {
	  	 rounds_this_lambda++;
	  }
	
	  /* ==== 3. CONVERGENCE CHECKS ==== */
	  bool notconverged = false;
	  
	  switch(convergence_mode) {
		case ACCURATE:
		  notconverged = (newzeros != 0);
				notconverged = notconverged || (maxchange > (k == 0 ? DESIRED_ACCURACY  : (DESIRED_ACCURACY + (DESIRED_ACCURACY*10-DESIRED_ACCURACY)/K*k) ));
		 		 break;
			case COMPRESSED_SENSING:
				notconverged = notconverged || (maxchange > (k == 0 ? DESIRED_ACCURACY  : (DESIRED_ACCURACY + (DESIRED_ACCURACY*10-DESIRED_ACCURACY)/K*k) ));
				break;
	} 
	  
	  // Faster convergence on high k
	  if ( notconverged && k>0 && rounds_this_lambda == std::min(100, (K-k)*2)) {
	  	 printf("Enough rounds this lambda: %d\n", rounds_this_lambda);
	  	 notconverged = false;
	  }
	  
	  
	  // Monitor for divergence
	  if (notconverged && k == 0  &&  mydata->process_id==0 && rounds_this_lambda > 500 && (rounds_this_lambda % 100 == 0)) {
	  		 int nx = Acols->size();
	  		 if (last_x == NULL) {
	  		 	last_x = (double *) malloc(sizeof(double) * nx);
	  		 	saved_x =(double *) malloc(sizeof(double) * nx);
	  		 }
			  memcpy(last_x, XV, sizeof(double) * nx);

		
	  		 int n_nonzeros;
    		 double least_sqr = 0, _lambda, penalty;
      		 double obj = compute_objective(last_x, nx, &n_nonzeros, &_lambda, &least_sqr, &penalty);  
      		 if (lastobjective < obj*0.999 && lastobjective > 0) {
      		 	printf("%d Looks we are diverging. Enable locking %lf %lf  maxchange:%lf\n", rounds_this_lambda, lastobjective, obj, maxchange);
      		  	enable_locking = true;
	    	    // Restore best known vector	
	    	    safelock.lock();
	    	    memcpy(XV, saved_x, sizeof(double) * nx);
				safelock.unlock();
      		 } else {
      		 	enable_locking = false;
		 	    lastobjective = obj;
    			memcpy(saved_x, last_x, sizeof(double) * nx);

      		 }
	  }
	  
	 if (!notconverged && !converged_list[mydata->process_id]) {
	 // We are converged. If everyone else is converged as well, decrease lambda  
		  int nconverged = convcounter.inc();
		  converged_list[mydata->process_id] = true;
		  printf("Converged: %d/%d ) (%d) %lf\n", nconverged, num_of_processing_nodes, mydata->process_id, maxchange);
		  rounds_this_lambda = 0;
		  if (nconverged == num_of_processing_nodes) {
			 
			 // Decrease lambda if all converged!
			 if (k>0) {
			   if(procid == 0) {
			          double _prevlambda = lambda;
					  lambda = lambda_min * pow(alpha, --k);  // Should use shared data
					  sendlock.lock();
					  set_lambda(lambda, k);
					  std::cout << ">>> ======= LAMBDA = " << lambda  << " k = " << k<< std::endl;
					  	  
					  /* Update max value */
					  MAX_VALUE = MAX_VALUE / _prevlambda * lambda;
					  std::cout << ">>> New max: " << MAX_VALUE << std::endl;
					  
					  schedule_remotes();
					  sendlock.unlock();
				  }
				  // Set set_to_zero to 0
				 	 memset(set_to_zero, 0, sizeof(int)*Acols->size());
				  
				  } else {					
					CONVERGED = true;
			  }
			  convcounter = 0; 
			  memset(converged_list, 0, num_of_processing_nodes); 
		   }
		
	  }
	 }
	 
	  printf("End: %d\n", mydata->process_id);
}

/**
 * aapolasso APPLICATION CLASS
 */
class aapolassoapp : public iapp<blob_graph> {
  
private:

  std::string inputfile;
  bool disable_objcalc;
  
public: 
  
  
  /** CONSTRUCTOR **/
  aapolassoapp(std::string inputfile, int _num_of_processing_nodes, graphlab::distributed_control * _dc, 
  		std::string mode_str, std::string disable_objcalcs, double obj_termination_threshold, int _K) {
    this->inputfile = inputfile;
    num_of_processing_nodes = _num_of_processing_nodes;
    dc = _dc;
    
    K = _K;
    DESIRED_ACCURACY = obj_termination_threshold;
    
    if (mode_str == "cs") {
    	convergence_mode = COMPRESSED_SENSING;
    }
    else if (mode_str == "accurate") {
    	convergence_mode = ACCURATE;
    } else {
    	assert(false);
    }
    
    if (disable_objcalcs == "false") { disable_objcalc = false; }
    else disable_objcalc = true;
    
    printf("Mode: %d = %s\n", convergence_mode, mode_str.c_str());
  }
  
  ~aapolassoapp() {
    delete engine;
    delete g;
  }
  
  // Removes \n from the end of line
	void FIXLINE(char * s) {
		int len = strlen(s)-1; 	  
		if(s[len] == '\n') s[len] = 0;
	}
	
  
	void create_processing_nodes() {
    for(unsigned int i=0; i<num_of_processing_nodes; i++) {
      proc_vertex_data pvd;
      pvd.process_id = i;
      pvd.myvars = NULL;
      pvd.workarray = NULL;
      pvd.num_of_iters = 0;
      g->add_vertex(blob(sizeof(proc_vertex_data), &pvd));
    }
    converged_list = (bool*) calloc(num_of_processing_nodes, sizeof(bool));
  }   
  
  void compact_sparse_array(sparse_array * sarr) {
  			int * newarr_idx = (int *) malloc(sizeof(int) * sarr->size);
			memcpy(newarr_idx, sarr->idxs, sizeof(int)*sarr->size);
			double * newarr_val = (double *) malloc(sizeof(double) * sarr->size);
			memcpy(newarr_val, sarr->values, sizeof(double)*sarr->size);
			
			free(sarr->idxs);
			sarr->idxs = newarr_idx;
			free(sarr->values);
			sarr->values = newarr_val;
		
			compacted += (sizeof(double)+sizeof(int)) * (sarr->capacity-sarr->size);
			Asize += sarr->size;
			sarr->capacity = sarr->size;
  }
  
  void add_to_sparse_array(sparse_array * sarr, int idx, double val) {
    sarr->values[sarr->size] = val;
		sarr->idxs[sarr->size++] = idx;
    
		if (sarr->size == sarr->capacity) {
			int * newarr_idx = (int *) malloc(sizeof(int) * sarr->capacity * 2);
			memcpy(newarr_idx, sarr->idxs, sizeof(int)*sarr->capacity);
			double * newarr_val = (double *) malloc(sizeof(double) * sarr->capacity * 2);
			memcpy(newarr_val, sarr->values, sizeof(double)*sarr->capacity);
			
			sarr->capacity *= 2;
			free(sarr->idxs);
			sarr->idxs = newarr_idx;
			free(sarr->values);
			sarr->values = newarr_val;
		}
  }
  
  /**
   * Load graph
   */
  void load() {
    std::cout << "Opening " << inputfile << std::endl;
    FILE * f = fopen(inputfile.c_str(), "r");
    if (f == NULL) {
      std::cout << "File not found! " << std::endl;
      exit(0);
    }
    
    numprocs = dc->numprocs();
    printf("Num procs: %d \n", numprocs);
    unsigned int first_x;
    unsigned int last_x;
    
    
    g = new blob_graph();
    
    char s[255];
	char delims[] = ",";	
    char *t = NULL;
    std::string mode;
	unsigned int n = 0; 
	unsigned int vid = 0;
	unsigned int nx;

    while(fgets(s, 255, f) != NULL) {
			FIXLINE(s);
			
			if (n == 0) {
				t = strtok(s, delims);
				mode = std::string(t);
				if (mode == "y" || mode == "d" || mode == "m") {
					t = strtok(NULL, delims);
					sscanf(t, "%d", &n);
					std::cout << mode << " length=" << n << std::endl;
					if (mode == "y") {
						ny = n;
						y.reserve(ny);
					}
					if (mode == "d") {
						/* Create x vertices */
						nx = n;
						XV = (double *) malloc(n * sizeof(double));
						running = (bool *) calloc(n, sizeof(bool));
						todo = (bool *) calloc(n, sizeof(bool));
            			set_to_zero =  (int *) calloc(n, sizeof(int));
						for(unsigned int j=0; j<n; j++) {
							vertex_data vd;
							vd.value = 0.0; // TODO: INITIAL VALUE RANDOM??
							XV[j] = 0.0; // 1.0;
							vd.C_j_size = 0;
                            vd.C_j_counter = 0;
                            vd.C_j_indices = NULL; //(idx_value *) malloc(vd.C_j_size * sizeof(idx_value));
              				vd.C_j_vals = NULL;
							g->add_vertex(blob(sizeof(vertex_data), &vd));
						}
						std::cout << "Created vertices" << std::endl;
            
                        first_x =  FIRST_X_OF_PROCESS(procid, numprocs);
                        last_x = LAST_X_OF_PROCESS(procid, numprocs);
                         printf("====== Process %u/: %u => %u \n", procid, first_x, last_x);
            
						create_processing_nodes();
            
						/* Create outbound edges: from processing nodes to X_i */
            
						for(unsigned int j=first_x; j<=last_x; j++) {
							//printf("%d => %d\n", j, PROC_NODE_FOR_X(j));
							for(unsigned int pv=0; pv<num_of_processing_nodes; pv++) {
								proc_to_x_edge_data ped;
								ped.my_vertex = (PROC_NODE_FOR_X(j) == pv+n); 
								g->add_edge(pv+n, j, blob(sizeof(proc_to_x_edge_data), &ped));
							}
						}
					}
					vid = 0;
				} else if (mode == "A") {
					t = strtok(NULL, delims);
					sscanf(t, "%d", &n);
					t = strtok(NULL, delims);
					sscanf(t, "%d", &nx);
					t = strtok(NULL, delims);
					sscanf(t, "%d", &ny);
					std::cout << mode << " length=" << n << std::endl;
					std::cout << "Num of x: " << nx << std::endl;
					
					/* Initialize rows and cols */
					std::cout << "Initialize cols/rows" << std::endl;
					Acols = new std::vector<sparse_array>();
					Acols->reserve(nx);
					Arows = new std::vector<sparse_array>();
					Arows->reserve(ny);
					for(unsigned int i=0; i<nx; i++) {
						// TODO: Clean up
						sparse_array sarr;
						sarr.size = 0;
						sarr.capacity = 10;
						sarr.idxs = (int *) malloc(sizeof(int)*sarr.capacity);
						sarr.values = (double *) malloc(sizeof(double)*sarr.capacity);
						Acols->push_back(sarr);
					}
					for(unsigned int i=0; i<ny; i++) {
						sparse_array sarr2;
						sarr2.size = 0;
						sarr2.capacity = 10;
						sarr2.idxs = (int *) malloc(sizeof(int)*sarr2.capacity);
						sarr2.values = (double *) malloc(sizeof(double)*sarr2.capacity);
						Arows->push_back(sarr2);
					}
					
					array_sizes = (int *) calloc(nx, sizeof(int));
				}
			} else {
				n--;
				if (mode == "y") {
					// do nothing
					double f;
					sscanf(s, "%lf", &f);
					y.push_back(f);
				} else if (mode == "d") {
					double f;
					sscanf(s, "%lf", &f);
					g->vertex_data(vid++).as_ptr<vertex_data>()->d_i = f;
				} else if (mode == "A") {
					t = strtok(s, delims);
					uint64_t idx;
					double value;
					sscanf(t, "%ld", &idx);
					idx--;  // Matlab uses 1-index
					t = strtok(NULL, delims);
					sscanf(t, "%lf", &value);
					
					
					// Create edge from x_i to processing node
					int row_id = idx%ny;
					int col_id = idx/ny;
					
					sparse_array * sarr = &((*Acols)[col_id]);
					add_to_sparse_array(sarr, row_id, value); 
					
					sarr = &((*Arows)[row_id]);
					add_to_sparse_array(sarr, col_id, value);
				}
			}
		}
		
		
		/*
		// Compact 
		for(unsigned int i = 0; i<Acols->size(); i++) {
			compact_sparse_array(&(*Acols)[i]);
		}
		for(unsigned int i = 0; i<Arows->size(); i++) {
			compact_sparse_array(&(*Arows)[i]);
		}
		std::cout << "Matrix A size: " << Asize << "  compacted " << compacted << std::endl;

		*/
		
		
		if (nx > 1e6) {
			memorySavingMode = true;
			std::cout << "************* MEMORY SAVING MODE ****************" << std::endl;
		}
		
		std::cout << "Num edges: " << g->num_edges() << std::endl;
		std::cout << "Num vertices: " << g->num_vertices() << std::endl;
		std::cout << "Acols:" << Acols->size() << std::endl;
  }
  
  
  /** 
   * ==== SETUP AND START
   */
  void start() {
    
    procid = dc->procid();
    
    if (num_of_processing_nodes == 0) {
      num_of_processing_nodes = opts.ncpus;
    }
    
    std::cout << "Num of distributed processors: " << procid << " / " << dc->numprocs() << std::endl;
    
    /**** CREATE GRAPH ****/
    load();
    
    char pfilename[255];
    sprintf(pfilename, "%s_%d_progress_%ld.txt", inputfile.c_str(), opts.ncpus, time(NULL));
    progressFile = fopen(pfilename, "w");
    mode = opts.ncpus;

    
    // Compute how much memory the A'A matrix would take!
    if (opts.extra == "computemem") {
        computemem();
        return;
    }    
    // printf("C[1,1] = %lf, C[127,128]=%lf\n", AtAij(0,0), AtAij(126,127)); 
    
    // Testing
    /**** GRAPHLAB INITIALIZATION *****/
    engine = create_engine(g, NULL);
    
    if ((opts.extra.length() > 0)) {
		sscanf(opts.extra.c_str(), "%lf", &lambda_min);
    }
    
    std::cout << "====== Lambda = " << lambda_min << std::endl;
    
    
    int n = g->num_vertices();
      /// Barrier to wait everyone has loaded their part
    dc->barrier();
    
    
    /// Timing
    timer t;
    t.start();
    
    /**** RUN GRAPHLAB IN LOOP OVER THE REGULARIZATION PATH **/
    
    // TODO: make configurable
    if (K == 0) {
    	K = 1+(int)(Acols->size() / 2000);  // QUESTIONNABLE!!
    }
    printf("K = %d\n", K);
    printf("obj termination threshold: %lf\n", DESIRED_ACCURACY);
    
    lambda_max = calc_maxlambda();
    
    alpha = pow(lambda_max/lambda_min, 1.0/(1.0*K));
    k = K-1; 
    lambda = lambda_min * pow(alpha, k);
 
    /**** CREATE INITIAL TASKS ******/
    first_proc_node = n-num_of_processing_nodes;
    last_proc_node = n;
    for(int vid=n-num_of_processing_nodes; vid<n; vid++) {
      engine->get_scheduler().add_task(gl_types::update_task(vid, aapolasso_update_function), 1.0);
    }
    
    /* Compute initial objective */ 
    {
      int n_nonzeros;
      double least_sqr = 0, _lambda, penalty;
      MAX_VALUE = compute_objective(XV, Acols->size(), &n_nonzeros, &_lambda, &least_sqr, &penalty);   
    }
    printf("MAX VALUE = %lf\n", MAX_VALUE);
    
    /**** RUN GRAPHLAB ****/
 	
 	if (opts.extra == "safe" && opts.ncpus > 1) {
 		enable_locking = true;
 	}
 	
 	if (numprocs > 1) {
 	    value_broadcaster * broadcaster = new value_broadcaster();
 	    thread * thrbr = new  thread(broadcaster);
 	    thrbr->start();
 	}
 	
 	/* Start convergence counter */
 	thread * thr = new  thread(new objective_computation());
 	if (disable_objcalc == false) {
 		thr->start();
 	}
 	if (opts.extra == "endgame" || (numprocs == 1 && opts.ncpus == 1)) {
 	    // Need to write for benchmark suite
 	    max_endgame_iterations = 1e9;
 	    int seqiters = 0;
 		for(k=K; k>=0; k--) {
 		 	
 			lambda = lambda_min * pow(alpha, k);
 			printf("Lambda %lf\n", lambda);
	 		seqiters += endgame();
 		}
 		
 		
 		FILE * F = fopen(".runstats.R", "w");
        fprintf(F, "engine=\"sequential\"\nncpus=0\ntaskcount=%d\nwork=%d\nresidual=0\nmemory_writes_mb=0\nmemory_reads_mb=0\n", endgame_iterations, endgame_iterations);
        fclose(F);
 	    
 	} else {
 	    while(k>0) {
    	    // Engine may stop prematurarely, but then it will be
        	// restarted by message handlers
        	engine->start();
    	}
    	
    	printf("========== WAITING AT BARRIER ========\n");
    	
    	/**** POST-PROCESSING *****/
   		 if (numprocs > 1) {
	    	dc->barrier();
    	}	
    	
    	
    	 if (convergence_mode == ACCURATE && numprocs == 1) {
    	    std::cout << "******* Endgame ********" << std::endl;
	        endgame();
        }
    }
    

   
    // If running in shared memory, run shooting sequentially
    // in the end
    
//    if (numprocs == 1)
		terminated = true;	
        
    double runtime = t.current_time();
    set_status("Node %d Finished in %lf", procid, runtime);
    write_benchmark_value("execution_time", runtime*1.0);

    if (procid == 0) {  	
      char outfile[255];
      sprintf(outfile, "%s_out_d%d_%ld.m", inputfile.c_str(), numprocs, opts.ncpus);
      
      FILE * outf = fopen(outfile, "w");
      
      fprintf(outf, "mrecon=[");
      
      /* Get graph data and do some output */
      for(unsigned int vid=0; vid<n-num_of_processing_nodes; vid++) {
        //vertex_data * vdata = (vertex_data*) g->vertex_blob(vid).data;
        fprintf(outf, " %lf", XV[vid]);
        if (vid % 30 == 0) {
          fprintf(outf, "...\n");
        }
      }
      
      fprintf(outf, "];");
      fclose(outf);
      
      /** Count iterations **/
     
      int iters = 0;
      for(int vid=n-num_of_processing_nodes; vid<n; vid++) {
        proc_vertex_data * pvd = g->vertex_data(vid).as_ptr<proc_vertex_data>();
        iters += pvd->num_of_iters;
        printf("%d : %d\n",vid, pvd->num_of_iters);
      }
         if (iters>0) { // avoid writing zero when running only endgame. this is becoming really dirty!
          write_benchmark_value("taskcount", endgame_iterations + iters*1.0);
      }
     }
      delete(thr);
	  
	  int n_nonzeros;
      double least_sqr = 0, _lambda, penalty;
      double obj = compute_objective(XV, Acols->size(), &n_nonzeros, &_lambda, &least_sqr, &penalty);      
      write_benchmark_value("residual", obj);
 	  
	  fclose(progressFile);
  }
  
}; // end class


#endif

