#ifndef PMF_H__	 

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>

/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU


*/
#include <vector>


int MAX_ITER=10; //maximal number of iterations to run
int BURN_IN =10; //burn-in priod (for MCMC sampling - optional)
int D=20;         //number of features
bool FLOAT=false; //is data in float format
double LAMBDA=1;//regularization weight

typedef double  sdouble; 
using namespace itpp;

/** Vertex and edge data types **/
struct vertex_data {
  vec pvec;
  double rmse;
  int num_edges;
  vertex_data(){
    pvec = zeros(D);
    rmse = 0;
    num_edges = 0;
  }
}; 

struct edge_data {
  double weight; 
  double time;
  double avgprd;
  edge_data(){ weight = 0; time = 0; avgprd = 0;}
};

//containiner for handling multiple edge
struct multiple_edges{
  std::vector<edge_data> medges;
};
 

/*
sdouble* ones(int n, sdouble val) {
	assert(n>0);
	sdouble * ret = new sdouble[n];
	for (int i=0; i< n; i++)
		ret[i] = val;
	return ret;
}
void ones(double * x, int len){
	for (int i=0; i<len; i++)
		x[i] = 1.0;
}
void ones(double * x, int len, double val){
	for (int i=0; i<len; i++)
		x[i] = val;
}*/


inline double rmse(vec& x1, vec& x2, vec& x3, int len, double val){

	double sum = 0;
	for (int i=0; i< len; i++)
	sum += (x1[i] * x2[i] * x3[i]);

	sum = pow(sum - val, 2);
	return sum;
}
      
inline double rmse(vec& x1, vec& x2, int len, double val, double & sum){
	
	sum = 0;
	double ret = 0;
	for (int i=0; i< len; i++){
	ret = (x1[i] * x2[i]);
	sum+= ret;
	}
	
	return pow(sum - val, 2);
}

inline double rmse(vec& x1, vec& x2, vec *x3, int len, double val, double &sum){
	if (x3 == NULL) //matrix	
		return rmse(x1,x2,len,val,sum);

	//assert(val>=0 && val <= 5);
	sum = 0;
	double ret = 0;
	for (int i=0; i< len; i++){
	ret = (x1[i] * x2[i] * x3->get(i));
	//assert(ret>0);  //TODO
	sum+= ret;
	}
	
	return pow(sum - val, 2);
}

inline sdouble square_sum(vec& x, int len){
	sdouble xs = 0.0;
	for(int i=0; i<len; i++) {
		xs += x[i] * x[i];
	}
	return  xs;
}

 double get_rmse(const vertex_data & v){
    return v.rmse;
 }



//run modes
#define ALS_MATRIX 0  //alternating least squares for matrix factorization
#define BPTF_MATRIX 1 //baysian matrix factorization
#define BPTF_TENSOR 2 //bayesian tensor factorization
#define BPTF_TENSOR_MULT 3 //bayesian tensor factorization, with support for multiple edges between pair of ndoes
#define ALS_TENSOR_MULT 4


#endif

