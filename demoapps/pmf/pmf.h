#ifndef PMF_H__	 

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>

/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU

See algrithm description and explanation in: Liang Xiong, Xi Chen, Tzu-Kuo Huang, Jeff Schneider, Jaime G. Carbonell, Temporal Collaborative Filtering with Bayesian Probabilistic Tensor Factorization. In Proceedings of SIAM Data Mining, 2010.

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
  vec pvec; //vector of learned values U,V,K
  double rmse; //root of mean square error
  int num_edges; //number of edges
  vertex_data(){
    pvec = zeros(D);
    rmse = 0;
    num_edges = 0;
  }


  void save(graphlab::oarchive& archive) const {  
    serialize(archive, pvec._data(), sizeof(double)*D); 
    archive << rmse << num_edges; 
  }  
   
  void load(graphlab::iarchive& archive) {  
    deserialize(archive, pvec._data(), sizeof(double)*D); 
    archive >> rmse >> num_edges;  
  }
};

struct edge_data {
  double weight;  //observation 
  double time; //time of observation (for tensor algorithms)
  double avgprd;
  edge_data(){ weight = 0; time = 0; avgprd = 0;}

 void save(graphlab::oarchive& archive) const {  
    archive << weight << time << avgprd;; 
  }  
   
  void load(graphlab::iarchive& archive) {  
    archive >> weight >> time >> avgprd;  
  }
};

//containiner for handling multiple edge
struct multiple_edges{
  std::vector<edge_data> medges;
 void save(graphlab::oarchive& archive) const {  
       archive << medges; 
  }  
   
  void load(graphlab::iarchive& archive) {  
    for (int i=0; i< medges.size(); i++)
    	archive >> medges;  
  }
};
 
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

//faster reset for vectors
void _zeros(vec & pvec, int d){
  assert(pvec.size() == d);
  memset(pvec._data(), 0, sizeof(double)*d);
}
 
void _zeros(mat & pmat, int rows, int cols){
  assert(pmat.size() == rows*cols);
  memset(pmat._data(), 0, sizeof(double)*rows*cols);
} 


enum testtype{
    TRAINING = 0,
    VALIDATION = 1,
    TEST = 2
};

const char * testtypename[] = {"TRAINING", "VALIDATION", "TEST"};
//run modes

enum runmodes{
   ALS_MATRIX = 0,
   BPTF_MATRIX = 1,
   BPTF_TENSOR = 2,
   BPTF_TENSOR_MULT = 3,
   ALS_TENSOR_MULT = 4
};
/*#define ALS_MATRIX 0  //alternating least squares for matrix factorization
#define BPTF_MATRIX 1 //baysian matrix factorization
#define BPTF_TENSOR 2 //bayesian tensor factorization
#define BPTF_TENSOR_MULT 3 //bayesian tensor factorization, with support for multiple edges between pair of ndoes
#define ALS_TENSOR_MULT 4
*/
const char * runmodesname[] = {"ALS_MATRIX", "BPTF_MATRIX", "BPTF_TENSOR", "BPTF_TENSOR_MULT", "ALS_TENSOR_MULT"};
#endif

