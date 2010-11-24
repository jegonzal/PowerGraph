// written by Danny Bickson, CMU #include <cstdio>
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>

#include <fstream>
#include <cmath>
#include <cstdio>
#include <cfloat>

#include <graphlab.hpp>
#include <graphlab/app_support/engine_options.hpp>

#include <graphlab/app_support/engine_options.hpp>
extern "C" {
 #include "../../apps/linear_algebra/random.hpp"
};

typedef double  sdouble; 

#include "../../apps/linear_algebra/vecops.hpp"
#include "../../apps/linear_algebra/vecutils.hpp"
#include "../../apps/linear_algebra/itppvecutils.hpp"
#include "../../apps/linear_algebra/prob.hpp"

//#include <itpp/stat/misc_stat.h>
#include <graphlab/schedulers/round_robin_scheduler.hpp>
#include <graphlab/macros_def.hpp>

//run modes
#define ALS_MATRIX 0  //alternating least squares for matrix factorization
#define BPTF_MATRIX 1 //baysian matrix factorization
#define BPTF_TENSOR 2 //bayesian tensor factorization
#define BPTF_TENSOR_MULT 3 //bayesian tensor factorization, with support for multiple edges between pair of ndoes
#define ALS_TENSOR_MULT 4


#define SAMPLE_U_NODE_OFFSET M+N+K
#define SAMPLE_V_NODE_OFFSET M+N+K+1
#define SAMPLE_T_NODE_OFFSET M+N+K+2
#define T_NODE_OFFSET M+N
#define U_NODE_OFFSET 0
#define V_NODE_OFFSET M

bool BPTF = true;
#define D 6 //diemnsion for U,V
int options;
timer gt;
using namespace itpp;
using namespace std;
bool debug = false;

std::vector<edge_id_t> * edges;
std::string infile;

extern bool finish; //defined in convergence.hpp
int iiter = 1;//count number of time zero node run

/* Variables for PMF */
int M,N,K,L;//training size
int Le = 0; //test size
double pU = 10; //regularization for matrix
double pT = 10;
double pV = 10;
double muT = 1;
double * _ones;
vec vones; 
mat eDT; 
mat dp;


/* variables for BPTF */
double nuAlpha = 1;
double Walpha = 1;
double mu0 = 0;
double mu0T = 1;
double nu0 = D;
double alpha = 0;
double beta = 1;
double beta0 = 1; //TODO
mat W0;
mat W0T;
double iWalpha;
mat iW0;
mat iW0T;
mat A_U, A_V, A_T;
vec mu_U, mu_V, mu_T;

bool record_history = false;
int BURN_IN =10;
bool tensor = true;
double counter[20];


/** Vertex and edge data types **/
struct vertex_data {
    double pvec[D];
    double rmse;
    int num_edges;
    vertex_data(){
       memset(pvec, 0, sizeof(double)*D);
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

struct multiple_edges{
   std::vector<edge_data> medges;
   bool control;
   multiple_edges(){ control = false;}
             
};
 
//vertex_data * times = NULL;

typedef graphlab::graph<vertex_data, multiple_edges> graph_type;
typedef graphlab::types<graph_type> gl_types;
gl_types::iengine * engine;
graph_type g;
graph_type g1;
gl_types::thread_shared_data sdm;


double get_rmse(const vertex_data & v){
     return v.rmse;
}

const size_t RMSE = 0;

 

static void apply_func(size_t index,
                                  const gl_types::ishared_data& shared_data,
                                  graphlab::any& current_data,
                                  const graphlab::any& new_data) {

  double rmse = new_data.as<double>();
  rmse = sqrt(rmse/L);
  std::cout << "Training RMSE is : " << rmse << std::endl;
  // write the final result into the shared data table
  current_data = (double)rmse;
}



void init_self_pot(){
  //assert(BPTF);

  W0 = eye(D);
  W0T = eye(D);
  iWalpha = 1.0/Walpha;
  iW0 = inv(W0);
  iW0T = inv(W0T);
  nu0 = D;

  A_U = eye(D);
  A_V = eye(D);
  A_T = eye(D);

  mu_U = zeros(D); mu_V = zeros(D); mu_T = zeros(D);
  printf("nuAlpha=%g, Walpha=%g, mu=%g, muT=%g, nu=%g, "
        "beta=%g, W=%g, WT=%g BURN_IN=%d\n", nuAlpha, Walpha, mu0, 
         mu0T, nu0, beta0, W0(1,1), W0T(1,1), BURN_IN);


  //test_randn(); 
  //test_wishrnd();
  //test_wishrnd2(); 
  //test_chi2rnd();
  //test_wishrnd3();
  //test_mvnrndex();
}

 
/**
 * pmf APPLICATION CLASS
 */
    
    /** CONSTRUCTOR **/
  void init_pmf() {
    eDT = itpp::eye(D)*pT;
    _ones = gl::ones(D);
    vones = itpp::ones(D);
  }
    
  void load_pmf_graph(const char* filename, graph_type * g, bool flag);    
  void calc_T(int id);    
  double calc_obj(/*double res*/);
  void last_iter();


void sample_alpha(double res2){
 double res = powf(sdm.get(RMSE).as<double>(),2) * L;
 //assert(res > 0.1);

 printf("res vs. res2 %g %g\n", res, res2); 
 if (res < 0.2)
     res = L * 2.5;

 //res = res2;
 assert(BPTF);
 if (nuAlpha > 0){
       double nuAlpha_ =nuAlpha+ L;
       mat iWalpha_(1,1);
       iWalpha_.set(0,0, iWalpha + res);
       mat iiWalpha_ = zeros(1,1);
       iiWalpha_ = inv(iWalpha_);
       alpha = wishrnd(iiWalpha_, nuAlpha_).get(0,0);
       assert(alpha != 0);

       if (debug)
          cout<<"Sampling from alpha" <<nuAlpha_<<" "<<iiWalpha_<<" "<<alpha<<endl;
       printf("sampled alpha is %g\n", alpha); 
 }
}


mat calc_MMT(int nodes, vec &Umean, gl_types::iscope & scope){

   int batchSize = 1000;
   mat U = zeros(batchSize, D);
   mat MMT = zeros(D,D);
   int cnt = 0;
   timer t;

   std::vector<edge_id_t> tedges;
   if (nodes == M)
 	tedges= scope.out_edge_ids();
   else tedges = scope.in_edge_ids();

   assert(tedges.size() == nodes);

   int i=0;
   foreach(graphlab::edge_id_t oedgeid, tedges) {

         if (i % batchSize == 0){
            U.zeros();
            cnt = 1;
         }

         vertex_data * data= &scope.neighbor_vertex_data(nodes == M ? scope.target(oedgeid) : scope.source(oedgeid));
         vec mean = *vec2vec(data->pvec, D);
         Umean += mean;
          //Q.set_col(k, ret);
         t.start(); 
          for (int s=0; s<D; s++)
              U(i%batchSize,s)=mean(s);
     	 if (debug && (i==0 || i == nodes-1))
              std::cout<<" clmn "<<(i+((nodes==N)?M:0))<< " vec: " << mean <<std::endl;

         if ((cnt  == batchSize) || (cnt < batchSize && i == nodes-1)){
             MMT = MMT+transpose(U)*U;
	 }
         counter[8] += t.current_time();
         cnt++;
         i++;
   }
   Umean /= (nodes);
   if (debug)
      cout<<"mean: "<<Umean<<endl;

   assert(MMT.rows() == D && MMT.cols() == D);
   assert(Umean.size() == D);
    return MMT;
}


void sample_U(gl_types::iscope &scope, gl_types::icallback &scheduler, gl_types::ishared_data* shared_data) {
 
   if (debug)
     printf("Entering sample_U node\n");

   assert(BPTF);
   assert((int)scope.out_edge_ids().size() == M);
 

   vec Umean = zeros(D);
   mat UUT = calc_MMT(M,Umean,scope);
  
   double beta0_ = beta0 + M;
   vec mu0_ = zeros(D);
   mu0_ = (beta0*mu0 + M*Umean)/beta0_;
   double nu0_ = nu0 +M;
   vec dMu = mu0 - Umean;
   if (debug)
     cout<<"dMu:"<<dMu<<"beta0: "<<beta0<<" beta0_ "<<beta0_<<" nu0_ " <<nu0_<<" mu0_ " << mu0_<<endl;
   mat UmeanT = zeros(D,D);
   UmeanT = M*(itpp::outer_product(Umean, Umean));
   assert(UmeanT.rows() == D && UmeanT.cols() == D);
   mat dMuT =  zeros(D,D);
   dMuT = (beta0*M/beta0_)*(itpp::outer_product(dMu, dMu));
  mat iW0_ = iW0 + UUT - UmeanT + dMuT;
   mat W0_; 
   bool ret =inv(iW0_, W0_);
   assert(ret);
   mat tmp = (W0_+transpose(W0_))*0.5;
   if (debug)
      cout<<iW0<<UUT<<UmeanT<<dMuT<<W0_<<tmp<<nu0_<<endl;
    A_U = wishrnd(tmp, nu0_);
   mat tmp2;  
   ret =  inv(beta0_ * A_U, tmp2);
   assert(ret);
   mu_U = mvnrndex(mu0_, tmp2, D);
   if (debug)
          cout<<"Sampling from U" <<A_U<<" "<<mu_U<<" "<<Umean<<" "<<W0_<<tmp<<endl;
}

/*int _iiter(){
   return ((round_robin_scheduler<graph_type>*)&engine->get_scheduler())->get_iterations();   
}*/

void sample_V(gl_types::iscope &scope, gl_types::icallback &scheduler, gl_types::ishared_data* shared_data) {

   if (debug)
     printf("Entering sample_V node\n");


   assert(BPTF);
   vec Vmean = zeros(D);
   mat VVT = calc_MMT(N, Vmean, scope);   

   assert((int)scope.in_edge_ids().size() == N);

   double beta0_ = beta0 + N;
   vec mu0_ = (beta0*mu0 + N*Vmean)/beta0_;
   double nu0_ = nu0 +N;
   vec dMu = mu0 - Vmean;
   if (debug)
   cout<<"dMu:"<<dMu<<"beta0: "<<beta0<<" beta0_ "<<beta0_<<" nu0_ " <<nu0_<<endl;
   mat VmeanT = zeros(D,D);
   VmeanT = N*(itpp::outer_product(Vmean, Vmean));
   assert(VmeanT.rows() == D && VmeanT.cols() == D);
   mat dMuT = zeros(D,D);
   dMuT =  (beta0*N/beta0_)*itpp::outer_product(dMu, dMu);
   mat iW0_ = iW0 + VVT - VmeanT + dMuT;
   mat W0_=zeros(D,D); 
   bool ret = inv(iW0_, W0_);
   assert(ret);
   mat tmp = (W0_+transpose(W0_))*0.5;
  if (debug)
      cout<<"iW0: "<<iW0<<" VVT: "<<VVT<<" VmeanT: "<<VmeanT<<" dMuT: " <<dMuT<<"W0_"<< W0_<<" tmp: " << tmp<<" nu0_: "<<nu0_<<endl;
    A_V = wishrnd(tmp, nu0_);
   mat tmp2; 
   ret = inv(beta0_*A_V, tmp2);
   assert(ret);
   mu_V = mvnrndex(mu0_, tmp2, D);
 if (debug)
          cout<<"Sampling from V: A_V" <<A_V<<" mu_V: "<<mu_V<<" Vmean: "<<Vmean<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;
}


mat calc_DT(vertex_data * ptimes){

  assert(tensor);

  mat T = zeros(D, K);
  for (int i=0; i<K; i++){
     vec tmp = *vec2vec(ptimes[i].pvec, D);
     T.set_col(i,tmp);
  }
  
  mat diff = zeros(D,K-1);
  for (int i=0; i<K-1; i++){
     diff.set_col(i , T.get_col(i) - T.get_col(i+1));
  }
  if (debug)
     cout<<"T:"<<T<<" diff: " << diff<<endl;
  
  return diff;

}

void sample_T(gl_types::iscope &scope, gl_types::icallback &scheduler, gl_types::ishared_data* shared_data) {
   assert(BPTF);
   assert(tensor);
   assert((int)scope.vertex() == SAMPLE_T_NODE_OFFSET);
 
  if (debug)
     printf("Entering sample_T node\n");


  
   vertex_data ptimes[K];
   std::vector<edge_id_t> outs = scope.out_edge_ids();
   assert((int)outs.size() == K);

   int c = 0;
   foreach(graphlab::edge_id_t oedgeid, outs) {
       const vertex_data & pdata = scope.neighbor_vertex_data(scope.target(oedgeid));
       assert(scope.target(oedgeid) >= T_NODE_OFFSET);
       memcpy(ptimes[c].pvec, pdata.pvec, D*sizeof(double));
       c++;
   }

   assert(c == K);

   double beta0_ = beta0 + 1;
   vec pvec = *vec2vec(ptimes[0].pvec,D);
   vec mu0_ = (pvec + beta0*mu0T)/beta0_;
   double nu0_ = nu0 +K;
   //vec dMu = mu0 - Umean;
   if (debug){
     cout<<"beta0_ " << beta0_ << " beta0: " << beta0 << " nu0_ " << nu0_ << endl;
   } 

   mat dT = calc_DT(ptimes);
   vec dTe = pvec - mu0T;
   mat iW0_ = iW0T + dT*transpose(dT) + (beta0/beta0_)*(itpp::outer_product(dTe,dTe));
  
   mat W0_;
   bool ret =inv(iW0_, W0_);
   assert(ret);
   mat tmp = W0_+transpose(W0_)*0.5;
   A_T = wishrnd(tmp, nu0_);

   mat tmp2 ;
   ret = inv(beta0_*A_T, tmp2);
   assert(ret);
   mu_T = mvnrndex(mu0_, tmp2, D);
   if (debug)
          cout<<"Sampling from T: A_T" <<A_T<<" mu_V: "<<mu_T<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;
   
}



void T_update_function(gl_types::iscope &scope, gl_types::icallback &scheduler, gl_types::ishared_data* shared_data) {

     assert(tensor);

     int i = scope.vertex() - M-N;
     assert(i >=0 && i < K);	
    
     if (debug && (i == 0 || i == K-1)){
       printf("entering time node  %d \n", i);   
     } 
     if (K == 1)
       return;

     bool debug = false;

     //std::vector<edge_id_t> tedges = scope.out_edge_ids();
     //tedges.insert(tedges.end(), scope.in_edge_ids().begin(), scope.in_edge_ids().end());
  //for (int i=0; i< K; i++){
     assert(i >=0 && i < K);
     //assert(tedges.size() > edges[i].size());
     if (debug /*&& (i==0 || i == K-1)*/)
         printf("node %d with Q size: %d %d\n", i, (int)edges[i].size(), (int)D);

     int k=0;
     int batchSize = 100;
     mat Q(batchSize, D); 
     vec vals(batchSize);
     timer t;
     t.start();
     mat QQ = zeros(D,D);
     vec RQ = zeros(D);
     counter[4] += t.current_time();
     int cnt =0;

     vertex_data ptimes[K];

     foreach (edge_id_t edge, edges[i]){
     

         if (k % batchSize == 0){
            Q.zeros();
            vals.zeros();
            cnt = 1;
         }

         //find the right edge which matches the current time
         multiple_edges * medges= &scope.edge_data(edge);
         edge_data * data;
         for (int j=0; j< (int)medges->medges.size(); j++){
             data = &medges->medges[j];
             if (data->time == i)
                  break;
         }

         assert(data->time == i);
         assert((int)scope.target(edge)>= M);
         assert((int)scope.source(edge)< M);
         vertex_data * v1 = &scope.neighbor_vertex_data(scope.target(edge));
         vertex_data * v2 = &scope.neighbor_vertex_data(scope.source(edge));
 	 
         vec ret = dot(v1->pvec, v2->pvec, D);
     
          //Q.set_col(k, ret);
         t.start(); 
          for (int s=0; s<D; s++)
              Q(k%batchSize,s)=ret(s);

     	 if (debug && (i==0 || i == K-1) && (k == 0 || k == (int)edges[i].size() - 1))
              std::cout<<" clmn "<<k<< " vec: " << ret<<std::endl;

         assert(data->weight >0);
         vals[k%batchSize] = data->weight;
         
         k++;

         if ((cnt  == batchSize) || (cnt < batchSize && k == (int)edges[i].size()-1)){
             QQ += transpose(Q)*Q;
             RQ += transpose(Q)*vals;
             assert(QQ.rows() == D && QQ.cols() == D);
	 }
         counter[5] += t.current_time();
         cnt++;

         
     }

    if (debug  && (i == 0 ||i== K-1 )){
        std::cout<<"QQ:"<<QQ<<std::endl;
        std::cout<<"RQ:"<<RQ<<std::endl;
     }

     int found = 0;  
     foreach (edge_id_t oedgeid, scope.out_edge_ids()){
          int other = scope.target(oedgeid);
          if (other >= T_NODE_OFFSET && other < SAMPLE_U_NODE_OFFSET){   
             vertex_data &vdata = scope.neighbor_vertex_data(scope.target(oedgeid));
	     memcpy(ptimes[other-(T_NODE_OFFSET)].pvec, vdata.pvec, D*sizeof(double)); 
             found++; 
          }
     }
     foreach (edge_id_t iedgeid, scope.in_edge_ids()){
          int other = scope.source(iedgeid);
          if (other >= T_NODE_OFFSET && other < SAMPLE_U_NODE_OFFSET){   
             vertex_data &vdata = scope.neighbor_vertex_data(scope.source(iedgeid));
	     memcpy(ptimes[other-(T_NODE_OFFSET)].pvec, vdata.pvec, D*sizeof(double)); 
             found++; 
          }
     }

    if (BPTF && tensor)
        if (i == 0 || i == K-1)
         	assert(found == 1);
    	else assert(found == 2); 
    
     assert(RQ.size() == D);
    assert((unsigned int)k == edges[i].size());
    vec sol;    
    vec out;

    t.start();
    if (i == 0){
      vec t1 = *vec2vec(ptimes[1].pvec, D);
      if (!BPTF){
        QQ = QQ+ 2*eDT;
        bool ret = itpp::ls_solve(QQ, pT*(t1 + vones*muT)+ RQ, out);
        assert(ret);
      }
      else {
        mat A_k = 2*A_T+alpha*QQ;
        mat iAk_;
        bool ret = inv(A_k, iAk_);
        assert(ret);
        vec muk_ = iAk_*(A_T*(t1 + mu_T)+alpha*RQ); //TODO
        out = mvnrndex(muk_, iAk_, D);
      }
    }
    else if (i == K-1){
      vec tk_2 = *vec2vec(ptimes[K-2].pvec, D);
      if (!BPTF){
        QQ = QQ + eDT;
        bool ret = itpp::ls_solve(QQ, pT*tk_2 + RQ, out);
        assert(ret); 
      }
      else {
        mat A_k = A_T+alpha*QQ;
        mat iAk_; 
        bool ret = inv(A_k, iAk_);
        assert(ret);
        vec muk_ = iAk_*(A_T*tk_2+alpha*RQ); //TODO
        out = mvnrndex(muk_, iAk_, D);
      }
    }
    else {
      vec tsum = *vec2vec(ptimes[i-1].pvec,D)+ *vec2vec(ptimes[i+1].pvec, D);
      if (!BPTF){
        QQ = QQ + 2*eDT;
        bool ret = itpp::ls_solve(QQ, pT*tsum + RQ, out);
        assert(ret);
      }
      else {
        mat A_k = 2*A_T+alpha*QQ;
        mat iAk_;
        bool ret = inv(A_k, iAk_);
        assert(ret);
        vec muk_ = iAk_*(A_T* tsum +alpha*RQ); //TODO
        out = mvnrndex(muk_, iAk_, D);
      }
    }

    //vec2vec(ptimes[i].pvec, out, D);
   
     vertex_data & data = scope.vertex_data();
     vec2vec(data.pvec, out, D);
     if (debug && (i == 0|| i == K-1)) 
           std::cout<<"data.pvec"<<*vec2vec(data.pvec,D)<<std::endl;
   
     assert(QQ.rows() == D && QQ.cols() == D);
     counter[6] += t.current_time();

     if (i == K-1){
          last_iter();
     }

}



double calc_rmse(graph_type * _g, bool test, double & res, bool fast){

    if (test && Le == 0)
       return NAN;
   
    fast = false;
 
     res = 0;
     double RMSE = 0;
     int e = 0;
     for (int i=0; i< M+N; i++){ //TODO: optimize to start from N?
       vertex_data * data = &g.vertex_data(i);
       foreach(edge_id_t iedgeid, _g->in_edge_ids(i)) {
       
         multiple_edges & edges = _g->edge_data(iedgeid);
         if (edges.control)
             continue;

         vertex_data * pdata = &g.vertex_data(_g->source(iedgeid)); 
         for (int j=0; j< (int)edges.medges.size(); j++){       

                 edge_data & edge = edges.medges[j];
           assert(edge.weight != 0);
         
           vertex_data * ptimes = NULL;
          if (tensor)
               ptimes = &g.vertex_data((int)edge.time +T_NODE_OFFSET);
 
            double sum = 0; 
           double add = gl::rmse(data->pvec, pdata->pvec, tensor? ptimes->pvec:NULL, D, edge.weight, sum);
           assert(sum != 0);         
           if (BPTF && iiter > BURN_IN)
              edge.avgprd += sum;

           if (debug && (i== M || i == M+N-1) && (e == 0 || e == (test?Le:L)))
              cout<<"RMSE:"<<i <<"u1"<< *vec2vec(data->pvec, D) << " v1 "<< *vec2vec(pdata->pvec, D)<<endl; 
          //assert(add<25 && add>= 0);
       
           if (BPTF && iiter > BURN_IN){
             add = pow((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
           }
         
           RMSE+= add;
           e++;
       }
     }
   }
  res = RMSE;
  if (!fast)
     assert(e == (test?Le:L));
  return sqrt(RMSE/(double)e);

}

void parse_edge(edge_data& edge, vertex_data & pdata, mat & Q, vec & vals, int i, vertex_data* ptime){
        
         double buf[D];  
         double * pbuf = buf; 

         assert(edge.weight != 0);

         if (tensor){
         	gl::dot2(pdata.pvec,  ptime[(int)edge.time].pvec, buf, D);  
         }
	 else {
                pbuf = pdata.pvec;
         }
 
         for (int j=0; j<D; j++){
         	Q.set(j,i, pbuf[j]); 
         }
         vals[i] = edge.weight;
}
 

int count_edges(const std::vector<edge_id_t> &es){
   int cnt = 0; 
   for (int j=0; j< (int)es.size(); j++){
       cnt += g.edge_data(es[j]).medges.size();
    }
   return cnt;
}


/***
 * UPDATE FUNCTION
 */
void pmf_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler,
                         gl_types::ishared_data* shared_data) {
    

  bool debug = false;
  /* GET current vertex data */

  vertex_data& vdata = scope.vertex_data();
 
  
  /* CALCULATE new value */
  if (debug&& (scope.vertex() == 0 || (int)scope.vertex() == M-1 || (int)scope.vertex() == M || (int)scope.vertex() == M+N-1)){
    printf("entering %s node  %u \n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex());   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , vdata.pvec, D);
  }

  assert((int)scope.vertex() < M+N);

      vdata.rmse = 0;

 
     //int i=0; 
     // update neighbors 
    std::vector<edge_id_t> outs = scope.out_edge_ids();
    std::vector<edge_id_t> ins = scope.in_edge_ids();
    int numedges = vdata.num_edges;

    if (numedges == 0){
            return;
    }

    vertex_data ptimes[K];

    //assert(outs.size() == 0 || ins.size() == 0);   
    assert(numedges > 0);

    timer t;
    t.start(); 
    mat Q(D,numedges);
    vec vals(numedges);
    counter[0] += t.current_time();

    int i=0;

    if ((int)scope.vertex() < M){
    
     //MOVIE NODES
      bool found = false;  
      foreach(graphlab::edge_id_t oedgeid, ins) {

           int other = scope.source(oedgeid);
           if (other <M+N || other == SAMPLE_U_NODE_OFFSET)
 		continue;

           const vertex_data & pdata = scope.neighbor_vertex_data(other); 
	   memcpy(ptimes[other-M-N].pvec, pdata.pvec, D*sizeof(double));
 	      //  scheduler.add_task(gl_types::update_task(edge.time+M+N,

              //                         T_update_function),1);
              found = true;
         }   
       if (BPTF&&tensor) assert(found);

       foreach(graphlab::edge_id_t oedgeid, outs) {


           multiple_edges &medges =scope.edge_data(oedgeid);
           vertex_data  pdata = scope.neighbor_vertex_data(scope.target(oedgeid)); 
           assert(scope.target(oedgeid) >= M && scope.target(oedgeid)<M+N);
	   
           for (int j=0; j< (int)medges.medges.size(); j++){
              edge_data& edge = medges.medges[j];
              parse_edge(edge, pdata, Q, vals, i, ptimes); 
     
              if (debug && ((int)scope.vertex() == 0 || (int)scope.vertex() == M-1) && (i==0 || i == numedges-1))
                 std::cout<<"set col: "<<i<<" " <<Q.get_col(i)<<" " <<std::endl;
              i++;
              //if (tensor) 
 	      //  scheduler.add_task(gl_types::update_task(edge.time+M+N,
              //                         T_update_function),1);
         }   
           
          // scheduler.add_task(gl_types::update_task(scope.target(oedgeid),
          //                             pmf_update_function),1);
         
            
     }
   }

   else {

      bool found = false;
      foreach(graphlab::edge_id_t iedgeid, outs) {

           int other = scope.target(iedgeid);
           if (other <M+N || other == SAMPLE_V_NODE_OFFSET)
 		continue;

           const vertex_data & pdata = scope.neighbor_vertex_data(other); 
	   memcpy(ptimes[other-M-N].pvec, pdata.pvec, D*sizeof(double));
 	      //  scheduler.add_task(gl_types::update_task(edge.time+M+N,
              //                         T_update_function),1);
           found = true;   
      }   
     if (BPTF && tensor) assert(found);

     //USER NODES
     foreach(edge_id_t iedgeid, ins) {

         vertex_data & pdata = scope.neighbor_vertex_data(scope.source(iedgeid)); 
         assert(scope.source(iedgeid) <M);
	 multiple_edges & medges =scope.edge_data(iedgeid);
	 for (int j=0; j< (int)medges.medges.size(); j++){
              edge_data& edge = medges.medges[j];
              parse_edge(edge, pdata, Q, vals, i, ptimes); 
     
              if (debug && ((int)scope.vertex() == M || (int)scope.vertex() == M+N-1) && (i==0 || i == numedges-1))
                 std::cout<<"set col: "<<i<<" " <<Q.get_col(i)<<" " <<std::endl;

              i++;
              //if (tensor) scheduler.add_task(gl_types::update_task(edge.time+M+N,
              //                         T_update_function),1);
              double sum;     
              double trmse = gl::rmse(vdata.pvec, pdata.pvec, tensor?ptimes[(int)edge.time].pvec:NULL, D, edge.weight, sum);
              assert(sum != 0);
              if (BPTF && iiter > BURN_IN){
              //edge.avgprd += sum;       
                 trmse = pow((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
              }
              vdata.rmse += trmse; 
 
     
           }
      
            //scheduler.add_task(gl_types::update_task(scope.target(iedgeid),
            //                           pmf_update_function),1);
           //vertex_data->rmse += gl::rmse(vertex_data->pvec, pdata->pvec, D, edge->weight);
     }

    }
      assert(i == numedges);
     

      vec result;
      
     if (!BPTF){
        t.start();
        bool ret = itpp::ls_solve(Q*itpp::transpose(Q)+eDT, Q*vals, result);
        if (debug && scope.vertex() == 0)
            cout<<"Q:"<<Q<<endl;
        //assert(result.size() == D);     
        assert(ret);
        counter[3] += t.current_time();
        }
    else {
        assert(Q.rows() == D);
        t.start();
        mat iAi_;
        bool ret =inv(((int)scope.vertex() < M? A_U : A_V) + alpha *  Q*itpp::transpose(Q), iAi_);
        counter[10]+= t.current_time();
        assert(ret);
        t.start();
        vec mui_ =  iAi_*((((int)scope.vertex() < M)? (A_U*mu_U) : (A_V*mu_V)) + alpha * Q * vals);
        counter[11]+= t.current_time();
       
         t.start();
        result = mvnrndex(mui_, iAi_, D); 
        assert(result.size() == D);
        counter[9] += t.current_time();
    }

    if (debug && ((int)scope.vertex() < 2 || (int)scope.vertex() == M-1 || (int)scope.vertex() == M || (int)scope.vertex() == M+N-1)){
	  std::cout <<(BPTF?"BPTF":"ALS")<<" Q: " << Q << " result: " << result << " edges: " << i << std::endl;
    }
      
     for (i=0; i<D; i++){
         vdata.pvec[i] =result[i];
         //if (record_history) 
         //   vertex_data->agg[i] = 0.8*vertex_data->agg[i] + 0.2*result[i];
     }

    if (!tensor && (int)scope.vertex() == M+N-1)
          last_iter();

}

void last_iter(){
/*         printf("Entering last iter with %d\n",
 	 	((round_robin_scheduler<graph_type>*)&engine->get_scheduler())->get_iterations());   */

	 double res,res2;
         double rmse = calc_rmse(&g, false, res, false);
         //rmse=0;
         printf("%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f TEST RMSE=%0.4f.\n", gt.current_time(), BPTF?"BPTF":"ALS", iiter,calc_obj(/*res*/),  rmse, calc_rmse(&g1, true, res2, true));
         iiter++;
         if (iiter == BURN_IN && BPTF){
                printf("Finished burn-in period. starting to aggregate samples\n");
          }
         
         if (iiter == 20){
 		engine->stop();
          }
          if (BPTF){
     		sample_alpha(res);
     		//sample_U();
     		//sample_V();
     		//if (tensor) 
                //    sample_T();
  	  }
}




double calc_obj(/*double res*/){
   
     double res = powf(sdm.get(RMSE).as<double>(),2) * L;

     double sumU = 0, sumV = 0, sumT = 0;
     timer t;
     t.start(); 
     for (int i=0; i< M; i++){
	vertex_data * data = &g.vertex_data(i);
 	sumU += gl::square_sum(data->pvec, D);
     } 

    for (int i=M; i< M+N; i++){
	vertex_data * data = &g.vertex_data(i);
 	sumV += gl::square_sum(data->pvec, D);
     } 


   mat T;
   if (tensor){
      T = zeros(D,K);
      for (int i=0; i<K; i++){
        vertex_data & data = g.vertex_data(i+T_NODE_OFFSET);
        sumT += pow(norm(*vec2vec(data.pvec,D) - vones, 2),2);
        vec tmp = *vec2vec(data.pvec,D);
        T.set_col(i, tmp);
     }
    sumT *= pT;
   }
   counter[7]+= t.current_time();
  
   double obj = (res + pU*sumU + pV*sumV + sumT + (tensor?trace(T*dp*T.transpose()):0)) / 2.0;
   return obj;
}






 
/** 
 * ==== SETUP AND START
 */
void start(int argc, char ** argv) {
      
  
  printf("loading data file %s\n", infile.c_str());
  load_pmf_graph(infile.c_str(), &g, false);

  printf("loading data file %s\n", (infile+"e").c_str());
  load_pmf_graph((infile+"e").c_str(),&g1, true);

   engine_options o;
    o.parse_command_line(argc, argv);
    engine = o.create_engine(g);
    engine->set_shared_data_manager(&sdm);
   //sdm.sync(PRIMAL_LOSS);   
   sdm.set_sync(RMSE,
                gl_types::sync_ops::sum<double, get_rmse>,
                apply_func,
               double(0),  20000); 

  dp = GenDiffMat(K)*pT;
  if (debug)
  	std::cout<<dp<<std::endl;

  /**** CREATE INITIAL TASKS ******/
  std::vector<vertex_id_t> um;
  std::vector<vertex_id_t> tv;
  
  for (int i=0; i< M+N; i++)
     um.push_back(i);
  
   engine->add_tasks(um, pmf_update_function, 1);
  
  if (tensor){
    for (int i=M+N; i< M+N+K; i++)
       tv.push_back(i);
    engine->add_tasks(tv, T_update_function, 1);
    
    tv.clear(); tv.push_back(SAMPLE_U_NODE_OFFSET);
    engine->add_tasks(tv, sample_U,1);
    tv.clear(); tv.push_back(SAMPLE_V_NODE_OFFSET);
    engine->add_tasks(tv, sample_V,1);
    tv.clear(); tv.push_back(SAMPLE_T_NODE_OFFSET);
    engine->add_tasks(tv, sample_T,1);

    assert(g.num_vertices() == M+N+K+3);

  }


  size_t start_node = 0;
  if (BPTF && tensor)
      start_node = SAMPLE_U_NODE_OFFSET;

  engine->set_options("start_vertex", start_node);
  
  printf("%s for %s (%d, %d, %d):%d.  D=%d\n", BPTF?"BPTF":"PTF_ALS", tensor?"tensor":"matrix", M, N, K, L, D);
  
  if (!BPTF)
     printf("pU=%g, pV=%g, pT=%g, muT=%g, D=%d\n", pU, pV, pT, muT,D);  

  //if (BPTF)
    init_self_pot(); //init anyway
   init_pmf();

  double /*res,*/ res2;
  //double rmse =  calc_rmse(&g, false, res, false);
  printf("complete. Obj=%g, TEST RMSE=%0.4f.\n", calc_obj(), calc_rmse(&g1, true, res2, false));


   if (BPTF){
     sample_alpha(1*L);
     //sample_U();
     //sample_V();
     //if (tensor) 
     //   sample_T();
  }

  /// Timing
  gt.start();

  // assert(false); // Code not checked in breaks build
  g.set_finalized();
  
  //g1.set_finalized();
  
  /**** START GRAPHLAB *****/
  engine->start();

  //rmse =  calc_rmse(&g, false, res, false);
  printf("Final result. Obj=%g, TEST RMSE= %0.4f.\n", calc_obj(),  calc_rmse(&g1, true, res2, false));
  
  sdm.sync(RMSE);
  printf("SDM rmse is: %g\n", sqrt(sdm.get(RMSE).as<double>() / L));


  /**** POST-PROCESSING *****/
  double runtime = gt.current_time();
  printf("Finished in %lf \n", runtime);
      
  //timing counters
  for (int i=0; i<20; i++){
    printf("Counters are: %d) %g\n",i, counter[i]); 
  }
     
}

   template<typename edgedata>
  int read_mult_edges(FILE * f, int nodes, graph_type * g, bool test){
     

     //typedef typename graph::edge_data_type edge_data;
     bool * flags = NULL;
     if (options == BPTF_TENSOR_MULT || options == ALS_TENSOR_MULT){
         flags = new bool[nodes];
         memset(flags, 0, sizeof(bool)*nodes);
     }
 
     unsigned int e;
     int rc = fread(&e,1,4,f);
     assert(rc == 4);
     printf("Creating %d edges...\n", e);
     assert(e>0);
     int total = 0;
     edgedata* ed = new edgedata[200000];
     int edgecount_in_file = e;
     while(true){
        //memset(ed, 0, 200000*sizeof(edata3));
       rc = (int)fread(ed, sizeof(edgedata), _min(200000, edgecount_in_file - total), f);
       total += rc;

      for (int i=0; i<rc; i++){
        multiple_edges edges;
        edge_data edge;
     	assert(ed[i].weight != 0); // && ed[i].weight <= 5);
    	assert((int)ed[i].from >= 1 && (int)ed[i].from <= nodes);
    	assert((int)ed[i].to >= 1 && (int)ed[i].to <= nodes);
    	assert((int)ed[i].to != (int)ed[i].from);
        edge.weight = (double)ed[i].weight;
        edge.time = (double)ed[i].time - 1;
        if (tensor && !test && BPTF){
          multiple_edges medges;
          medges.control = true;
          medges.medges.push_back(edge);
          if ((int)ed[i].from - 1 < M){
          	g->add_edge((int)ed[i].time - 1 + T_NODE_OFFSET, (int)ed[i].from -1, medges);    	
          	g->add_edge((int)ed[i].to - 1, (int)ed[i].time -1 + T_NODE_OFFSET, medges);    	
          }
          else {
                g->add_edge((int)ed[i].from - 1, (int)ed[i].time -1 + T_NODE_OFFSET, medges);    	
                g->add_edge((int)ed[i].time - 1 + T_NODE_OFFSET, (int)ed[i].to -1, medges);    	
          }
        }
         std::pair<bool, edge_id_t> ret;
         if (options != BPTF_TENSOR_MULT && options != ALS_TENSOR_MULT){//no support for specific edge returning on different times
              ret.first = false;
         }
         else if (flags[(int)ed[i].from-1] == true && flags[(int)ed[i].to-1] == true){
              ret = g->find((int)ed[i].from-1, (int)ed[i].to-1);
         }
         else ret.first = false;

         if (ret.first == false){
                edges.medges.push_back(edge); 
        	g->add_edge((int)ed[i].from-1, (int)ed[i].to-1, edges); // Matlab export has ids starting from 1, ours start from 0
                if (options == BPTF_TENSOR_MULT || options == ALS_TENSOR_MULT){
                   flags[(int)ed[i].from-1] = true;
                   flags[(int)ed[i].to-1] = true;
                }
         }
         else {
              g->edge_data(ret.second).medges.push_back(edge);
         }
	} 
        printf(".");
        fflush(0);
        if (rc == 0 || total >= edgecount_in_file)
	   break;
    }
    assert(total == (int)e);
    delete [] ed; ed = NULL;
    if (flags != NULL)
        delete[] flags;
  
    if (!test && tensor && BPTF){
         assert(g->num_vertices() == M+N+K+3);

         edge_data ed;
         multiple_edges medges;
         medges.control = true;
         medges.medges.push_back(ed);
         // SAMPLE T TO ALL TIMES
         for (int i=0; i<K; i++){
              assert(g->out_edge_ids(i+T_NODE_OFFSET).size() > 0 || g->in_edge_ids(i+T_NODE_OFFSET).size() > 0);
              g->add_edge(SAMPLE_T_NODE_OFFSET, T_NODE_OFFSET+i, medges);
         }
         // SAMPLE U to all U
         for (int i=0; i<M; i++){
              g->add_edge(SAMPLE_U_NODE_OFFSET, i, medges);
         }
         // V to all SMAPLE_V
         for (int i=M; i<N+M; i++){
              g->add_edge(i, SAMPLE_V_NODE_OFFSET, medges);
         }
         // Chain of T's
         for (int i=0; i<K-1; i++){
              g->add_edge(T_NODE_OFFSET+i, T_NODE_OFFSET+i+1, medges);
         }

     } 

    return e;
  }



void load_pmf_graph(const char* filename, graph_type * g, bool test) {

  printf("Loading %s %s\n", filename, test?"test":"train");
  FILE * f = fopen(filename, "r");
  if (test && f == NULL){
     printf("skipping test data\n");
     return;
  }

  assert(f!= NULL);

  fread(&M,1,4,f);//movies
  fread(&N,1,4,f);//users/
  fread(&K,1,4,f);//time
  assert(K>=1);
  assert(M>=1 && N>=1); 

 
  vertex_data vdata;
  gl::ones(vdata.pvec, D, 0.1);

  //ADD MOVE NODES
  for (int i=0; i<M; i++){
 	//gl::rand(vdata.pvec, D);%TODO
    // g->add_vertex(vdata);
        g->add_vertex(vdata);
        if (debug && (i<= 5 || i == M-1))
            debug_print_vec("U: ", vdata.pvec, D);
 }
   
  //ADD USER NODES
  for (int i=0; i<N; i++){
 	//gl::rand(vdata.pvec, D);
       g->add_vertex(vdata);
        if (debug && (i<=5 || i==N-1))
            debug_print_vec("V: ", vdata.pvec, D);
}
  
 //ADD TIME NODES
 if (!test && tensor){
  //init times
  vertex_data * times = new vertex_data[K];
  vec tones = itpp::ones(D)*(K==1?1:0.1);
  for (int i=0; i<K; i++){
  	 vec2vec(times[i].pvec ,tones, D);
         g->add_vertex(times[i]);
         if (debug && (i <= 5 || i == K-1))
		debug_print_vec("T: ", times[i].pvec, D);
  }
  }

  //ADD SAMPLE U,V,T nodes
  if (!test & tensor && BPTF){
       vertex_data vdata;
       for (int i=0; i<3; i++)
           g->add_vertex(vdata);
  }

  int val = read_mult_edges<edata2>(f, M+N, g, test);
  if (!test)
    L = val;
  else Le = val;

  if (!test && tensor && K>1) 
  	edges = new std::vector<edge_id_t>[K]();

  bool normalize = false;
  double normconst = 1;
  if (!strcmp(filename, "netflow") || !strcmp(filename, "netflowe")){
        normconst = 138088872;
        normalize = true;
  }
        

 //verify edges
  for (int i=M; i < M+N; i++){
      foreach(graphlab::edge_id_t eid, g->in_edge_ids(i)){          
          multiple_edges * tedges= &g->edge_data(eid);
          int from = g->source(eid);
          int to = g->target(eid);
          assert(from < M);
          assert(to >= M && to < M+N);

          for (int j=0; j< (int)tedges->medges.size(); j++){
              edge_data * data= &tedges->medges[j];
              assert(data->weight != 0);  
              if (normalize)
                data->weight /= normconst;    
    //assert(data->weight == 1 || data->weight == 2 || data->weight == 3 || data->weight ==4 || data->weight == 5);
                assert(data->time < K);
  
              if (K > 1 && !test && tensor)
	  	   edges[(int)data->time].push_back(eid);
              }
          }
  }
     
  if (!test){
     for (int i=0; i<M+N; i++){
        vertex_data &vdata = g->vertex_data(i);
        if (i < M)
        	vdata.num_edges = count_edges(g->out_edge_ids(i));
        else
                vdata.num_edges = count_edges(g->in_edge_ids(i));
     }
  }
   
  if (!test && tensor && K>1){
    int cnt = 0;
    for (int i=0; i<K; i++){
       cnt+= edges[i].size();
    }
    assert(cnt == L);
  }
  fclose(f);
}


 
int main(int argc,  char *argv[]) {

  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);



  double weight = -1;

  infile = argv[1];

  if (infile == "" || argc <= 3) {
    std::cout << "PMF <intput file> <weight> <options>\n";
    return 0;
  }
  weight = atof(argv[2]);
   
  assert(weight>0);
  printf("setting regularization %e\n", weight);
  pV = pU = weight;  
  assert(pV > 0);

 options = atoi(argv[3]);
 	 printf("setting run mode %d\n", options);
         switch(options){
              case ALS_MATRIX:
                 tensor = false; BPTF = false;
                 break;
              case BPTF_TENSOR:
                 tensor = true; BPTF = true;
                 break;
              case BPTF_MATRIX:
                 tensor = false; BPTF = true;
                 break;
              case BPTF_TENSOR_MULT:
                 tensor = true; BPTF = true;
                 break;
              case ALS_TENSOR_MULT:
                 tensor = true; BPTF = false;
                 break;
               default:
                 assert(0);
         }

  logger(LOG_INFO,(BPTF?"BPTF starting\n": "PMF starting\n"));
   start(argc, argv);
}



#include <graphlab/macros_undef.hpp>
