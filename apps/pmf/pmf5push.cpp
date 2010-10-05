// written by Danny Bickson, Aapo Kyrola CMU 
#include <cstdio>
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>

#include <fstream>
#include <cmath>
#include <cstdio>
#include <cfloat>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab.hpp>
#include <graphlab/distributed/graph/distributed_graph.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_fullsweep_sdm.hpp>
#include <graphlab/distributed/pushy_distributed_engine.hpp>
#include <graphlab/distributed/distributed_round_robin_scheduler.hpp>


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

#define TIME_OFFSET 1
#define A_U_OFFSET 2
#define A_V_OFFSET 3
#define A_T_OFFSET 4
#define ALPHA_OFFSET 5

#define IS_ROUND_ROBIN_CONSTANT 6
#define MAX_ITERATIONS_CONSTANT 7

const int NUM_ITERATIONS_TO_RUN = 20;

bool BPTF = true;
#define D 30 //diemnsion for U,V
int options;
timer gt;
using namespace itpp;
using namespace std;
bool debug = false;

distributed_control * __dc;

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
//double alpha = 0;
double beta = 1;
double beta0 = 1; //TODO
mat W0;
mat W0T;
double iWalpha;
mat iW0;
mat iW0T;
//mat A_U, A_V, A_T;
//vec mu_U, mu_V, mu_T;

bool record_history = false;
int BURN_IN =10;
bool tensor = true;
double counter[20];

int myprocid;





struct command_line_options {
  size_t ncpus;
  bool prepart;
  std::string scope;
  std::string scheduler;
  std::string basefile;
};

bool parse_command_line(int argc, char** argv, ::command_line_options& opts) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Denoise a randomly generated image using Gibbs Sampling.");
  // Set the program options
  desc.add_options()
    ("ncpus",  boost_po::value<size_t>(&(opts.ncpus))->default_value(2),
     "Number of cpus to use.")
    ("prepartition",  boost_po::value<bool>(&(opts.prepart))->default_value(false),
     "prepartiton")
    ("scope",
     boost_po::value<std::string>(&(opts.scope))->default_value("edge"),
     "Options are {vertex, edge, full}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("round_robin"),
     "Options are multiqueue_fifo/round_robin")
    ("basefile",
     boost_po::value<std::string>(&(opts.basefile))->default_value(""),
     "Base input/output graphname");
  // Parse the arguments
  boost_po::variables_map vm;
  boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
  boost_po::notify(vm);
  if(vm.count("help")) {
    std::cout << desc << std::endl;
    return false;
  }
  return true;
} // end of parse command line arguments





/** Vertex and edge data types **/
struct vertex_data {
  double pvec[D];
  double rmse;
  int num_edges;
  int rounds;
  vertex_data(){
    memset(pvec, 0, sizeof(double)*D);
    rmse = 0;
    num_edges = 0;
    rounds = 0;
  }
  void save(graphlab::oarchive &oarc) const{
    serialize(oarc, this, sizeof(vertex_data));  
  }
  void load(graphlab::iarchive &iarc) {
    deserialize(iarc, this, sizeof(vertex_data));
  }
}; 

struct edge_data {
  double weight; 
  double time;
  double avgprd;
  edge_data(){ weight = 0; time = 0; avgprd = 0;}
  void save(graphlab::oarchive &oarc) const{
    serialize(oarc, this, sizeof(edge_data));  
  }
  void load(graphlab::iarchive &iarc) {
    deserialize(iarc, this, sizeof(edge_data));  
  }
};

struct multiple_edges{
  std::vector<edge_data> medges;
  void save(graphlab::oarchive &oarc) const{
    oarc << medges;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> medges;
  }
};

struct QQR{
  mat QQ;
  vec QR;
  QQR(){ QQ = zeros(D,D); QR = zeros(D); };
  void save(graphlab::oarchive &oarc) const{
    oarc << QQ;
    oarc << QR;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> QQ;
    iarc >> QR;
  }
};

struct mult_QQR{
  std::vector<QQR> vals;
  double rmse;
  mult_QQR(){ rmse = 0;  }
  void save(graphlab::oarchive &oarc) const{
    oarc << vals;
    oarc << rmse;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> vals;
    iarc >> rmse;
  }
  
};

struct mult_vec{
  std::vector<vec> vals;
  double rmse;
 
  mult_vec(){ rmse = 0; }
  void save(graphlab::oarchive &oarc) const{
    oarc << vals;
    oarc << rmse;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> vals;
    iarc >> rmse;
  }
};



void agg_MM(const vertex_data & pdata, mat & QQ, vec & QR);
void agg_QQ(const edge_data& edge, const vertex_data & pdata, const vertex_data & other, mat & QQ, vec & QR, int i);
typedef graphlab::graph<vertex_data, multiple_edges> graph_type;
typedef graphlab::distributed_graph<vertex_data, multiple_edges> graph_dtype;
typedef graphlab::types<graph_type> gl_types;
typedef graphlab::types<graph_dtype> gl_dtypes;
graph_type g;
graph_dtype * dg = NULL;
graph_type g1;
//gl_types::thread_shared_data sdm;


double get_rmse(const vertex_data & v){
  return v.rmse;
}


void sync_MMT(size_t index, const gl_dtypes::ishared_data& shared_data,
              gl_dtypes::iscope& scope, graphlab::any& new_data) {


  timer t;
  t.start();
  bool debug = false;
  /* GET current vertex data */

  
  const vertex_data& vdata = scope.const_vertex_data();
 
  QQR& newdata = new_data.as<QQR>();

  
  /* CALCULATE new value */
  if (debug&& ((int)scope.vertex() == 0 || (int)scope.vertex() == M-1 || (int)scope.vertex() == M || (int)scope.vertex() == M+N-1)){
    printf("MMT entering %s node  %d for time: %d\n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex(), (int)index);   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , vdata.pvec, D);
  }

  assert((int)scope.vertex() < M+N);


  agg_MM(vdata, newdata.QQ, newdata.QR);
            
  if (debug && ((int)scope.vertex() == 0 || (int)scope.vertex() == M-1 || (int)scope.vertex() == M || (int)scope.vertex() == M+N-1))
    std::cout<<"AG MM set col: "<<scope.vertex() <<" " <<newdata.QQ<<" " <<newdata.QR<<" " <<std::endl;

  new_data = newdata;
  ASSERT_GT(sumsum(newdata.QQ), 0);
  counter[3] += t.current_time();
}


void sync_QQR(size_t index, const gl_dtypes::ishared_data& shared_data,
              gl_dtypes::iscope& scope, graphlab::any& new_data) {


  timer t2;
  t2.start();
  bool debug = false;
  /* GET current vertex data */

  assert((int)index>=0 && (int)index<K);

  const vertex_data& vdata = scope.const_vertex_data();
 
  mult_QQR & newdata = new_data.as<mult_QQR>();
  
  //printf("Vertex id %d rmse = %lf / %lf\n", scope.vertex(), vdata.rmse, newdata.rmse);
  
  newdata.rmse += vdata.rmse;
  
  /* CALCULATE new value */
  if (debug&& ((int)scope.vertex() == 0 || (int)scope.vertex() == M-1 || (int)scope.vertex() == M || (int)scope.vertex() == M+N-1) && ((int)index ==0 || (int)index == K-1)){
    printf("QQR entering %s node  %u for time: %d\n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex(), (int)index);   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , vdata.pvec, D);
  }

  assert((int)scope.vertex() < M+N);

  //int i=0; 
  // update neighbors 
  edge_list outs = scope.out_edge_ids();
  edge_list ins = scope.in_edge_ids();
  int numedges = vdata.num_edges;

  if (numedges == 0){
    return;
  }

  assert(numedges > 0);

  timer t;
  t.start(); 
  counter[0] += t.current_time();

  int i=0;

  if ((int)scope.vertex() < M){
    
    //MOVIE NODES
    foreach(graphlab::edge_id_t oedgeid, outs) {

      const multiple_edges &medges =scope.const_edge_data(oedgeid);
      const vertex_data & pdata = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
     
      for (int j=0; j< (int)medges.medges.size(); j++){
        const edge_data& edge = medges.medges[j];
        //parse_edge(edge, pdata, Q, vals, i); 
        // if (edge.time != index)
        //      continue;     
        assert(edge.time >= 0 && edge.time < K);
        agg_QQ(edge, vdata, pdata, newdata.vals[(int)edge.time].QQ, newdata.vals[(int)edge.time].QR, i);
            
        if (debug && ((int)scope.vertex() == 0 || (int)scope.vertex() == M-1) && (i==0 || i == numedges-1) && ((int)index ==0 || (int)index == K-1))
          std::cout<<"set col: "<<i<<" " <<newdata.vals[(int)edge.time].QQ<<" " <<newdata.vals[(int)edge.time].QR<<" " <<std::endl;
        i++;
      }   
    }
  }

  else {

    //USER NODES
    foreach(edge_id_t iedgeid, ins) {

      const multiple_edges & medges =scope.const_edge_data(iedgeid);
      const vertex_data  &pdata = scope.const_neighbor_vertex_data(scope.source(iedgeid)); 
   
      for (int j=0; j< (int)medges.medges.size(); j++){
        const edge_data& edge = medges.medges[j];
        //if (edge.time != index)
        //    continue;    

        assert(edge.time >= 0 && edge.time < K);
        agg_QQ(edge, vdata, pdata, newdata.vals[(int)edge.time].QQ, newdata.vals[(int)edge.time].QR, i);
 
        if (debug && ((int)scope.vertex() == M || (int)scope.vertex() == M+N-1) && (i==0 || i == numedges-1) && ((int)index ==0 || (int)index == K-1))
          std::cout<<"set col: "<<i<<" " <<newdata.vals[(int)edge.time].QQ<<" " <<newdata.vals[(int)edge.time].QR<<" " <<std::endl;

        i++;
      }
    }

  }
     
  new_data = newdata;
  //for (int i=0; i<K; i++){ 
  //    assert(sumsum(newdata.vals[i].QQ)>0);
  //}
  counter[4] += t2.current_time();
}



size_t RMSE = 0;
static void merge_QQR(size_t index, const gl_dtypes::ishared_data& shared_data,
                      graphlab::any& current_data, const graphlab::any& new_data) {
  const QQR &newd = new_data.as<QQR>();
  QQR &curd = current_data.as<QQR>();
    
  curd.QQ += newd.QQ;
  curd.QR += newd.QR;
  //  std::cout << "*** Merge QQR index:" << index << " curd = " << curd.QR << std::endl;
}
static void merge_mult(size_t index, const gl_dtypes::ishared_data& shared_data,
                       graphlab::any& current_data, const graphlab::any& new_data) {
  const mult_QQR &newd = new_data.as<mult_QQR>();
  mult_QQR &curd = current_data.as<mult_QQR>();
  assert(newd.vals.size() == curd.vals.size() && newd.vals.size() == (unsigned int) K);  
 
  for (unsigned int i=0; i< newd.vals.size(); i++){ 
    curd.vals[i].QQ += newd.vals[i].QQ;
    curd.vals[i].QR += newd.vals[i].QR;
  }
     
  curd.rmse += newd.rmse;
     
  // HACK
  for (int i=0; i<20; i++){
    char  ss_[255];
    sprintf(ss_, "counter_%d", i);
    std::string ss(ss_);
    distributed_metrics::instance(__dc)->set_value(ss, counter[i]);
  }
}


void reduce_RMSE(size_t index, const gl_dtypes::ishared_data& shared_data,
                 gl_dtypes::iscope& scope, graphlab::any& new_data) {

  double curval = new_data.as<double>();
  //  double x =  double(curval + scope.const_vertex_data().rmse);
  new_data = double(curval + scope.const_vertex_data().rmse);
     
  //  std::cout << "*** Reduce RMSE index:" << index << " new = " << x << std::endl;
}

static void merge_RMSE(size_t index, const gl_dtypes::ishared_data& shared_data,
                       graphlab::any& current_data, const graphlab::any& new_data) {
  const double & newval = new_data.as<double>();
  double x =  double(newval + current_data.as<double>());
  current_data = double(newval + current_data.as<double>());
  std::cout << "*** merge RMSE index:" << index << " current = " << x<< std::endl;

} 

static void apply_MMT_U(size_t index, const gl_dtypes::ishared_data& shared_data,
                        graphlab::any& current_data, const graphlab::any& new_data) {
  printf("Apply %d: Apply MMT_U\n", (int) index);   
  assert(BPTF);
  timer t;
  t.start();
  const vec &Umean = new_data.as<QQR>().QR / M;
   
  const mat &UUT = new_data.as<QQR>().QQ;
  
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
  mat A_U = wishrnd(tmp, nu0_);
  mat tmp2;  
  ret =  inv(beta0_ * A_U, tmp2);
  assert(ret);
  vec mu_U = mvnrndex(mu0_, tmp2, D);
  if (debug)
    cout<<"Sampling from U" <<A_U<<" "<<mu_U<<" "<<Umean<<" "<<W0_<<tmp<<endl;

  QQR retm;
  retm.QQ = A_U;
  retm.QR = mu_U;
  current_data = retm;
  counter[5] += t.current_time();
}

static void apply_MMT_V(size_t index, const gl_dtypes::ishared_data& shared_data,
                        graphlab::any& current_data, const graphlab::any& new_data) {

  printf("Apply %d: Apply MMT_V\n", (int) index);     
  assert(BPTF);
  const vec &Vmean = new_data.as<QQR>().QR/N;
  const mat &VVT = new_data.as<QQR>().QQ;
 
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
  mat A_V = wishrnd(tmp, nu0_);
  mat tmp2; 
  ret = inv(beta0_*A_V, tmp2);
  assert(ret);
  vec mu_V = mvnrndex(mu0_, tmp2, D);
  if (debug)
    cout<<"Sampling from V: A_V" <<A_V<<" mu_V: "<<mu_V<<" Vmean: "<<Vmean<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;


  QQR retm;
  retm.QQ = A_V;
  retm.QR = mu_V;
  current_data = retm;



}


 

static void apply_func(size_t index,
                       const gl_dtypes::ishared_data& shared_data,
                       graphlab::any& current_data,
                       const graphlab::any& new_data) {

  double rmse = new_data.as<double>();
  rmse = sqrt(rmse/L);
  std::cout << myprocid << ":  Training RMSE is : " << rmse << std::endl;
  
  ASSERT_GT(rmse,0);
  
  // write the final result into the shared data table
  current_data = (double)rmse;
}


static void apply_QQR(size_t index,  const gl_dtypes::ishared_data& sdm,
                      graphlab::any& current_data, const graphlab::any& new_data) {

  printf("Apply %d: Apply QQR at proc %d\n", (int) index, myprocid);   
  if (debug && ((int)index == 0 || (int)index == K-1)){
    printf("entering time node  %d \n", (int)index);   
  } 

  mult_QQR data = new_data.as<mult_QQR>();
  mult_vec &ret = current_data.as<mult_vec>();
  //mult_vec othert = ret; 

  for (int i=0; i<K; i++){
    mat  QQ = data.vals[i].QQ;
    vec  QR = data.vals[i].QR;
 
    assert(i>=0 && i<K);
    if (debug && (i==0 || i == K-1))
      printf("node %d with Q size: %d\n", i, (int)D);



    if (debug && (i == 0 ||i== K-1 )){
      std::cout<<"QQ:"<<data.vals[i].QQ<<std::endl;
      std::cout<<"QR:"<<data.vals[i].QR<<std::endl;
    }

    assert(data.vals[i].QR.size() == D && data.vals[i].QQ.rows() == D && data.vals[i].QQ.cols() == D);
    vec sol;    
    vec out;

    QQR A_T = sdm.get(A_T_OFFSET).as<QQR>();
    assert(A_T.QQ.rows() == D && A_T.QR.size() == D);
    double alpha = sdm.get(ALPHA_OFFSET).as<double>();
    assert(alpha > 0.0001 && alpha < 1000);

    timer t;
    t.start();
    if (i == 0){
      vec t1 = ret.vals[1];
      if (!BPTF){
        QQ = QQ+ 2*eDT;
        bool ret = itpp::ls_solve(QQ, pT*(t1 + vones*muT)+ QR, out);
        assert(ret);
      }
      else {
        mat A_k = 2*A_T.QQ+alpha*QQ;
        mat iAk_;
        bool ret = inv(A_k, iAk_);
        assert(ret);
        vec muk_ = iAk_*(A_T.QQ*(t1 + A_T.QR)+alpha*QR); //TODO
        out = mvnrndex(muk_, iAk_, D);
      }
    }
    else if (i == K-1){
      vec tk_2 = ret.vals[K-2];
      if (!BPTF){
        QQ = QQ + eDT;
        bool ret = itpp::ls_solve(QQ, pT*tk_2 + QR, out);
        assert(ret); 
      }
      else {
        mat A_k = A_T.QQ+alpha*QQ;
        mat iAk_; 
        bool ret = inv(A_k, iAk_);
        assert(ret);
        vec muk_ = iAk_*(A_T.QQ*tk_2+alpha*QR); //TODO
        out = mvnrndex(muk_, iAk_, D);
      }
    }
    else {
      vec tsum = ret.vals[i-1] + ret.vals[i+1];
      if (!BPTF){
        QQ = QQ + 2*eDT;
        bool ret = itpp::ls_solve(QQ, pT*tsum + QR, out);
        assert(ret);
      }
      else {
        mat A_k = 2*A_T.QQ+alpha*QQ;
        mat iAk_;
        bool ret = inv(A_k, iAk_);
        assert(ret);
        vec muk_ = iAk_*(A_T.QQ* tsum +alpha*QR); //TODO
        out = mvnrndex(muk_, iAk_, D);
      }
    }

    //vec2vec(times[i].pvec, out, D);
    //sdm.atomic_set(i, out);

    if (debug && (i == 0|| i == K-1)) 
      std::cout<<out<<std::endl;
    
    assert(QQ.rows() == D && QQ.cols() == D);
    counter[6] += t.current_time();
    //}
    ret.vals[i] = out;

  }

  double rmse = data.rmse;
  rmse = sqrt(rmse/L);
  std::cout << "Training RMSE is : " << rmse << std::endl;
  // write the final result into the shared data table
  ret.rmse = rmse;
  
  distributed_metrics::instance(__dc)->set_value("residual", rmse);

}



void init_self_pot(graphlab::distributed_fullsweep_sdm<gl_dtypes::graph> &sdm){
  //assert(BPTF);

  W0 = eye(D);
  W0T = eye(D);
  iWalpha = 1.0/Walpha;
  iW0 = inv(W0);
  iW0T = inv(W0T);
  nu0 = D;


  printf("nuAlpha=%g, Walpha=%g, mu=%g, muT=%g, nu=%g, "
         "beta=%g, W=%g, WT=%g BURN_IN=%d\n", nuAlpha, Walpha, mu0, 
         mu0T, nu0, beta0, W0(1,1), W0T(1,1), BURN_IN);

  if (tensor || BPTF){
    QQR A_V, A_U, A_T;
    A_U.QQ = eye(D);
    A_V.QQ = eye(D);
    A_T.QQ = eye(D);

    A_U.QR = zeros(D); 
    A_V.QR = zeros(D); 
    A_T.QR = zeros(D);


    sdm.atomic_set(A_U_OFFSET, A_U);
    sdm.atomic_set(A_V_OFFSET, A_V);
    sdm.atomic_set(A_T_OFFSET, A_T);
  }
 
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
void load_pmf_distgraph(const char* filename, graph_dtype * g, bool test, distributed_control& dc);
//void calc_T(int id,  gl_types::iscope &scope, gl_types::icallback &scheduler, gl_types::ishared_data* shared_data);    
double calc_obj(gl_dtypes::ishared_data &sdm);
void last_iter(gl_dtypes::ishared_data &sdm);


void sample_alpha(double res2, gl_dtypes::ishared_data &sdm){
  double res;
  if (!tensor)
    res = powf(sdm.get(RMSE).as<double>(),2) * L;
  else res = powf(sdm.get(TIME_OFFSET).as<mult_vec>().rmse,2)*L;
  //assert(res > 0.1);

  printf("res vs. res2 %g %g\n", res, res2); 
  if (res < 1000)
    res = L * 3;

  // res = res2;
  assert(BPTF);
 
  std::cout << "sample alpha, nualpha = " << nuAlpha << std::endl;
 
 
  if (nuAlpha > 0){
    double nuAlpha_ =nuAlpha+ L;
    mat iWalpha_(1,1);
    iWalpha_.set(0,0, iWalpha + res);
    mat iiWalpha_ = zeros(1,1);
    iiWalpha_ = inv(iWalpha_);
    double alpha = wishrnd(iiWalpha_, nuAlpha_).get(0,0);
    assert(alpha != 0);

    if (debug)
      cout<<"Sampling from alpha" <<nuAlpha_<<" "<<iiWalpha_<<" "<<alpha<<endl;
    printf("sampled alpha is %g\n", alpha); 
    sdm.atomic_set(ALPHA_OFFSET, alpha);
  }
 

}


mat calc_MMT(int start_pos, int end_pos, vec &Umean){

  int batchSize = 1000;
  mat U = zeros(batchSize, D);
  mat MMT = zeros(D,D);
  int cnt = 0;
  timer t;

  for (int i=start_pos; i< end_pos; i++){
    if ((i-start_pos) % batchSize == 0){
      U.zeros();
      cnt = 1;
    }

    vertex_data * data= &g.vertex_data(i);
    vec* tmp = vec2vec(data->pvec, D);
    vec mean = *tmp;
    delete tmp;
    Umean += mean;
    //Q.set_col(k, ret);
    t.start(); 
    for (int s=0; s<D; s++)
      U(i%batchSize,s)=mean(s);
    if (debug && (i==start_pos || i == end_pos-1))
      std::cout<<" clmn "<<i<< " vec: " << mean <<std::endl;

    if ((cnt  == batchSize) || (cnt < batchSize && i == end_pos-1)){
      MMT = MMT+transpose(U)*U;
    }
    counter[8] += t.current_time();
    cnt++;
  }
  Umean /= (end_pos-start_pos);
  if (debug)
    cout<<"mean: "<<Umean<<endl;

  assert(MMT.rows() == D && MMT.cols() == D);
  assert(Umean.size() == D);
  return MMT;
}





mat calc_DT(gl_dtypes::ishared_data &sdm){

  assert(tensor);

  mat T = zeros(D, K);
  mult_vec tvec = sdm.get(TIME_OFFSET).as<mult_vec>();
  for (int i=0; i<K; i++){
    vec tmp = tvec.vals[i];
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

void sample_T(gl_dtypes::ishared_data& sdm){
  assert(BPTF);
  assert(tensor);

  double beta0_ = beta0 + 1;
  mult_vec tvec = sdm.get(TIME_OFFSET).as<mult_vec>();
  vec pvec = tvec.vals[0];
  vec mu0_ = (pvec + beta0*mu0T)/beta0_;
  double nu0_ = nu0 +K;
  //vec dMu = mu0 - Umean;
  if (debug){
    cout<<"beta0_ " << beta0_ << " beta0: " << beta0 << " nu0_ " << nu0_ << endl;
  } 

  mat dT = calc_DT(sdm);
  vec dTe = pvec - mu0T;
  mat iW0_ = iW0T + dT*transpose(dT) + (beta0/beta0_)*(itpp::outer_product(dTe,dTe));
  
  mat W0_;
  bool ret =inv(iW0_, W0_);
  assert(ret);
  mat tmp = W0_+transpose(W0_)*0.5;
  QQR A_T;
  A_T.QQ = wishrnd(tmp, nu0_);

  mat tmp2 ;
  ret = inv(beta0_*A_T.QQ, tmp2);
  assert(ret);
  A_T.QR = mvnrndex(mu0_, tmp2, D);
  if (debug)
    cout<<"Sampling from T: A_T" <<A_T.QQ<<" mu_V: "<<A_T.QR<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;
  
  sdm.atomic_set(A_T_OFFSET, A_T);
}


/*
  void T_update_function(gl_types::iscope &scope, gl_types::icallback &scheduler, gl_types::ishared_data* shared_data) {

  assert(tensor);

  int id = scope.vertex() - M-N;
  assert(id >=0 && id < K); 
  if (debug && (id == 0 || id == K-1)){
  printf("entering time node  %d \n", id);   
  } 
  if (K > 1)
  calc_T(id, scope, scheduler, shared_data); 
  }*/

double calc_rmse(graph_type * _g, bool test, double & res, gl_dtypes::ishared_data& sdm){

  if (test && Le == 0)
    return NAN;
   
  res = 0;
  double RMSE = 0;
  int e = 0;

  mult_vec tvec;
  if (tensor) 
    tvec = sdm.get(TIME_OFFSET).as<mult_vec>();

  for (int i=M; i< M+N; i++){ //TODO: optimize to start from N?
    vertex_data * data = &dg->vertex_data(i);
    foreach(edge_id_t iedgeid, _g->in_edge_ids(i)) {
         
      multiple_edges & edges = _g->edge_data(iedgeid);
      vertex_data * pdata = &dg->vertex_data(_g->source(iedgeid)); 
      for (int j=0; j< (int)edges.medges.size(); j++){       
 
        edge_data & edge = edges.medges[j];
        assert(edge.weight != 0);
        assert(edge.time < K);
        double sum = 0; 
        double add ;
        if (tensor) {
          double* tmp = vec2vec(&tvec.vals[(int)edge.time]);
          add = gl::rmse(data->pvec, pdata->pvec, tmp, D, edge.weight, sum);
          delete [] tmp;
        }
        else  {
          add = gl::rmse(data->pvec, pdata->pvec, NULL, D, edge.weight, sum);
        }
           
        assert(sum != 0);         
        //if (BPTF && iiter > BURN_IN)
        //   edge.avgprd += sum;

        if (debug && (i== M || i == M+N-1) && (e == 0 || e == (test?Le:L)))
          cout<<"RMSE:"<<i <<"u1"<< *vec2vec(data->pvec, D) << " v1 "<< *vec2vec(pdata->pvec, D)<<endl; 
        //assert(add<25 && add>= 0);
       
        //if (BPTF && iiter > BURN_IN){
        //  add = pow((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
        //}
         
        RMSE+= add;
        e++;
      }
    }
  }
  res = RMSE;
  assert(e == (test?Le:L));
  return sqrt(RMSE/(double)e);

}


void parse_edge(const edge_data& edge, const vertex_data & pdata, mat & Q, vec & vals, int i, const mult_vec & tvec){
        
  double buf[D];  
  double * pbuf = buf; 

  assert(edge.weight != 0);

  if (tensor){
    dot2(const_cast<double*>(pdata.pvec),  tvec.vals[(int)edge.time], buf, D);  
  }
  else {
    pbuf = const_cast<double*>(pdata.pvec);
  }
 
  for (int j=0; j<D; j++){
    Q.set(j,i, pbuf[j]); 
  }
  vals[i] = edge.weight;
}
 
void agg_QQ(const edge_data& edge, const vertex_data & pdata, const vertex_data & other, mat & QQ, vec & QR, int i){
        

  assert(tensor);
  assert(edge.weight != 0);

  vec tt = dot(pdata.pvec, other.pvec, D);  
  QQ += outer_product(tt,tt); 
  QR += tt*edge.weight;
  assert(QQ.rows() == D && QQ.cols() == D && QR.size() == D);
}
void agg_MM(const vertex_data & pdata, mat & QQ, vec & QR){
        

  assert(BPTF);
  vec* tmp = vec2vec((double*)pdata.pvec, D);
  vec tt = *tmp;
  delete tmp;
  QQ += outer_product(tt,tt); 
  QR += tt;
        
  assert(QQ.rows() == D && QQ.cols() == D && QR.size() == D);
}
 
int count_edges(graphlab::edge_list es){
  int cnt = 0; 
  for (int j=0; j< (int)es.size(); j++){
    cnt += dg->edge_data(es[j]).medges.size();
  }
  return cnt;
}


/***
 * UPDATE FUNCTION
 */
void pmf_update_function(gl_dtypes::iscope &scope, 
                         gl_dtypes::icallback &scheduler,
                         gl_dtypes::ishared_data* sdm) {
    

  bool debug = false;
  /* GET current vertex data */

  vertex_data& vdata = scope.vertex_data();
 
  
  /* CALCULATE new value */
  if (debug&& (scope.vertex() == 0 || (int)scope.vertex() == M-1 || (int)scope.vertex() == M || (int)scope.vertex() == M+N-1)){
    printf("entering %s node  %u \n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex());   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , vdata.pvec, D);
  }

  ASSERT_LE(scope.vertex(), M+N);

  vdata.rmse = 0;

 
  //int i=0; 
  // update neighbors 
  edge_list outs = scope.out_edge_ids();
  edge_list ins = scope.in_edge_ids();
  int numedges = vdata.num_edges;

  if (numedges == 0){
    return;
  }

  //assert(outs.size() == 0 || ins.size() == 0);   
  assert(numedges > 0);

  timer t;
  t.start(); 
  mat Q(D,numedges);
  vec vals(numedges);
  counter[0] += t.current_time();
  mult_vec tvec;
  if (tensor)
    tvec = sdm->get(TIME_OFFSET).as<mult_vec>();

  int i=0;

  if ((int)scope.vertex() < M){
    
    //MOVIE NODES
    foreach(graphlab::edge_id_t oedgeid, outs) {

      const multiple_edges &medges =scope.const_edge_data(oedgeid);
      const vertex_data & pdata = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
     
      for (int j=0; j< (int)medges.medges.size(); j++){
        const edge_data& edge = medges.medges[j];
        parse_edge(edge, pdata, Q, vals, i, tvec); 
     
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


    //USER NODES
    foreach(edge_id_t iedgeid, ins) {

      const vertex_data&  pdata = scope.const_neighbor_vertex_data(scope.source(iedgeid)); 
      multiple_edges & medges =scope.edge_data(iedgeid);
      for (int j=0; j< (int)medges.medges.size(); j++){
        edge_data& edge = medges.medges[j];
        parse_edge(edge, pdata, Q, vals, i, tvec); 
     
        if (debug && ((int)scope.vertex() == M || (int)scope.vertex() == M+N-1) && (i==0 || i == numedges-1))
          std::cout<<"set col: "<<i<<" " <<Q.get_col(i)<<" " <<std::endl;

        i++;
        //if (tensor) scheduler.add_task(gl_types::update_task(edge.time+M+N,
        //                         T_update_function),1);
        double sum;     
        double trmse;
              
        if (tensor) {
          double *tmp = vec2vec(&tvec.vals[(int)edge.time]);
          trmse = gl::rmse(vdata.pvec, const_cast<double*>(pdata.pvec), tmp, D, edge.weight, sum);
          delete [] tmp;
        }
        else {
          trmse = gl::rmse(vdata.pvec, const_cast<double*>(pdata.pvec), NULL, D, edge.weight, sum);
        }
              
        assert(sum != 0);
        if (BPTF && iiter > BURN_IN){
          edge.avgprd += sum;       
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

    QQR A = sdm->get((int)scope.vertex() < M ? A_U_OFFSET: A_V_OFFSET).as<QQR>();
    assert(Q.rows() == D && sumsum(A.QQ) > 0);
    t.start();
    mat iAi_;
    double alpha = sdm->get(ALPHA_OFFSET).as<double>();
    counter[16]+= t.current_time();

    assert(alpha > 0.00001 && alpha < 10000);
    bool ret =inv(A.QQ + alpha *  Q*itpp::transpose(Q), iAi_);
    counter[10]+= t.current_time();
    assert(ret);
    t.start();
    vec mui_ =  iAi_*(A.QQ*A.QR + alpha * Q * vals);
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

  if (scope.vertex() == 1)
    last_iter(*sdm);
  vdata.rounds++;
    
  if (sdm->get_constant(IS_ROUND_ROBIN_CONSTANT).as<bool>() == false) {
    if (vdata.rounds < (int) sdm->get_constant(MAX_ITERATIONS_CONSTANT).as<size_t>()) {
      gl_dtypes::update_task task(scope.vertex(), pmf_update_function);      
      scheduler.add_task(task, 100.0);
    }
  }
}

void last_iter(gl_dtypes::ishared_data &sdm){

  if (!tensor)
    sdm.trigger_sync(RMSE);
  if (tensor)
    sdm.trigger_sync(TIME_OFFSET);
  if (BPTF){
    sdm.trigger_sync(A_U_OFFSET);
    sdm.trigger_sync(A_V_OFFSET);
  }

  double res;
  double rmse = NAN;
  //double rmse = calc_rmse(&g, false, res, false);
  //rmse=0;
  //printf("%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f TEST RMSE=%0.4f.\n", gt.current_time(), BPTF?"BPTF":"ALS", iiter,calc_obj(sdm),  rmse, calc_rmse(&g1, true, res2, true, sdm));
  printf("%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f.\n", gt.current_time(), BPTF?"BPTF":"ALS", iiter,calc_obj(sdm),  rmse);
  iiter++;
  if (iiter == BURN_IN && BPTF){
    printf("Finished burn-in period. starting to aggregate samples\n");
  }
         
  if (BPTF){
    sample_alpha(res, sdm);
    //sample_U();
    //sample_V();
    if (tensor) 
      sample_T(sdm);
  }
}



double calc_obj(gl_dtypes::ishared_data &sdm){
   
  double res;  
  if (!tensor)
    res = powf(sdm.get(RMSE).as<double>(),2) * L;
  else res = powf(sdm.get(TIME_OFFSET).as<mult_vec>().rmse,2)*L;

  double sumU = 0, sumV = 0, sumT = 0;
  timer t;
  t.start(); 
  /*
    for (int i=0; i< M; i++){
    vertex_data * data = &g.vertex_data(i);
    sumU += gl::square_sum(data->pvec, D);
    } 

    for (int i=M; i< M+N; i++){
    vertex_data * data = &g.vertex_data(i);
    sumV += gl::square_sum(data->pvec, D);
    } */

  mat T;
  if (tensor){
    T = zeros(D,K);
    const mult_vec  tval = sdm.get(TIME_OFFSET).as<mult_vec>();
    for (int i=0; i<K; i++){
      sumT += pow(norm(tval.vals[i] - vones, 2),2);
      T.set_col(i, tval.vals[i]);
    }
    sumT *= pT;
  }
  counter[7]+= t.current_time();
  
  double obj = (res + pU*sumU + pV*sumV + sumT + (tensor?trace(T*dp*T.transpose()):0)) / 2.0;
  return obj;
}




gl_dtypes::iengine* tengine = NULL;
void barrier_fn(){ 
  std::cout << "Sync Barrier." << std::endl;
  // assert(false);
  /* ((graphlab::pushy_distributed_engine<gl_dtypes::graph,
     graphlab::distributed_scheduler_wrapper<gl_dtypes::graph, 
     graphlab::distributed_round_robin_scheduler<gl_dtypes::graph> >,
     general_scope_factory<gl_dtypes::graph> > *)(tengine))->sync_barrier();
  */
}


 
/** 
 * ==== SETUP AND START
 */
void start(int argc, char ** argv, distributed_control & dc) {
  ::command_line_options o;
  parse_command_line(argc, argv, o);
  __dc = &dc;
  myprocid = dc.procid();
  
 // IMPORTANT: all nodes need to seed random number generator similarly
  graphlab::random::seed(1);
 
  if (o.prepart) { 
    ASSERT_EQ(dc.numprocs(), 1);
    printf("loading data file %s\n", infile.c_str());
    load_pmf_graph(infile.c_str(), &g, false);

    printf("loading data file %s\n", (infile+"e").c_str());
    load_pmf_graph((infile+"e").c_str(),&g1, true);

    distributed_graph<vertex_data,multiple_edges>::partition_graph_tofile(g, o.ncpus, 
                                                                          graphlab::partition_method::PARTITION_RANDOM, infile);
    exit(0);
  }

  if (dc.procid() == 0){
    printf("loading test data file by node 0: %s\n", (infile+"e").c_str());
    load_pmf_graph((infile+"e").c_str(),&g1, true);
  }


  
  if (infile.find("netflix-r") != string::npos){
    M=480189; N=17770; K=27; L=99072112; Le=1408395;
  } else if (infile.find("netflix") != string::npos) {
    M=95526; N=3561; K=27; L=3298163; Le = 545177;
  }
  
  
  dg = new graph_dtype(dc);
  dg->set_constant_edges(true);

  load_pmf_distgraph(infile.c_str(), dg, false, dc);

  //dg.load(infile, dc);
  
  assert(dg->num_vertices() == (unsigned int) M+N);

  std::cout << M + N << " " << dg->num_vertices() << std::endl;

  dc.barrier();
  //g.distribute(dc);

  graphlab::distributed_fullsweep_sdm<gl_dtypes::graph> sdm(dc, o.ncpus, *dg);

  
  if (o.scheduler == "round_robin") {
    graphlab::pushy_distributed_engine<gl_dtypes::graph,
      graphlab::distributed_scheduler_wrapper<gl_dtypes::graph, 
      graphlab::distributed_round_robin_scheduler<gl_dtypes::graph> >,
      general_scope_factory<gl_dtypes::graph> > *d_engine = new  
      graphlab::pushy_distributed_engine<gl_dtypes::graph,
    graphlab::distributed_scheduler_wrapper<gl_dtypes::graph, 
    graphlab::distributed_round_robin_scheduler<gl_dtypes::graph> >,
    general_scope_factory<gl_dtypes::graph> >(dc, *dg, o.ncpus);
    // d_engine->set_caching(true);
    //d_engine->set_vertex_scope_pushed_updates(true);
    d_engine->set_default_scope(scope_range::EDGE_CONSISTENCY);
    // d_engine->set_max_backlog(1000);
    
    
    d_engine->set_shared_data_manager(&sdm);  
    tengine = d_engine;
  }
  else {
    assert(false);
  }
  //engine=tengine;
  if(tengine == NULL) {
    std::cout << "Unable to construct engine!" << std::endl;
    assert(false);
  }


  dc.barrier();
    

  //sdm.sync(PRIMAL_LOSS);   
  // compute QQ 
  // if (!tensor){ 

  //}
  if (tensor){ 
     
    mult_QQR mult;
    for (int i=0; i<K; i++)
      mult.vals.push_back(QQR());
    mult_vec  multret;
    for (int i=0; i<K; i++)
      multret.vals.push_back(ones(D)*0.1);
    //for (int i=0; i<K; i++){

    sdm.set_fullsweep_sync(TIME_OFFSET, sync_QQR, apply_QQR,merge_mult, mult, 100000000, 
                           scope_range::READ_CONSISTENCY,M, M+N);
    sdm.atomic_set(TIME_OFFSET, multret);
  }
  if (BPTF){
    sdm.set_fullsweep_sync(A_U_OFFSET, sync_MMT, apply_MMT_U,merge_QQR, QQR(), 400000000,
                           scope_range::VERTEX_READ_CONSISTENCY, 0, M);
    sdm.set_fullsweep_sync(A_V_OFFSET, sync_MMT, apply_MMT_V,merge_QQR, QQR(), 40000000,
                           scope_range::VERTEX_READ_CONSISTENCY, M, M+N);
  }

 
  if (tensor) 
    dp = GenDiffMat(K)*pT ;
 
  if (debug)
    std::cout<<dp<<std::endl;

  /**** CREATE INITIAL TASKS ******/
  
  /*std::vector<vertex_id_t> um;
    std::vector<vertex_id_t> tv;

    
    for (int i=0; i< M+N; i++)
    um.push_back(i);*/
  
  printf(" === %d: num of vertices: %d ====\n", 
         myprocid, int(dg->my_vertices().size()));
  
  // Check max and sum of edges
  size_t sumedges = 0, maxedges = 0;
  foreach(vertex_id_t v, dg->my_vertices()) {
    size_t n = dg->in_edge_ids(v).size() + dg->out_edge_ids(v).size();
    sumedges += n;
    maxedges = std::max(n, maxedges);
  }
  
  
  std::cout << myprocid << ": sum edges: " << sumedges << " max edges: " << maxedges << std::endl;
  
  
  //   tengine->get_scheduler().add_tasks(um, pmf_update_function, 1);
  if (o.scheduler == "round_robin") {
    vertex_id_t zero = 0;
    vertex_id_t Mv = M;

    tengine->get_scheduler().add_task_to_all(pmf_update_function, 1.0);

    tengine->get_scheduler().set_option(scheduler_options::MAX_ITERATIONS, (void*)NUM_ITERATIONS_TO_RUN);
    tengine->get_scheduler().set_option(scheduler_options::DISTRIBUTED_CONTROL, (void*)(&dc));
    tengine->get_scheduler().set_option(scheduler_options::BARRIER, &zero);
    tengine->get_scheduler().set_option(scheduler_options::BARRIER, &Mv);

    sdm.set_constant(IS_ROUND_ROBIN_CONSTANT, bool(true));
  }
  
  
 
  if (dc.procid() == 0){
    /* if (tensor){
       for (int i=M+N; i< M+N+K; i++)
       tv.push_back(i);
       engine->get_scheduler().add_tasks(tv, T_update_function, 1);
       }*/

    printf("%s for %s (%d, %d, %d):%d.  D=%d\n", BPTF?"BPTF":"PTF_ALS", tensor?"tensor":"matrix", M, N, K, L, D);
  
    if (!BPTF) {
      printf("pU=%g, pV=%g, pT=%g, muT=%g, D=%d\n", pU, pV, pT, muT,D);  
    }
  }


  dc.barrier();

  //if (BPTF)
  init_self_pot(sdm); //init anyway
  init_pmf();

  // double res, res2;
  //double rmse =  calc_rmse(&g, false, res, false);
  //printf("complete. Obj=%g, TEST RMSE=%0.4f.\n", calc_obj(sdm), calc_rmse(&g1, true, res2, false, sdm));
  printf("complete. Obj=%g\n", calc_obj(sdm));

  printf("BPTF: %d procid %d \n", (int) BPTF, (int) dc.procid());


  // Have to declare this from all procs
  distributed_metrics::instance(&dc)->set_value("residual", 0.0);
  distributed_metrics::instance(&dc)->set_value("custom_output_1", 0.0);

  if (dc.procid() == 0) {
    if (BPTF){
      sample_alpha(1*L, sdm);
      //sdm.sync( dg,A_U_OFFSET);
      //sdm.sync( dg,A_V_OFFSET);
      if (tensor) 
        sample_T(sdm);
    }
  }
  
  dc.barrier();
 

  if (BPTF){
    sdm.sync_from_local(A_U_OFFSET);
    sdm.sync_from_local(A_V_OFFSET);
  }
  if (!tensor) sdm.sync_from_local(RMSE);
  dc.barrier();
  /// Timing
  gt.start();



  /**** START GRAPHLAB *****/
  tengine->start();
  dc.barrier();
  if (dc.procid() == 0){
    double runtime = gt.current_time();
    double res2;
    double test_rmse = calc_rmse(&g1, true, res2, sdm);
    printf("Final result. Obj=%g, TEST RMSE= %0.4f.\n", calc_obj(sdm),  test_rmse);
    distributed_metrics::instance(&dc)->set_value("custom_output_1", test_rmse);

    /*
      printf("Final result. Obj=%g, TEST RMSE= %0.4f.\n", calc_obj(sdm),  calc_rmse(&g1, true, res2, sdm));
      sdm.sync(RMSE);
      if (!tensor)
      printf("SDM rmse is: %g\n", sqrt(sdm.get(RMSE).as<double>() / L));
      else  printf("SDM rmse is: %g\n", sqrt(sdm.get(TIME_OFFSET).as<mult_vec>().rmse / L));
    */
    /**** POST-PROCESSING *****/
    printf("Finished in %lf \n", runtime);
      
  
  } 
  dc.barrier();
  
  //timing counters
  for (int i=0; i<20; i++){
    printf("Counters are: %d/%d) %g\n", myprocid, i, counter[i]); 
    // char  ss_[255];
    // sprintf(ss_, "counter_%d", i);
    // std::string ss(ss_);
    // distributed_metrics::instance(&dc)->set_value(ss, counter[i]);

  }


}

template<typename edgedata>
int read_mult_edges(FILE * f, int nodes, graph_type * g, bool symmetry = false){
     

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

  for (int i=0; i<M; i++){
    //gl::rand(vdata.pvec, D);%TODO
    // g->add_vertex(vdata);
    g->add_vertex(vdata);
    if (debug && (i<= 5 || i == M-1))
      debug_print_vec("U: ", vdata.pvec, D);
  }
   
  for (int i=0; i<N; i++){
    //gl::rand(vdata.pvec, D);
    g->add_vertex(vdata);
    if (debug && (i<=5 || i==N-1))
      debug_print_vec("V: ", vdata.pvec, D);
  }
 
  /* 
     if (!test && tensor){
     //init times
     times = new vertex_data[K];
     vec tones = itpp::ones(D)*(K==1?1:0.1);
     for (int i=0; i<K; i++){
     vec2vec(times[i].pvec ,tones, D);
     g->add_vertex(times[i]);
     if (debug && (i <= 5 || i == K-1))
     debug_print_vec("T: ", times[i].pvec, D);
     }
     }
  */
  int val = read_mult_edges<edata2>(f, M+N, g);
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

//
// DISTRIBUTED GRAPH LOAD
//

template<typename edgedata>
int read_mult_edges_dist(FILE * f, int nodes, graph_dtype * g, bool symmetry = false){
     

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
             if (g->owner(ed[i].from-1) == myprocid ||g->owner(ed[i].to-1) == myprocid ) {
                g->edge_data(ret.second).medges.push_back(edge);
            }
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
  return e;
}




void load_pmf_distgraph(const char* filename, graph_dtype * g, bool test, distributed_control& dc) {
  int myprocid = dc.procid();
  int numprocs = dc.numprocs();
  
  printf("DISTRIBUTED (%d/%d) Loading %s %s\n", myprocid, numprocs, filename, test?"test":"train");
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

  for (int i=0; i<M; i++){
    //gl::rand(vdata.pvec, D);%TODO
    // g->add_vertex(vdata);
    g->add_vertex(graphlab::random::rand_int(numprocs-1), vdata);
    if (debug && (i<= 5 || i == M-1))
      debug_print_vec("U: ", vdata.pvec, D);
  }
   
  for (int i=0; i<N; i++){
    //gl::rand(vdata.pvec, D);
    g->add_vertex(graphlab::random::rand_int(numprocs-1), vdata);
    if (debug && (i<=5 || i==N-1))
      debug_print_vec("V: ", vdata.pvec, D);
  }
  
  int val = read_mult_edges_dist<edata2>(f, M+N, g);
  if (!test)
    L = val;
  else Le = val;

  if (!test && tensor && K>1) 
    edges = new std::vector<edge_id_t>[K]();

  // finalize
  g->finalize_dist(true);

  bool normalize = false;
  double normconst = 1;
  if (!strcmp(filename, "netflow") || !strcmp(filename, "netflowe")){
    normconst = 138088872;
    normalize = true;
  }
        

  //verify edges
  for (int i=M; i < M+N; i++){
    if (g->owner(i) == myprocid) {
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
  }
     
  if (!test){
    for (int i=0; i<M+N; i++){
      if (g->owner(i) == myprocid) {
          vertex_data &vdata = g->vertex_data(i);
          if (i < M)
            vdata.num_edges = count_edges(g->out_edge_ids(i));
          else
            vdata.num_edges = count_edges(g->in_edge_ids(i));
           
           g->update_vertex(i);
     
     }
    }
  }
  /*
  if (!test && tensor && K>1){
    int cnt = 0;
    for (int i=0; i<K; i++){
      cnt+= edges[i].size();
    }
    assert(cnt == L);
  } */
  fclose(f);
}



 
int main(int argc,  char *argv[]) {

  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::distributed_control dc(&argc, &argv);
  dc.init_message_processing(4);
  dc.barrier();
 

  double weight = -1;

  infile = argv[1];

  if (infile == "" || argc <= 3) {
    std::cout << "PMF <intput file> <weight> <options>\n";
    return EXIT_FAILURE;
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
  dc.barrier();
  start(argc, argv,dc);
   
}



#include <graphlab/macros_undef.hpp>
