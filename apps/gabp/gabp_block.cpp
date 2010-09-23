/*
 * \author akyrolaa, bickson
 */

#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include <itpp/base/smat.h>


//#include <cmath>
#include <cstdio>
#include <graphlab.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/app_support/iapp.hpp>
#include <graphlab/graph/graph.hpp>

typedef double sdouble; 

/*struct vdata{
  double mean;
  double real;
  double prec;
};*/


double counter[7];

#include <graphlab/util/timer.hpp>


// #include "../custom_schedulers/linear_scheduler.hpp"
#include "../linear_algebra/vecutils.hpp"
#include "../linear_algebra/blockutils.hpp"

#include <graphlab/macros_def.hpp>

using namespace graphlab;
using namespace itpp;


edge_id_t ** ees = NULL;
int TOTAL = 0;
sdouble THRESHOLD;

bool support_null_variance = false;
bool round_robin = false;
bool finish; //defined in convergence.hpp
int iiter = 0;//count number of time zero node run
int D = 1; //block size     

typedef Sparse_Mat<double> ssmat;
typedef Sparse_Vec<double> ssvec;
#define GABP_VEC_SIZE 4
#define GABP_MAT_SIZE 3

void pr(mat &pmat){
   std::cout<<pmat<<std::endl;
}
void pr(vec & pvec){
   std::cout<<pvec<<std::endl;
}


/** Vertex and edge data types **/
struct gabp_block_vertex_data {
  
  vec ** data;
   mat ** pmat;
   gabp_block_vertex_data(){
      data = new vec*[GABP_VEC_SIZE];
      for (int i=0; i<GABP_VEC_SIZE; i++){
         data[i] = new vec(D); data[i]->zeros();
      }
      pmat = new mat*[GABP_MAT_SIZE];
      for (int i=0; i< GABP_MAT_SIZE; i++){
         pmat[i] = new mat(D,D);
         pmat[i]->zeros();
      }
  }
}; 
 
#define GABP_PRIOR_MEAN_OFFSET 0
#define GABP_REAL_OFFSET 1
#define GABP_CUR_MEAN_OFFSET 2
#define GABP_PREV_MEAN_OFFSET  3 

#define GABP_PRIOR_PREC_OFFSET 0
#define GABP_CUR_PREC_OFFSET 1
#define GABP_PREV_PREC_OFFSET  2
struct gabp_block_edge_data {
  mat * weight; /// Example - can be changed
  vec * mean; /// Example - can be changed
  mat * prec;

   gabp_block_edge_data(){
     weight = new mat(D,D); weight->zeros();
     mean = new vec(D); //*mean = zeros(D);
     prec = new mat(D,D); prec->zeros();
  }
};
  

typedef graphlab::graph<gabp_block_vertex_data, gabp_block_edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;
gl_types::iengine * engine;
graph_type g;
gl_types::thread_shared_data sdm;


#define GABP_EDGE_WEIGHT_OFFSET 0

typedef std::pair<int,int> mypair;
typedef std::map<mypair, const gabp_block_edge_data*> mymap;


double CONV_THRESHOLD = 1e-3;
 
/*class convergence : public igraphlablistener {
        
  public:
    unsigned long output_priority_freq_tasks;
    unsigned long task_counter;
  
    //int n;
    //int m;
    int cur_mean_offset;
    int prev_mean_offset;
    int real_offset;
    int start_node_offset;
    int end_node_offset;
    double conv_threshold;
    double first_norm;
    bool first_time;
    bool error;

    convergence(unsigned long output_priority_freq_tasks = 1000) {
      // I.e how often to output max priority
      this->output_priority_freq_tasks = output_priority_freq_tasks;
      task_counter = 0;
      conv_threshold = CONV_THRESHOLD;
      first_time = true;
      error = false;
    }

    ~convergence(){
      //if (real != NULL) delete [] real;
    }
   
   void restart(){
     finish = false;
     first_time = true;
    
   }
  
    virtual void init(iengine* engine) {
       app_set_status("====== CONVERGENCE MONITOR STARTED ===== \n");
    }
    
    void scheduler_task_scheduled(update_task task, double current_max_priority) {
      if (finish)
          return;

      if (thread::GET_THREAD_ID() != 0)
          return;
 
      // Not thread safe!!!
      if (++task_counter % output_priority_freq_tasks == 0) {
                  
	  double real_norm = collect_diff(start_node_offset, end_node_offset, g, real_offset, cur_mean_offset,D);
          double relative_norm = collect_diff(start_node_offset, end_node_offset, g, prev_mean_offset, cur_mean_offset,D); 
          printf("task_counter %d relative error is %e change with last round %e \n", (int)task_counter, real_norm, relative_norm);

	 if ((relative_norm<conv_threshold && first_time) || (!first_time && (fabs(relative_norm/first_norm)<conv_threshold))){
               printf("Detected convergnce! aborting rest of queue\n");
               graphlab->get_scheduler().abort();
               finish = true;
         }
         if (relative_norm == NAN || real_norm == NAN || relative_norm == INFINITY || real_norm == INFINITY){
              printf("ALgorithm diverged! aborting run!\n");
              graphlab->get_scheduler().abort();
	      finish = true;
              error = true;
         }
         if (first_time){
             first_time = false;
             first_norm = relative_norm;
         }
      }

    }
    
    
    void set_status_str(char * str) {
       printf("Status %s\n", str); 
    }
};*/

/*void dispatch_vecb(int start_pos, int end_pos, int offset, graph * g, double val, int d = 1);
void dispatch_vecb(int start_pos, int end_pos, int offset, graph * g, double * _vec, bool free_vec = false, int d = 1);
void dispatch_mat(int start_pos, int end_pos, int offset, graph * g, double * _vec, bool free_vec = false, int d = 1);
double* collect_vec(int start_pos, int end_pos, int offset, graph * g, int d );
void read_nodesb(FILE * f, int len, int nodes, graph * g, int offset);
edge_id_t** read_block_edges(FILE * f, int len, int offset, int nodes, graph * g, bool symmetry = false, int d=1, ssmat * sm=NULL, int M =0, bool doparse = true, int num=0);
*/


/**
 * gabp APPLICATION CLASS
 */
class gabp_block_app : public iapp {
    
  public:
  int n,e,m;
  char ** argv;
  int argc;

/** CONSTRUCTOR **/
  gabp_block_app() {
    support_null_variance = false;
    THRESHOLD = 1e-10;
  }
    
  ~gabp_block_app() {
  }
    
   
  void start(); 
   
    
  void load_gabp_graph(const char* filename);    
    
}; // end class
  

/***
 * UPDATE FUNCTION
 */
void gabp_block_update_function(gl_types::iscope &scope, gl_types::icallback &scheduler,
                          gl_types::ishared_data * shared_data) {
    
  bool debug = false;
  timer t; 
  /* GET current vertex data */
  gabp_block_vertex_data * vertex_data = &scope.vertex_data();
 
  t.start(); 
  *vertex_data->data[GABP_PREV_MEAN_OFFSET] = *vertex_data->data[GABP_CUR_MEAN_OFFSET];
  *vertex_data->pmat[GABP_PREV_PREC_OFFSET] = *vertex_data->pmat[GABP_CUR_PREC_OFFSET];
  vec mu_i = *vertex_data->data[GABP_PRIOR_MEAN_OFFSET];
  mat J_i = *vertex_data->pmat[GABP_PRIOR_PREC_OFFSET];
  counter[0] += t.current_time();
   //if (!support_null_variance)
  //  assert(trace(J_i) != 0);

  /* CALCULATE new value */
  if (debug){
    printf("entering node  %d\n", (int)scope.vertex());   
    std::cout<< *vertex_data->pmat[GABP_PRIOR_PREC_OFFSET] << " " << *vertex_data->data[GABP_PRIOR_MEAN_OFFSET]<<std::endl;    
  }
  
  //std::vector<edge_id_t> inedgeid = scope.in_edge_ids();
  //assert(scope.in_edge_ids().size() == scope.out_edge_ids().size());
 
  //foreach(edge_id_t eid, inedgeid) {
 t.start();
  for (int T=0; T< TOTAL; T++){    
    if (ees[T][scope.vertex()] != 0xFFFFFFFF){

    const gabp_block_edge_data * edata = &scope.edge_data(ees[T][scope.vertex()]);
    if (debug)
     std::cout<<"edge from "<<scope.vertex() <<" prec: "<< T<< (*edata->prec) <<" weight: "<< (*edata->weight) <<"  mean: "<<*edata->mean<<std::endl;

    assert(edata != NULL);
    mu_i += (*edata->mean);
    J_i += *edata->prec;
    if (debug) std::cout<< (J_i) << " " << mu_i<< std::endl;   
  }
  } 
  counter[1] += t.current_time();
 
 
  if (debug) {
    printf("%d) summing up all messages\n", (int)scope.vertex());
    std::cout<< (J_i) << " " << mu_i<< std::endl;   
  }

  if (support_null_variance && trace((J_i)) == 0){
    *vertex_data->data[GABP_CUR_MEAN_OFFSET] = mu_i;
    *vertex_data->pmat[GABP_CUR_PREC_OFFSET] = (J_i);
  } else { 
    
   t.start(); 
   bool ret = itpp::ls_solve((J_i), mu_i, *vertex_data->data[GABP_CUR_MEAN_OFFSET]);
    counter[2]+= t.current_time();
    assert(ret);
    //assert(vertex_data->data[GABP_CUR_MEAN_OFFSET] != NAN);
    *vertex_data->pmat[GABP_CUR_PREC_OFFSET] = (J_i);
  }
  //assert(vertex_data->data[GABP_CUR_MEAN_OFFSET] != NAN);
  /* SEND new value and schedule neighbors */
  sdouble residual = norm(*vertex_data->data[GABP_PREV_MEAN_OFFSET]- *vertex_data->data[GABP_CUR_MEAN_OFFSET]) + 
       (norm(*vertex_data->pmat[GABP_PREV_PREC_OFFSET] - *vertex_data->pmat[GABP_CUR_PREC_OFFSET]));
  if (residual > THRESHOLD) {
    int i=0;
      
    for (int T=0; T<TOTAL; T++){
    //foreach(edge_id_t oedgeid, scope.out_edge_ids()) {
      // Get the reverse edge
      //graphlab::edge_id_t ineid = scope.reverse_edge(oedgeid);
      // Ensure that the edges are in the same order
      //gabp_block_edge_data * in_edge = (gabp_block_edge_data*)scope.edge_blob(inedgeid[i]).data;
      //gabp_block_edge_data * out_edge = (gabp_block_edge_data*)scope.mutable_edge_blob(oedgeid).data;
     if (ees[T][scope.vertex()]!= 0xFFFFFFFF){ 
      edge_id_t e = ees[T][scope.vertex()];
       gabp_block_edge_data * in_edge = &scope.edge_data(e);
      assert(ees[scope.vertex()][T] != 0xFFFFFFFF);
      gabp_block_edge_data * out_edge = &scope.edge_data(ees[scope.vertex()][T]);
       i++;   
        
      vec mu_i_j = mu_i - *in_edge->mean; 
      mat J_i_j = J_i - *in_edge->prec;

      //if (!support_null_variance)
      //  assert(trace(J_i_j) != 0);            

      assert(sumsum(abs((*out_edge->weight))) != 0);
       
     // if (support_null_variance && sumsum(J_i_j) == 0){
     //   *out_edge->mean = zeros(D);;
     //   *out_edge->prec = J_i_j;
     // } else {

       t.start(); 
       mat Jinv = -(out_edge->weight->transpose()) * inv(J_i_j);
       //mat Jinv = ls_solve_chol(J_i_j, -out_edge->weight->transpose());
        counter[3]+= t.current_time();
       t.start(); 
        
         *out_edge->mean = Jinv * mu_i_j;
        //bool ret = ls_solve(J_i_j, mu_i_j, *out_edge->mean);
        //assert(ret);
        //*out_edge->mean = -out_edge->weight->transpose() * *out_edge->mean;
         counter[4]+= t.current_time();
       
        t.start(); 
        *out_edge->prec = Jinv * *out_edge->weight;
        //ret = ls_solve_chol(J_i_j, *out_edge->weight, *out_edge->prec);
        //assert(ret);
        //*out_edge->prec = -out_edge->weight->transpose() * *out_edge->prec;
        counter[5]+= t.current_time();
     // }
      if (!finish && !round_robin)
        scheduler.add_task(gl_types::update_task(T, gabp_block_update_function),1);
      if (debug){
        printf("sending to %d  \n", T);
        std::cout<<"mean: " << *out_edge->mean<<" " << " prec: " << (*out_edge->prec) << " weight: " << (*out_edge->weight) <<std::endl;   
      }
      // break;

      //} else continue;          
    }
    }
     
     //shared_data->post(0, scope.vertex(), reducible_double(vertex_data->data[GABP_CUR_MEAN_OFFSET]-vertex_data->prev_mean));
    //} 
  }
  else {
    if (debug)
      printf("%d finished because residual is %e\n", (int)scope.vertex(), residual);
  }
 
  if (scope.vertex() == 0)
     iiter++;
  /* Write new value to vertex */
  //if (!support_null_variance)
  //  assert(vertex_data->prev_prec != 0);
}
 
/** 
 * ==== SETUP AND START
 */
void gabp_block_app::start() {
     
  bool debug = true; 
  /**** CREATE GRAPH ****/
  assert(argc>= 4);

    printf("setting block size to %d\n", atoi(argv[3]));
    D = atoi(argv[3]); 
      
  printf("loading data file %s\n", argv[1]);
  load_gabp_graph(argv[1]);
 
  printf("setting threshold %e\n", atof(argv[2]));
  THRESHOLD = atof(argv[2]); 
 
    engine_options o;
    o.parse_command_line(argc, argv);
    engine = o.create_engine(g);
    engine->set_shared_data_manager(&sdm);

  /*convergence * checker = new convergence(100);
  checker->g = g;
  checker->start_node_offset = 0;
  checker->end_node_offset = n+m;
  checker->real_offset = GABP_REAL_OFFSET;
  checker->prev_mean_offset = GABP_PREV_MEAN_OFFSET;
  checker->cur_mean_offset = GABP_CUR_MEAN_OFFSET;
    */  
  //graphlab = create_engine(g, checker);
  //checker->graphlab = graphlab;     
 
  /**** CREATE INITIAL TASKS ******/
  int k = g.num_vertices();
  engine->get_scheduler().add_task_to_all(gabp_block_update_function, 1);
      
      
  /// Timing
  timer t;
  t.start();
      
  /**** START GRAPHLAB *****/
  engine->start();
      
  /**** POST-PROCESSING *****/
  double runtime = t.current_time();
  printf("Finished in %lf\n", runtime);
      
  double diff(0);
  for (int i=0; i<k; i++){
    gabp_block_vertex_data * vdata = &g.vertex_data(i);
    diff += pow(norm(*vdata->data[GABP_REAL_OFFSET] - *vdata->data[GABP_CUR_MEAN_OFFSET]),2);
  }
  printf("gabp converged to an accuracy of %e\n", sqrt(diff));
  if (debug){
  	double * vec = collect_vec(0,n,GABP_CUR_MEAN_OFFSET,&g,D);
        debug_print_vec("result: ", vec, n*D);
        delete [] vec; 
  }
   
}


 
void gabp_block_app::load_gabp_graph(const char* filename) {

  bool debug = true;

  set_status("Loading %s", filename);
  FILE * f = fopen(filename, "r");
  assert(f!= NULL);

  fread(&n,1,4,f);
  fread(&m,1,4,f);
  //assert(m==0);
  printf("creating a model with %d nodes\n", n);
 
  n = ceil(n/(double)D);
  printf("divided into %d row blocks\n", D); 

  if (m != 0){
  	m = ceil(n/(double)D);
  	printf("X m nodes %d col \n", m*D); 
  }

  //g.set_finalized();
  
  read_nodesb(f, sizeof(gabp_block_vertex_data)/sizeof(sdouble), m==0?n:m,&g, GABP_PRIOR_MEAN_OFFSET,D);
  if (debug){
     double * ret = collect_vec(0,n,GABP_PRIOR_MEAN_OFFSET,&g,D);
     debug_print_vec("y:", ret, n*D);  
  }
   double * prec = read_vec(f, n*D);
  if (debug)
    debug_print_vec("prec:", prec, n*D);

  dispatch_mat(0,n,GABP_PRIOR_PREC_OFFSET, &g, prec, true,D);
  double * real = read_vec(f, n*D);
  if (debug)
    debug_print_vec("real: ", real, n*D);

  dispatch_vecb(0,n,GABP_REAL_OFFSET, &g, real, true,D);
  read_block_edges(f, sizeof(gabp_block_edge_data)/sizeof(sdouble), 0, D*(m+n),&g, true, D, NULL, 0, true, n);
  if (debug)
      debug_print_graph("graph", &g, 0, n, 0 , n, n,D);

  fclose(f);
}
  
  int main(int argc,  char *argv[]) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "Block gabp starting\n");

  assert(argc > 2);  
  gabp_block_app * app = new gabp_block_app();  
  app->parse_args(argc, argv);
  app->argv = argv; 
  app->argc = argc;
  app->start();
}

#include <graphlab/macros_undef.hpp>
 

