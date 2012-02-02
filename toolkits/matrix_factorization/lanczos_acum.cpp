/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

#include <graphlab.hpp>
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include <graphlab/graph/graph2.hpp>
#include <graphlab/graph/graph3.hpp>
#include <graphlab/graph/graph.hpp>
#include <google/malloc_extension.h>
#include <graphlab/macros_def.hpp>
using namespace graphlab;
using namespace std;


/**
 *
 *  Implementation of the Lanczos algorithm, as given in:
 *  http://en.wikipedia.org/wiki/Lanczos_algorithm
 * 
 *  Code written by Danny Bickson, CMU, June 2011
 * */



//LANCZOS VARIABLES
vec lancbeta;
vec lancalpha;
bipartite_graph_descriptor info;
int max_iter = 10;
bool debug;
bool fix_init;
bool measure_save_time = false;
bool edge_data_flag = true;
bool is_square = false;

struct vertex_data {
  vec pvec;
  double value;
  vertex_data(){
     pvec.resize(max_iter+1);
  }
  void add_self_edge(double value) { }

  void set_val(double value, int field_type) { 
    pvec[field_type] = value;
  }
void save(graphlab::oarchive &oarc) const {
  }

  void load(graphlab::iarchive &iarc) {
  }

  //double get_output(int field_type){ return pred_x; }
}; // end of vertex_data

struct edge_data {
  real_type weight;
  edge_data(double weight = 1) : weight(weight) { }
void save(graphlab::oarchive &oarc) const {
    oarc << weight;
  }


  void load(graphlab::iarchive &iarc) {
    iarc >> weight;
  }

 };

#define USE_GRAPH_VER 3

#if USE_GRAPH_VER == 1
typedef graphlab::graph<vertex_data, edge_data> graph_type;
#elif USE_GRAPH_VER == 2
typedef graphlab::graph2<vertex_data, edge_data> graph_type;
#elif USE_GRAPH_VER == 3
typedef graphlab::graph3<vertex_data, edge_data> graph_type;
#endif


/**
 *
 * [n,k] = size(A);
   V = zeros(k,m+1);
   V(:,2) = rand(k,1);
   V(:,2)=V(:,2)/norm(V(:,2),2);
   beta(2)=0;
 *
 * */
void init_lanczos(graph_type * g, bipartite_graph_descriptor & info){

   int m = max_iter;
   assert(m > 0);
   lancbeta = zeros(m+3);
   lancalpha = zeros(m+3);
   double sum = 0;

  for (int i = info.get_start_node(false); i< info.get_end_node(false); i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec = zeros(m+3);
    data->pvec[1] = fix_init ? 0.5 : rand();
    sum += data->pvec[1]*data->pvec[1];
  }

  sum = sqrt(sum);
  for (int i=info.get_start_node(false); i< info.get_end_node(false); i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec[1] /= sum;
    if (debug && i- info.get_start_node(false) < 20)
      std::cout<<"Initial V(:,2) is " << data->pvec[1] << std::endl;
  }
  for (int i = info.get_start_node(true); i< info.get_end_node(true); i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec = zeros(m+3);
  }
}
void print_v(bool rows, int offset, graph_type * g){

  int start= info.get_start_node(rows);;
  int end = info.get_end_node(rows);;
  vec v = zeros(end-start);
  for (int i=start; i< end; i++){ 
    vertex_data * data = &g->vertex_data(i);
    v[i - start] = data->pvec[offset];
  }
  cout<<"v is: " << mid(v,0,20) << endl;
  if (end - start > 40)
    cout<<"v end is: " << mid(v, v.size()-20, 20) << endl;
}

/***
 * UPDATE FUNCTION (ROWS)
 */

struct lanczos_update :
   public graphlab::iupdate_functor<graph_type, lanczos_update>{
   void operator()(icontext_type &context){


  /* GET current vertex data */
  vertex_data& user = context.vertex_data();
  int m = context.get_global<int>("m");
  int id = context.vertex_id(); 
  timer t; t.start(); 
  user.value = 0;

  if (info.is_row_node(id)){
 
  int offset = context.get_global<int>("offset");
 
  /* print statistics */
  if (debug && info.toprint(id)){
    printf("Lanczos ROW Axb: entering  node  %d \n",  id);   
 //   debug_print_vec("V" , user.pvec, m);
  }

  const edge_list_type outs= context.out_edges();
  if (outs.size() == 0)
     return;

  // for (size_t i=0; i< outs.size(); i++) {
      // const edge_type& e = outs[i];
  foreach (const edge_type& e, outs) {
      const edge_data& edge = context.const_edge_data(e);
      const vertex_data  & movie = context.const_vertex_data(e.target());
      user.value += edge.weight * movie.pvec[offset];
  }

  if (debug && info.toprint(id)){
    printf("Lanczos ROWS computed value  %d %g \n",  id, user.value);   
  }

  } else 
  { //column node

  
  if (debug && info.toprint(id)){
    printf("Lanczos COLS: entering  node  %d \n",  id);   
  //  debug_print_vec("V" , user.pvec, m);
  }
  
  int offset2 = context.get_global<int>("offset2");
  int offset3 = context.get_global<int>("offset3");

  const edge_list_type ins = context.in_edges();
  if (ins.size() == 0)
    return;

   // for (size_t i=0; i< ins.size(); i++) {
      // const edge_type& e = ins[i];
   foreach (const edge_type& e, ins) {
      const edge_data& edge = context.const_edge_data(e);
      const vertex_data  & movie = context.const_vertex_data(e.source());
       user.value += edge.weight * movie.value;
   }
   
   assert(offset2 < m+2 && offset3 < m+2);
   user.value -= lancbeta[offset2] * user.pvec[offset3];

  if (debug && info.toprint(id)){
    printf("Lanczos COLS computed value  id:%d val:%g beta: %g v %g (%g)\n",  id, 
        user.value,lancbeta[offset2],  user.pvec[offset3], lancbeta[offset2]*user.pvec[offset3]);   
  }

  }
}

  void operator+=(const lanczos_update& other) { 
  }

  void finalize(iglobal_context_type& context) {
  } 
};
  
double wTV(int j, graph_type *g){

  double lancalpha = 0;
  for (int i=info.get_start_node(false); i< info.get_end_node(false); i++){ 
    const vertex_data * data = &g->vertex_data(i);
    lancalpha+= data->value*data->pvec[j];
  }
  if (debug)
	cout<<"alpha: " << lancalpha<<endl;

  return lancalpha;
}

void substruct(int curoffset, int j, graph_type* g, double alpha){
  assert(j >= 0 && j < curoffset);
  assert(alpha != 0);

  for (int i=info.get_start_node(false); i< info.get_end_node(false); i++){ 
    vertex_data * data = &g->vertex_data(i);
    data->value -= alpha * data->pvec[j];
  }
}


void orthogolonize_vs_all(int curoffset, graph_type *g){
  for (int i=1; i< curoffset-1; i++){
     double alpha = wTV(i, g);
     if (alpha != 0)
       substruct(curoffset, i, g, alpha);
  }
}

double w_norm_2(bool rows, graph_type * g){
  double norm = 0;
  for (int i=info.get_start_node(rows); i< info.get_end_node(rows); i++){ 
    vertex_data * data = &g->vertex_data(i);
    norm += data->value*data->value;
  }
  return sqrt(norm);
}



void w_minus_lancalphaV(int j, graph_type * g){
  
  if (debug)
	cout << "w: " ;
  for (int i=info.get_start_node(false); i< info.get_end_node(false); i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->value -= lancalpha[j]*data->pvec[j];
  }
}

void update_V(int j, graph_type * g){

  if (debug)
	cout<<"V: ";

  for (int i=info.get_start_node(false); i< info.get_end_node(false); i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->value /= lancbeta[j]; //DB: NEW
    data->pvec[j] = data->value;
    if (debug && i-info.get_start_node(false)<20)
	cout<<data->pvec[j]<<" ";
  }

  if (debug)
	cout<<endl;
}

mat calc_V(graph_type * g){

  mat V = zeros(info.num_nodes(false), max_iter+1);
  for (int i=info.get_start_node(false); i< info.get_end_node(false); i++){ 
    const vertex_data * data = (vertex_data*)&g->vertex_data(i);
    set_row(V, i-info.get_start_node(false), mid(data->pvec, 1, max_iter+1));
  }
  return V;
}


void print_w(bool rows, graph_type * g){

  int size = info.num_nodes(rows);
  int start= info.get_start_node(rows);
  int end = info.get_end_node(rows);
  vec v = zeros(size);
  for (int i=start; i< end; i++){ 
    const vertex_data * data = (vertex_data*)&g->vertex_data(i);
    v[i - start] = data->value;
  }
  cout<<"w is: " << mid(v,0,20) << endl;
  if (end - start > 40)
    cout<<"w end is: " << mid(v, v.size()-20, 20) << endl;
}


void compute_residual(const vec & eigenvalues, const mat & eigenvectors, graph_type * g, graphlab::core<graph_type, lanczos_update> & glcore){


  if (debug)
   cout<<"Eigenvectors (in resid: ) " << eigenvectors << endl;
 
  for (int j=1; j< max_iter; j++){
      glcore.set_global("offset", j);
      glcore.set_global("offset3", j);
      glcore.set_global("offset2",j);
      lancbeta[j] = eigenvalues[j-1];
      if (debug)
          printf("Residual eigenvalue %d is: %g\n", j, eigenvalues[j-1]);
      for (int i= info.get_start_node(false); i< info.get_end_node(false); i++){
         g->vertex_data(i).pvec[j] = get_val( eigenvectors, i-info.get_start_node(false), eigenvectors.cols() - j);
         //printf("%g ", g->vertex_data(i).pvec[j]);
      }  
      //printf("\n");
      //glcore.schedule_all(lanczos_update());
      //glcore.start();
      cout<<"Residual for eigenvalue " << j << endl;
      print_v(false, j ,&glcore.graph());
     glcore.aggregate_now("sync");
      print_w(true, &glcore.graph());
      print_w(false, &glcore.graph());
       double sum = 0;
      for (int i= info.get_start_node(false); i< info.get_end_node(false); i++){
        sum += pow(g->vertex_data(i).value,2);

      }
      printf("Residual for eigenvalue %d %g is: %g\n", j, eigenvalues[j-1], sqrt(sum));

  }


}


void lanczos(graphlab::core<graph_type, lanczos_update> & glcore, 
             bipartite_graph_descriptor & info, timer & mytimer){
   
   glcore.set_global("m", max_iter);

   init_lanczos(&glcore.graph(), info);

   lanczos_update lupdate;
   glcore.add_aggregator("sync", lupdate, 1000);
   
   //for j=2:m+2
   for (int j=1; j<= max_iter+1; j++){
        //w = A*V(:,j) 
        glcore.set_global("offset", j);
        glcore.set_global("offset3", j-1);
        glcore.set_global("offset2",j);
        glcore.aggregate_now("sync");

        if (debug){
          print_w(true,&glcore.graph());
          print_w(false,&glcore.graph());
        }
       
        //lancalpha(j) = w'*V(:,j);
        //w =  w - lancalpha(j)*V(:,j);
      	lancalpha[j] = wTV(j, &glcore.graph());
        w_minus_lancalphaV(j, &glcore.graph());
        orthogolonize_vs_all(j+1,&glcore.graph());
    
        if (debug)
          print_w(false,&glcore.graph());

        //lancbeta(j+1)=norm(w,2)a
        lancbeta[j+1] = w_norm_2(false, &glcore.graph());

        //V(:,j+1) = w/lancbeta(j+1);
        update_V(j+1, &glcore.graph()); 
        logstream(LOG_INFO) << "Finished iteration " << j << " in time: " << mytimer.current_time() << std::endl;

   } 
  /* 
 * T=sparse(m+1,m+1);
 * for i=2:m+1
 *     T(i-1,i-1)=lancalpha(i);
 *     T(i-1,i)=lancbeta(i+1);
 *     T(i,i-1)=lancbeta(i+1);
 * end 
 * T(m+1,m+1)=lancalpha(m+2);
 * V = V(:,2:end-1);
 */
 int m = max_iter;
 mat T=zeros(m+1,m+1);
 for (int i=1; i<=m; i++){
   set_val(T,i-1,i-1,lancalpha[i]);
   set_val(T,i-1,i,lancbeta[i+1]);
   set_val(T,i,i-1,lancbeta[i+1]);
 }
 set_val(T,m,m,lancalpha[m+1]);
 if (debug && m < 100){
    cout<<"Matrix T is: " << T << endl;
 }

 mat Vectors=calc_V(&glcore.graph()); 
   
 vec eigenvalues; 
 mat eigenvectors;
 assert(::eig_sym(T, eigenvalues, eigenvectors));
 cout << "Here are the computed eigenvalues, from larger to smaller" << endl;
 for (int i=0; i< std::min((int)eigenvalues.size(),20); i++)
	cout<<"eigenvalue " << i << " val: " << sqrt(eigenvalues[i]) << endl;



 mat U=Vectors*eigenvectors;
 //compute_residual(eigenvalues, U, &glcore.graph(), glcore);
 if (debug)
   cout<<"Eigen vectors are:" << U << endl << "V is: " << Vectors << endl << " Eigenvectors (u) are: " << eigenvectors;
 mat V=zeros(eigenvalues.size(),1);
 set_col(V,0,eigenvalues); 

}

int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile, outdir;
  std::string format = "matrixmarket";
  size_t sync_interval = 10000;
  int unittest = 0;
  int num_rows = 0;

  bool load_bin = false;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix input file");
  clopts.add_positional("data");
  clopts.attach_option("outdir", &outdir, outdir, "Output directory");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("syncinterval", 
                       &sync_interval, sync_interval, 
                       "sync interval (number of update functions before convergen detection");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("max_iter", &max_iter, max_iter, "max iterations");
  clopts.attach_option("fix_init", &fix_init, fix_init, "fix random vector init to be const"); 
  clopts.attach_option("num_rows", &num_rows, num_rows, "number of matrix rows");
  clopts.attach_option("measure_save_time", &measure_save_time, measure_save_time, "Measure save time and exit");
  clopts.attach_option("load_bin", &load_bin, load_bin, "Load binary format");
  clopts.attach_option("edge_data_flag", &edge_data_flag, edge_data_flag, "Allow edge data");
  clopts.attach_option("is_square", &is_square, is_square, "square matrix?");
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab V2 matrix factorization library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Lanczos" << std::endl;


  // Create a core
  graphlab::core<graph_type, lanczos_update> core;
#if USE_GRAPH_VER == 2
  core.graph().set_use_vcolor(false);
#endif
  core.set_options(clopts); // Set the engine options
  core.set_scope_type("vertex");
  //unit testing
  if (unittest == 1){
    datafile = "lanczos2";
  }

  logstream(LOG_WARNING)<<"Using Graph Version: " << USE_GRAPH_VER << std::endl;

  std::cout << "Load matrix " << datafile << std::endl;
  timer mytimer; mytimer.start(); 
#if USE_GRAPH_VER != 3
  if (load_bin) {
    core.graph().load(datafile+".out");
  } else {
    load_graph(datafile, format, info, core.graph());
  }
  core.graph().finalize();
#if USE_GRAPH_VER == 2
  if (measure_save_time){
    logstream(LOG_INFO) << "Going to save " << mytimer.current_time() << endl;
    core.graph().save(datafile + ".out");
    logstream(LOG_INFO) << "Saved in time " << mytimer.current_time() << endl;
    graph_type graph0;
    logstream(LOG_INFO) << "Going to load " << mytimer.current_time() << endl;
    graph0.load(datafile + ".out");
    logstream(LOG_INFO) << "Load in time " << mytimer.current_time() << endl;
    exit(1);
  }
#endif
  //save_to_bin("/usr0/bickson/" + datafile, core.graph(), true);
   // save_to_bin("/usr1/haijieg/netflix3/" + datafile, core.graph(), true);
   // exit(1);
#else
  core.graph().load_directed(outdir+ datafile, false, !edge_data_flag);
#endif
  info.nonzeros = core.graph().num_edges();
  if (is_square)
     info.cols = num_rows;
  else 
  info.cols = core.graph().num_vertices() - num_rows;
  info.rows = num_rows;

  logstream(LOG_INFO)<<"Loadded a matrix of size " << info.rows << "x" << info.cols << "  nnz: " << info.nonzeros << endl;
  logstream(LOG_INFO)<<"Time taken to load graph: " << mytimer.current_time() << std::endl; 

   // After loading. Google TMalloc Profile
   //
   size_t value;
   MallocExtension::instance()->GetNumericProperty("generic.heap_size", &value);
   std::cout << "Heap Size: " << (double)value/(1024*1024) << "MB" << "\n";
   MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &value);
   std::cout << "Allocated Size: " << (double)value/(1024*1024) << "MB" << "\n";


  std::cout << "Schedule all vertices" << std::endl;

  core.schedule_all(lanczos_update());

  if (sync_interval < core.graph().num_vertices()){
    sync_interval = core.graph().num_vertices(); 
    logstream(LOG_WARNING) << "Sync interval is lower than the number of nodes: setting sync interval to " 
			   << sync_interval << std::endl;
  }
 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
  core.add_global("offset", int(0));
  core.add_global("offset2", int(0));
  core.add_global("offset3", int(0));
  core.add_global("m", int(0));


  lanczos(core, info, mytimer);
  std::cout << "Lanczos finished in " << mytimer.current_time() << std::endl;
  std::cout << "\t Updates: " << core.last_update_count() << " per node: " 
     << core.last_update_count() / core.graph().num_vertices() << std::endl;

  //vec ret = fill_output(&core.graph(), bipartite_graph_descriptor, JACOBI_X);

  //write_output_vector(datafile + "x.out", format, ret);


  if (unittest == 1){
  }

   return EXIT_SUCCESS;
}

#include <graphlab/macros_undef.hpp>

