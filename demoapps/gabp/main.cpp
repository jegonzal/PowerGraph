/**
 * GRAPHLAB implementation of Gaussiabn Belief Propagation Code See
 * algrithm description and explanation in: Danny Bickson, Gaussian
 * Belief Propagation: Theory and Application. Ph.D. Thesis. The
 * Hebrew University of Jerusalem. Submitted October 2008.
 * http://arxiv.org/abs/0811.2518 By Danny Bickson, CMU. Send any bug
 * fixes/reports to bickson@cs.cmu.edu Code adapted to GraphLab by
 * Joey Gonzalez, CMU July 2010
 *
 * Functionality: The code solves the linear system Ax = b using
 * Gaussian Belief Propagation. (A is either square matrix or
 * skinny). A assumed to be full column rank.  Algorithm is described
 * in Algorithm 1, page 14 of the above Phd Thesis.
 *
 * If you are using this code, you should cite the above reference. Thanks!
 */

#include <cmath>
#include <cstdio>
#include "linear.h"
#include "gabp.hpp"
#include "jacobi.hpp"
#include "cg.hpp"
#include <graphlab/macros_def.hpp>


uint32_t n = 0; // number of rows of A
uint32_t m = 0; // number of cols of A (only used for non square matrix. In squre matrix the number is n)
uint32_t e = 0; // number of edges

bool square=true;

//read node data from file and add the nodes into the graph.
#define BUFSIZE 500000
template<typename graph>
void read_nodes(FILE * f, int len, int offset, int nodes,
                graph * g){

  assert(offset>=0 && offset < len);
  assert(nodes>0);

  int toread = nodes;
  if (nodes > BUFSIZE)
    toread = BUFSIZE;
  int remain = nodes;

  while(remain > 0){
    double * temp = new double[remain];
    toread = (remain < BUFSIZE)?remain:toread;
    int rc = (int)fread(temp, sizeof(double), toread, f);
    //assert(rc == toread);
    remain -= rc;
    vertex_data data;
    for (int i=0; i< rc; i++){
      if (offset == GABP_PRIOR_MEAN_OFFSET)
      	data.prior_mean = temp[i];
      else if (offset == GABP_REAL_OFFSET)
        data.real = temp[i];
      else assert(false);
      g->add_vertex(data);
    }
    delete [] temp;
  }
}

//read a vector from file and return an array
double * read_vec(FILE * f, size_t len){
  double * vec = new double[len];
  assert(vec != NULL);
  fread(vec, len, sizeof(double), f);
  return vec;
}
//write an output vector to file
void write_vec(FILE * f, int len, double * array){
  assert(f != NULL && array != NULL);
  fwrite(array, len, sizeof(double), f);
}


template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset,
                  graph * g, double * vec, int len,
                  bool free_vec = false, int d = 1){
  assert(end_pos - start_pos > 1);
  assert(len == end_pos - start_pos);
  assert(vec != NULL);
  assert(offset>=0 && offset < len);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++){
      data[offset+d*j] = vec[d*(i-start_pos)+j];
    }
  }
  if (free_vec)
    delete [] vec;
}

template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset,
                  graph * g,sdouble val,int d = 1){
  assert(end_pos - start_pos > 1);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++)
      data[offset+d*j] = val;
  }
}



//struct for holding edge data in file
struct edata{
  int from;
  int to;
  double weight;
};


//read edges from file into the graph
template<typename graph>
int read_edges(FILE * f, int len, int offset, int nodes,
               graph * g, bool symmetry = false){
  assert(offset>=0 && offset < len);

  typedef typename graph::edge_data_type etype;

  unsigned int e,g0;
  int rc = fread(&e,1,4,f);
  assert(rc == 4);
  rc = fread(&g0,1,4,f); //zero pad
  assert(rc == 4);
  assert(g0 == 0);

  printf("Creating %d edges...\n", e);
  assert(e>0);
  int total = 0;
  edata* ed = new edata[200000];
  printf("symmetry: %d\n", symmetry);
  int edgecount_in_file = e;
  if (symmetry) edgecount_in_file /= 2;
  while(true){
    memset(ed, 0, 200000*sizeof(edata));
    rc = (int)fread(ed, sizeof(edata),
                    std::min(200000, edgecount_in_file - total), f);
    total += rc;

    edge_data tmp;
    for (int i=0; i<rc; i++){
      //memset(tmp, 0, len/sizeof(sdouble));
      tmp.weight =  ed[i].weight;
      assert(ed[i].weight != 0);
      assert(ed[i].from >= 1 && ed[i].from <= nodes);
      assert(ed[i].to >= 1 && ed[i].to <= nodes);
      assert(ed[i].to != ed[i].from);
      // Matlab export has ids starting from 1, ours start from 0
      g->add_edge(ed[i].from-1, ed[i].to-1, tmp);
      if (!square) { //add the reverse edge as well
        // Matlab export has ids starting from 1, ours start from 0
        g->add_edge(ed[i].to-1, ed[i].from-1, tmp);
      }
    }
    printf(".");
    fflush(0);
    if (rc == 0 || total >= edgecount_in_file)
      break;
  }
  delete [] ed; ed = NULL;
  return e;
}


//helper function to compute the norm between the true solution and the current computed mean
double get_real_norm(const vertex_data& v){
  return pow(v.real - v.cur_mean, 2);
}

//helper function to compute the norm between last round iteration and this round of the current mean
double get_relative_norm(const vertex_data& v){
  return pow(v.cur_mean - v.prev_mean, 2);
}


/**
 * Engine terminates when this returns true
 */
bool termination_condition() {
  double ret = RELATIVE_NORM_KEY.get_val();

  //std::cout<<"I was in term"<<std::endl;
  if (ret < THRESHOLD_KEY.get()){
    std::cout << "Aborting since relative norm is: " << ret << std::endl;
    return true;
  }
  return false;
}


/*
 * aggregate norm of real answer and the current solution
 */
static void apply_func_real(graphlab::any& current_data,
                            const graphlab::any& new_data) {

  double ret = new_data.as<double>();
  ret = sqrt(ret);
  //std::cout << "Real Norm is: " << ret << std::endl;
  // write the final result into the shared data table
  current_data = ret;
}

/*
 * aggregate norm of difference in messages
 */

static void apply_func_relative(graphlab::any& current_data,
                                const graphlab::any& new_data) {

  double ret = new_data.as<double>();
  ret = sqrt(ret);
  std::cout << "Relative Norm is: " << ret << std::endl;
  // write the final result into the shared data table
  current_data = (double)ret;

}


void load_square_matrix(const char* filename, graph_type& graph) {

  printf("Loading %s\n", filename);
  FILE * f = fopen(filename, "r");
  assert(f!= NULL);

  fread(&n, 1, 4, f);
  fread(&m, 1, 4, f);
  assert(m==0);
  read_nodes(f, sizeof(vertex_data)/sizeof(sdouble), GABP_PRIOR_MEAN_OFFSET,
             n, &graph);

  double * real = read_vec(f, n);
  dispatch_vec(0,n,GABP_REAL_OFFSET, &graph, real, n, true);

  double * prec = read_vec(f, n);
  dispatch_vec(0,n,GABP_PRIOR_PREC_OFFSET, &graph, prec, n, true);

  e = read_edges(f, sizeof(edge_data)/sizeof(sdouble), 0, n, &graph);
  fclose(f);
}

void load_non_square_matrix(const char* filename, graph_type& graph) {
  printf("Loading %s\n", filename);
  FILE * f = fopen(filename, "r");
  assert(f!= NULL);

  fread(&m,1,4,f);
  fread(&n,1,4,f);
  assert( n > 0);
  assert( m > 0);
  assert(m!=n); 
  printf("Loading a non-square matrix A of size %d x %d\n", m,n);
  read_nodes(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_PRIOR_MEAN_OFFSET,m,&graph);
  read_nodes(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_REAL_OFFSET,n,&graph);

  double * prec = read_vec(f, n+m);
  dispatch_vec(0,n+m,GABP_PRIOR_PREC_OFFSET, &graph, prec, n+m, true);
  dispatch_vec(0,n+m,GABP_PREV_MEAN_OFFSET, &graph, 1);
  dispatch_vec(0,n+m,GABP_PREV_PREC_OFFSET, &graph, 1);
  e = read_edges(f, sizeof(edge_data), 0, n+m, &graph);
  fclose(f);
}



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");
  // Setup additional command line arguments for the GABP program
  std::string datafile;
  double threshold = 1e-5;
  bool support_null_variance = false;
  bool finish = false;
  bool debug = false;
  size_t iter = 0;
  int syncinterval = 10000;
  int algorithm;

  clopts.attach_option("algorithm", &algorithm, "Algorithm 0=Gaussian BP, 1= Jacobi");
  clopts.add_positional("algorithm");
 
  clopts.attach_option("data", &datafile, "Binary input file (as created by the save_c_gl.m script)");
  clopts.add_positional("data");
  clopts.attach_option("threshold", &threshold, threshold, "termination threshold.");
  clopts.add_positional("threshold");
  clopts.attach_option("nullvar", &support_null_variance, support_null_variance,
                       "(optional) support invalid covariance matrices - with null variance.");
  clopts.attach_option("square", &square, square, "is the matrix square? ");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("syncinterval", &syncinterval, syncinterval, "sync interval (number of update functions before convergen detection");
  
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  // Ensure that a data file is provided
  if(!clopts.is_set("data")) {
    logstream(LOG_ERROR)<<"No data file provided!" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
 // Ensure that a data file is provided
  if(!clopts.is_set("algorithm")) {
    logstream(LOG_ERROR)<< "No algorithm provided! Choose from: 0) GaBP 1) Jacobi" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  // Collect some additional information from the command line options
  bool round_robin = false;
  if (clopts.get_scheduler_type() == "round_robin") round_robin = true;

  // Create a core
  gl_types::core core;
  core.set_engine_options(clopts); // Set the engine options


  // Load the graph --------------------------------------------------
  if (!square){
     load_non_square_matrix(datafile.c_str(), core.graph());
     if (algorithm == JACOBI){
        logstream(LOG_ERROR)<<" Jacobi can not run with non-square mastrix. Run with --sqaure=true and provide a square mastrix in the input file!" << std::endl;
        return EXIT_FAILURE;
     }
                        
  }
  else //square matrix
        load_square_matrix(datafile.c_str(), core.graph());



  // Initialize the shared data --------------------------------------
  // Set syncs
  //
  if (syncinterval > 0){
    core.set_sync(REAL_NORM_KEY,
                gl_types::glshared_sync_ops::sum<double, get_real_norm>,
                apply_func_real,
                double(0),  syncinterval,
                gl_types::glshared_merge_ops::sum<double>);
  
    core.set_sync(RELATIVE_NORM_KEY,
                gl_types::glshared_sync_ops::sum<double, get_relative_norm>,
                apply_func_relative,
                double(0),  syncinterval,
                gl_types::glshared_merge_ops::sum<double>);
  }
  // Create an atomic entry to track iterations (as necessary)
  ITERATION_KEY.set(0);
  // Set all cosntants
  THRESHOLD_KEY.set(threshold);
  SUPPORT_NULL_VARIANCE_KEY.set(support_null_variance);
  ROUND_ROBIN_KEY.set(round_robin);
  FINISH_KEY.set(finish);
  DEBUG_KEY.set(debug);
  MAX_ITER_KEY.set(iter);

  core.graph().compute_coloring();

  // CREATE INITIAL TASKS ******
  // TOOD is this the correct starting priority?
  double initial_priority = 1.0;
  switch(algorithm){
    case GaBP: 
	  core.add_task_to_all(gabp_update_function, initial_priority); break;

    case JACOBI:
          core.add_task_to_all(jacobi_update_function, initial_priority); break;
    
    default:
         logstream(LOG_ERROR) << "Unknown algorithm" << std::endl;
         clopts.print_description(); 
         return EXIT_FAILURE;
  } 

  // Add the termination condition to the engine
  if (syncinterval > 0)
  	core.engine().add_terminator(termination_condition);

  // START GRAPHLAB *****
  double runtime = core.start();

  // POST-PROCESSING *****
  std::cout << algorithmnames[algorithm] << " finished in " << runtime << std::endl;

  std::vector<double> means(m+n);
  std::vector<double> precs(m+n);
  double diff = 0;
  for (size_t i = m; i < core.graph().num_vertices(); i++){
    const vertex_data& vdata = core.graph().vertex_data(i);
    diff += ((vdata.real - vdata.cur_mean)*
             (vdata.real - vdata.cur_mean));
     means[i] = vdata.cur_mean;
     precs[i] = vdata.cur_prec;
  }
  std::cout << "Assuming the linear system is Ax=y, and the correct solution is x*," << algorithmnames[algorithm] << " converged to an accuracy norm(x-x*) of " << diff
            << " msg norm is: " << RELATIVE_NORM_KEY.get_val() << std::endl;

   FILE * f = fopen((datafile+".out").c_str(), "w");
   assert(f!= NULL);

   std::cout<<"Writing result to file: "<<datafile<<".out"<<std::endl;
   std::cout<<"You can read the file in Matlab using the load_c_gl.m matlab script"<<std::endl;
   write_vec(f, means.size(), &means[0]);
   write_vec(f, means.size(), &precs[0]);

   fclose(f);
   return EXIT_SUCCESS;
}

#include <graphlab/macros_undef.hpp>

