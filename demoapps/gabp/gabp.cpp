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
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>


//define sdouble as either float or double as needed
typedef double sdouble;


enum constant_offsets {GABP_PRIOR_MEAN_OFFSET = 0, //prior mean (b_i / A_ii)
                       GABP_PRIOR_PREC_OFFSET = 1, //prior precision P_ii = A_ii
                       GABP_REAL_OFFSET = 2, // the real solution (if known) x = inv(A)*b for square matrix or x = pinv(A) for skinny matrix
                       GABP_CUR_MEAN_OFFSET = 3, //intermediate value to store the mean \mu_i
                       GABP_CUR_PREC_OFFSET = 4, //current precision value P_i
                       GABP_PREV_MEAN_OFFSET = 5,// mean value from previous round (for convergence detection)
                       GABP_PREV_PREC_OFFSET = 6}; // precision value from previous round (for convergence detection)




/** Vertex and edge data types **/
struct vertex_data {
  sdouble prior_mean;  //prior mean (b_i / A_ii)
  sdouble prior_prec;   //prior precision P_ii = A_ii
  sdouble real;   //the real solution (if known) x = inv(A)*b for square matrix or x = pinv(A) for skinny matrix
  sdouble cur_mean; //intermediate value to store the mean \mu_i
  sdouble cur_prec; // //current precision value P_i
  sdouble prev_mean; //  mean value from previous round (for convergence detection)
  sdouble prev_prec; //precision value from previous round (for convergence detection)

  vertex_data():prev_mean(1000000){ };

};

//edge is a scalar  non zero entry A_{ij} in the matrix A (row i, col j)
struct edge_data {
  sdouble weight; //edge value
  sdouble mean; // message \mu_ij
  sdouble prec; // message P_ij
};


  uint32_t n = 0; // number of rows of A
  uint32_t m = 0; // number of cols of A (only used for non square matrix. In squre matrix the number is n)
  uint32_t e = 0; // number of edges


typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;


gl_types::glshared<double> REAL_NORM_KEY;
gl_types::glshared<double> RELATIVE_NORM_KEY;
gl_types::glshared<size_t> ITERATION_KEY;
gl_types::glshared_const<double> THRESHOLD_KEY;
gl_types::glshared_const<bool> SUPPORT_NULL_VARIANCE_KEY;
gl_types::glshared_const<bool> ROUND_ROBIN_KEY;
gl_types::glshared_const<bool> FINISH_KEY;
gl_types::glshared_const<bool> DEBUG_KEY;
gl_types::glshared_const<size_t> MAX_ITER_KEY;
//read node data from file and add the nodes into the graph.
#define BUFSIZE 500000
template<typename graph>
void read_nodes(FILE * f, int len, int offset, int nodes,
                graph * g){
  typedef typename graph::vertex_data_type vtype;

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
      sdouble * pdata = (sdouble*)&data;
      pdata[offset] = temp[i];
      g->add_vertex(*(vtype*)&data);
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
                  graph * g,float * vec, int len,
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

    sdouble tmp[len/sizeof(sdouble)];
    for (int i=0; i<rc; i++){
      memset(tmp, 0, len/sizeof(sdouble));
      tmp[offset] =  ed[i].weight;
      assert(ed[i].weight != 0);
      assert(ed[i].from >= 1 && ed[i].from <= nodes);
      assert(ed[i].to >= 1 && ed[i].to <= nodes);
      assert(ed[i].to != ed[i].from);
      // Matlab export has ids starting from 1, ours start from 0
      g->add_edge(ed[i].from-1, ed[i].to-1, *(etype*)&tmp);
      if (symmetry) { //add the reverse edge as well
        // Matlab export has ids starting from 1, ours start from 0
        g->add_edge(ed[i].to-1, ed[i].from-1, *(etype*)&tmp);
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



/***
 * UPDATE FUNCTION
 * \todo briefly describe what this function is doing?
 */
void gabp_update_function(gl_types::iscope &scope,
                          gl_types::icallback &scheduler) {


  /* GET current vertex data */
  vertex_data& vdata = scope.vertex_data();
  graphlab::edge_list inedgeid = scope.in_edge_ids();
  graphlab::edge_list outedgeid = scope.out_edge_ids();

  const double& threshold = THRESHOLD_KEY.get();
  const bool& support_null_variance  = SUPPORT_NULL_VARIANCE_KEY.get();
  const bool& round_robin = ROUND_ROBIN_KEY.get();
  const bool& finish = FINISH_KEY.get();
  const bool& debug = DEBUG_KEY.get();



  //store last round values
  vdata.prev_mean = vdata.cur_mean;
  vdata.prev_prec = vdata.cur_prec;

  //initialize accumlated values
  sdouble mu_i = vdata.prior_mean;
  sdouble J_i = vdata.prior_prec;
  if (!support_null_variance) assert(J_i != 0);

  /* CALCULATE new value */
  if (debug) {
    std::cout << "entering node " << scope.vertex()
              << " P=" << vdata.prior_prec
              << " u=" << vdata.prior_mean
              << std::endl;
  }

  //accumlate all messages (the inner summation in section 4 of Algorithm 1)
  foreach(gl_types::edge_id_t eid, inedgeid) {
    const edge_data& edata = scope.edge_data(eid);
    mu_i += edata.mean;
    J_i +=  edata.prec;
  }

  if (debug) {
    std::cout << scope.vertex() << ") summing up all messages "
              << mu_i << " " << J_i << std::endl;
  }

  // optional support for null variances
  if (support_null_variance && J_i == 0){
    vdata.cur_mean = mu_i;
    vdata.cur_prec = 0;
  } else {
    assert(J_i != 0);
    vdata.cur_mean = mu_i / J_i;
    assert(vdata.cur_mean != NAN);
    vdata.cur_prec = J_i;
  }
  assert(vdata.cur_mean != NAN);

  /* SEND new value and schedule neighbors */
  sdouble residual =
    fabs(vdata.prev_mean- vdata.cur_mean) +
    fabs(vdata.prev_prec - vdata.cur_prec);

    for(size_t i = 0; i < inedgeid.size(); ++i) {
      assert(scope.source(inedgeid[i]) == scope.target(outedgeid[i]));
      edge_data& in_edge = scope.edge_data(inedgeid[i]);
      edge_data& out_edge = scope.edge_data(outedgeid[i]);
      graphlab::vertex_id_t target = scope.target(outedgeid[i]);

      //substruct the sum of message sent from node j
      sdouble mu_i_j = mu_i - in_edge.mean;
      sdouble J_i_j  = J_i - in_edge.prec;

      if (!support_null_variance)  assert(J_i_j != 0);
      assert(out_edge.weight != 0);

      if (support_null_variance && J_i_j == 0){
        out_edge.mean = 0;
        out_edge.prec = 0;
      } else {
        //compute the update rule (Section 4, Algorithm 1)
        out_edge.mean = -(out_edge.weight * mu_i_j / J_i_j);
        out_edge.prec = -((out_edge.weight * out_edge.weight) / J_i_j);//matrix is assumed symmetric!
      }

      if (!finish && !round_robin) {
        gl_types::update_task task(target, gabp_update_function);
        double priority = fabs(vdata.cur_prec) + 1e-5;
        scheduler.add_task(task, priority);
      }

      if (debug) {
        std::cout << "Sending to " << target << " "
                  << out_edge.mean << " "
                  << out_edge.prec << " wdge weight "
                  << out_edge.weight << std::endl;
      }
    }


}






void load_gabp_graph(const char* filename, graph_type& graph) {

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

void load_gabp_graph2(const char* filename, graph_type& graph) {
  printf("Loading %s\n", filename);
  FILE * f = fopen(filename, "r");
  assert(f!= NULL);

  fread(&m,1,4,f);
  fread(&n,1,4,f);
  assert( n > 0);
  assert( m > 0);
  printf("Loading a graph of size %d x %d\n", n,m);
  read_nodes(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_PRIOR_MEAN_OFFSET,m,&graph);
  read_nodes(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_REAL_OFFSET,n,&graph);

  double * prec = read_vec(f, n);
  dispatch_vec(0,n,GABP_PRIOR_PREC_OFFSET, &graph, prec, n, true);
  dispatch_vec(0,n+m,GABP_PREV_MEAN_OFFSET, &graph, 1);
  dispatch_vec(0,n+m,GABP_PREV_PREC_OFFSET, &graph, 1);
  e = read_edges(f, sizeof(edge_data), 0, n+m, &graph);
  fclose(f);
}



int main(int argc,  char *argv[]) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);


  graphlab::command_line_options clopts("Gaussian Belief Propagation");

  // Setup additional command line arguments for the GABP program
  std::string datafile;
  double threshold = 1e-5;
  bool support_null_variance = false;
  bool finish = false;
  bool debug = false;
  bool square = true;
  size_t iter = 0;
  int syncinterval = 10000;

  clopts.attach_option("data", &datafile, "Binary input file (as created by the save_c_gl.m script)");
  clopts.add_positional("data");
  clopts.add_positional("threshold");
  clopts.attach_option("threshold", &threshold, threshold, "termination threshold.");
  clopts.attach_option("nullvar",
                       &support_null_variance, support_null_variance,
                       "(optional) support invalid covariance matrices - with null variance.");
  //clopts.attach_option("finish", &finish, finish, "?finish?.");
  clopts.attach_option("square",
                       &square, square,
                       "is the matrix square? ");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("iter", &iter, iter, "maximum allowed iterations (optional).");
  clopts.attach_option("syncinterval", &syncinterval, syncinterval, "sync interval (number of update functions before convergen detection");
  //clopts.scheduler_type = "round_robin";

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  // Ensure that a data file is provided
  if(!clopts.is_set("data")) {
    std::cout << "No data file provided!" << std::endl;
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
  if (!square)
         load_gabp_graph2(datafile.c_str(), core.graph());
  else //square matrix
        load_gabp_graph(datafile.c_str(), core.graph());



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

  /**** CREATE INITIAL TASKS ******/
  // TOOD is this the correct starting priority?
  double initial_priority = 1.0;
  core.add_task_to_all(gabp_update_function, initial_priority);

  // Add the termination condition to the engine
  if (syncinterval > 0)
  	core.engine().add_terminator(termination_condition);

  /**** START GRAPHLAB *****/
  double runtime = core.start();

  /**** POST-PROCESSING *****/
  std::cout << "Finished in " << runtime << std::endl;

  double means[m+n];
  double precs[m+n];
  double diff = 0;
  for (size_t i = 0; i < core.graph().num_vertices(); i++){
    const vertex_data& vdata = core.graph().vertex_data(i);
    diff += ((vdata.real - vdata.cur_mean)*
             (vdata.real - vdata.cur_mean));
     means[i] = vdata.cur_mean;
     precs[i] = vdata.cur_prec;
  }
  std::cout << "gabp converged to an accuracy of "
            << diff << " after " << REAL_NORM_KEY.get_val() << " relative norm: " << RELATIVE_NORM_KEY.get_val() << std::endl;

   FILE * f = fopen((datafile+".out").c_str(), "w");
   assert(f!= NULL);

   std::cout<<"Writing result to file: "<<datafile<<".out"<<std::endl;
   std::cout<<"You can read the file in Matlab using the load_c_gl.m matlab script"<<std::endl;
   write_vec(f, m+n, means);
   write_vec(f, m+n, precs);

   fclose(f);

  return EXIT_SUCCESS;
}

#include <graphlab/macros_undef.hpp>

