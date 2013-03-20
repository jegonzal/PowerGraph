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


/**
 * Functionality: The code solves the linear system Ax = b using
 * The Jacobi algorithm. (A is a square matrix). 
 * A assumed to be full column rank.  Algorithm is described
 * http://en.wikipedia.org/wiki/Jacobi_method
 * Written by Danny Bickson 
 */
#include "../collaborative_filtering/eigen_wrapper.hpp"
#include "../collaborative_filtering/types.hpp"
#include "../collaborative_filtering/eigen_serialization.hpp"
#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>



enum jacobi_fields{
  JACOBI_X = 0,
  JACOBI_REAL_X = 1,
  JACOBI_Y = 2,
  JACOBI_PREV_X = 3,
  JACOBI_PREC = 4
};

int actual_vector_len = 5;
int data_size = 5;
bool final_residual = true;
bool zero = false;  //allow for zero entries in sparse matrix market format
bool update_function = false;
double ortho_repeats = 3;
std::string vecfile;
int rows = 0, cols = 0;
int max_iter = 10;
double tol = 1e-5;
int quiet = 0;
int unittest = 0;

struct vertex_data {
  vec pvec;
  double A_ii;
  //real_type y, Aii;
  //real_type pvec[JACOBI_X], pvec[JACOBI_REAL_X], pvec[JACOBI_PREV_X];
  vertex_data(): A_ii(1) { //: y(0), Aii(1), pvec[JACOBI_X](0), pvec[JACOBI_REAL_X](0), 
                 // pvec[JACOBI_PREV_X](-1) 
     pvec = zeros(data_size);
     pvec[JACOBI_PREV_X] = -1;
  }
  void save(graphlab::oarchive& arc) const { 
    arc << pvec << A_ii;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) { 
    arc >> pvec >> A_ii;
  }
}; // end of vertex_data

class gather_type {
  public:
    vec pvec;
    double training_rmse;
    double validation_rmse;
    gather_type() { training_rmse= validation_rmse = 0; }
    void save(graphlab::oarchive& arc) const { arc << pvec << training_rmse << validation_rmse; }
    void load(graphlab::iarchive& arc) { arc >> pvec >> training_rmse >> validation_rmse; }  
    gather_type& operator+=(const gather_type& other) {
      pvec += other.pvec;
      training_rmse += other.training_rmse;
      validation_rmse += other.validation_rmse;
      return *this;
    } 

};

//gather_type ret;



struct edge_data : public graphlab::IS_POD_TYPE {
  double obs;
  int role;
  enum data_role_type { TRAIN, VALIDATE, PREDICT  };

  edge_data(double obs = 1, data_role_type role = TRAIN) :
    obs(obs), role(role) { }

}; // end of edge data


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
graph_type * pgraph;

/**
 * \brief Given a vertex and an edge return the other vertex in the
 * edge.
 */
inline graph_type::vertex_type
get_other_vertex(graph_type::edge_type& edge, 
    const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex

//typedef double gather_type;
typedef double message_type;

void start_engine();
#include "../collaborative_filtering/math.hpp"


void verify_values(int unittest, double residual){
   if (unittest == 1)
     assert(residual < 1e-5);
}
/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph, 
    const std::string& filename,
    const std::string& line) {

  //no need to parse
  if (boost::algorithm::ends_with(filename ,vecfile))
    return true;

  ASSERT_FALSE(line.empty()); 
  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  
  // Parse the line
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  strm >> source_id >> target_id;
  source_id--; target_id--;
  if (source_id >= (uint)rows)
    logstream(LOG_FATAL)<<"Row number: " << source_id << " sould be < rows " << rows << " [ line: " << line << " ] " << std::endl;
  if (target_id >= (uint)cols)
    logstream(LOG_FATAL)<<"Col number: " << target_id << " sould be < cols " << cols << " [ line: " << line << " ] " << std::endl;
  strm >> obs;
  if (!info.is_square())
  target_id = rows + target_id;

  if (source_id == target_id){
      vertex_data data;
      data.A_ii = obs;
      data.pvec[JACOBI_PREC] = obs;
      graph.add_vertex(source_id, data);
  }
  // Create an edge and add it to the graph
  else graph.add_edge(source_id, target_id, edge_data(obs, role)); 
  return true; // successful load
} // end of graph_loader


#include "../collaborative_filtering/math.hpp" //uses vertex_data and edge_data so has to be included here
#include "../collaborative_filtering/printouts.hpp" // the same
typedef graphlab::omni_engine<Axb> engine_type;
engine_type * pengine = NULL;

struct linear_model_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;

  int pos;
  linear_model_saver(int pos): pos(pos) {}

  std::string save_vertex(const vertex_type& vertex) const {
     assert(pos >= 0 && pos < vertex.data().pvec.size());
     std::string ret;
     ret = boost::lexical_cast<std::string>(vertex.id() + 1) + " ";
     ret += boost::lexical_cast<std::string>(vertex.data().pvec[pos]) + "\n";
     return ret;
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
}; 


void start_engine(){
  vertex_set nodes = pgraph->select(selected_node);
  pengine->signal_vset(nodes);
  pengine->start();
}

int main(int argc, char** argv) {
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Solve a linear system using Jacobi method";
  graphlab::command_line_options clopts(description);
  std::string input_dir, output_dir;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
      "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("initial_vector", vecfile,"optional initial vector");
  clopts.attach_option("debug", debug, "Display debug output.");
  clopts.attach_option("unittest", unittest,  
      "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("max_iter", max_iter, "max iterations");
  clopts.attach_option("regularization", regularization, "regularization");
  clopts.attach_option("tol", tol, "convergence threshold");
  clopts.attach_option("rows", rows, "number of rows");
  clopts.attach_option("cols", cols, "number of cols");
  clopts.attach_option("quiet", quiet, "quiet mode (less verbose)");
  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  if (quiet){
    global_logger().set_log_level(LOG_ERROR);
    debug = false;
  }

  if (rows <= 0 || cols <= 0 || rows != cols)
    logstream(LOG_FATAL)<<"Please specify number of rows/cols of the input matrix" << std::endl;
    
 
  info.rows = rows;
  info.cols = cols;

  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; 
  graph_type graph(dc, clopts);  
  graph.load(input_dir, graph_loader); 
  pgraph = &graph;
  dc.cout() << "Loading graph. Finished in " 
    << timer.current_time() << std::endl;
  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " 
    << timer.current_time() << std::endl;


  dc.cout() 
    << "========== Graph statistics on proc " << dc.procid() 
    << " ==============="
    << "\n Num vertices: " << graph.num_vertices()
    << "\n Num edges: " << graph.num_edges()
    << "\n Num replica: " << graph.num_replicas()
    << "\n Replica to vertex ratio: " 
    << float(graph.num_replicas())/graph.num_vertices()
    << "\n --------------------------------------------" 
    << "\n Num local own vertices: " << graph.num_local_own_vertices()
    << "\n Num local vertices: " << graph.num_local_vertices()
    << "\n Replica to own ratio: " 
    << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
    << "\n Num local edges: " << graph.num_local_edges()
    //<< "\n Begin edge id: " << graph.global_eid(0)
    << "\n Edge balance ratio: " 
    << float(graph.num_local_edges())/graph.num_edges()
    << std::endl;

  dc.cout() << "Creating engine" << std::endl;
  engine_type engine(dc, graph, exec_type, clopts);
  pengine = &engine;

  init_math(&graph, info, ortho_repeats, update_function);

  if (vecfile.size() > 0){
    std::cout << "Load b vector from file" << input_dir << vecfile << std::endl;
    FILE * file = fopen((input_dir + vecfile).c_str(), "r");
    if (file == NULL)
      logstream(LOG_FATAL)<<"Failed to open initial vector"<< std::endl;
    vec input = vec::Zero(rows);
    double val = 0;
    for (int i=0; i< rows; i++){
      int rc = fscanf(file, "%lg\n", &val);
      if (rc != 1)
        logstream(LOG_FATAL)<<"Failed to open initial vector"<< std::endl;
      input[i] = val;
    }
    fclose(file);
    DistVec v0(info, JACOBI_Y, false, "v0");
    v0 = input;
  }  

  dc.cout() << "Running Jacobi" << std::endl;
  dc.cout() << "(C) Code by Danny Bickson, CMU " << std::endl;
  dc.cout() << "Please send bug reports to danny.bickson@gmail.com" << std::endl;
  timer.start();

  DistMat A(info);
  DistVec b(info, JACOBI_Y,true, "b");
  DistVec x(info, JACOBI_X,true, "x", JACOBI_PREV_X);
  DistVec A_ii(info, JACOBI_PREC, true, "A_ii");

  PRINT_VEC(b);
  PRINT_VEC(x);
  PRINT_VEC(A_ii);
  for (int i=0; i < max_iter; i++){
    mi.use_diag = false;
    x = (b - A*x)/A_ii;
    PRINT_VEC(x);
  }
 
  dc.cout() << "Jacobi finished in " << runtime << std::endl;
  dc.cout() << "\t Updates: " << engine.num_updates() << std::endl;

    DistVec p(info, JACOBI_PREV_X, true, "p");
    mi.use_diag = true;
    p = A*x -b;
    PRINT_VEC(p);
    DistDouble ret = norm(p);
    dc.cout() << "Solution converged to residual: " << ret.toDouble() << std::endl;
 
  //vec ret = fill_output(&core.graph(), info, JACOBI_X);
  //write_output_vector(datafile + "x.out", format, ret, false);
  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
    << std::endl
    << "Final Runtime (seconds):   " << runtime 
                                        << std::endl
                                        << "Updates executed: " << engine.num_updates() << std::endl
                                        << "Update Rate (updates/second): " 
                                          << engine.num_updates() / runtime << std::endl;

  graph.save("x.out", linear_model_saver(JACOBI_X), false, true, false, 1);
  graphlab::mpi_tools::finalize();

   return EXIT_SUCCESS;
}










