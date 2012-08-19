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
 *
 *  Implementation of the Lanczos algorithm, as given in:
 *  http://en.wikipedia.org/wiki/Lanczos_algorithm
 * 
 *  Code written by Danny Bickson, CMU, June 2011
 * */


#include "eigen_wrapper.hpp"
#include "types.hpp"
#include "eigen_serialization.hpp"
#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>



//when using negative node id range, we are not allowed to use
//0 and 1 so we add 2.
const static int SAFE_NEG_OFFSET=2;
int iter = 0;
//LANCZOS VARIABLES
int max_iter = 10;
bool no_edge_data = false;
int actual_vector_len;
int nv = 0;
int nsv = 0;
double tol = 1e-8;
bool finished = false;
double ortho_repeats = 3;
bool update_function = false;
bool save_vectors = false;
std::string datafile; 
std::string vecfile;
int unittest;
std::string format = "matrixmarket";
int nodes = 0;
int data_size = 0;

void start_engine();

struct vertex_data {
  /** \brief The number of times this vertex has been updated. */
  uint32_t nupdates;
  /** \brief The most recent L1 change in the pvec value */
  float residual; //! how much the latent value has changed
  /** \brief The latent pvec for this vertex */
  vec pvec;
  double A_ii;

  /** 
   * \brief Simple default constructor which randomizes the vertex
   *  data 
   */
  vertex_data() : nupdates(0), residual(1), A_ii(0) { randomize(); } 
  /** \brief Randomizes the latent pvec */
  void randomize() { pvec.resize(data_size); pvec.setRandom(); }
  /** \brief Save the vertex data to a binary archive */
  void save(graphlab::oarchive& arc) const { 
    arc << nupdates << residual << pvec << A_ii;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) { 
    arc >> nupdates >> residual >> pvec >> A_ii;
  }
}; // end of vertex data


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data svdo stores the most recent error estimate.
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  /**
   * \brief The type of data on the edge;
   *
   * \li *Train:* the observed value is correct and used in training
   * \li *Validate:* the observed value is correct but not used in training
   * \li *Predict:* The observed value is not correct and should not be
   *        used in training.
   */
  enum data_role_type { TRAIN, VALIDATE, PREDICT  };

  /** \brief the observed value for the edge */
  double obs;

  /** \brief The train/validation/test designation of the edge */
  data_role_type role;

  /** \brief basic initialization */
  edge_data(double obs = 0, data_role_type role = PREDICT) :
    obs(obs), role(role) { }

}; // end of edge data


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/**
 * \brief Given a vertex and an edge return the other vertex in the
 * edge.
 */
inline graph_type::vertex_type
get_other_vertex(graph_type::edge_type& edge, 
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex

typedef double gather_type;
typedef double message_type;


#include "math.hpp" //uses vertex_data and edge_data so has to be included here
#include "printouts.hpp" // the same
typedef graphlab::omni_engine<Axb> engine_type;
engine_type * pengine = NULL;





struct linear_model_saver_U {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
  */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() > 0){
      std::string ret = boost::lexical_cast<std::string>(vertex.id()) + ") ";
      //TODO for (uint i=0; i< vertex_data::NLATENT; i++)
      //  ret += boost::lexical_cast<std::string>(vertex.data().pvec[i]) + " ";
      //  ret += "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
}; 

struct linear_model_saver_V {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
  */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() == 0){
      std::string ret = boost::lexical_cast<std::string>(-vertex.id()-SAFE_NEG_OFFSET) + ") ";
     //TODO for (uint i=0; i< vertex_data::NLATENT; i++)
     //   ret += boost::lexical_cast<std::string>(vertex.data().pvec[i]) + " ";
     //   ret += "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
}; 



/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph, 
                         const std::string& filename,
                         const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
  else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;
  // Parse the line
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  strm >> source_id >> target_id;

  // for test files (.predict) no need to read the actual rating value.
  if(role == edge_data::TRAIN || role == edge_data::VALIDATE){
    strm >> obs;
  }
  target_id = -(graphlab::vertex_id_type(target_id + SAFE_NEG_OFFSET));
                          
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(obs, role)); 
  return true; // successful load
} // end of graph_loader


void init_lanczos(graph_type * g, bipartite_graph_descriptor & info){

  if (g->num_vertices() == 0)
     logstream(LOG_FATAL)<<"Failed to load graph. Aborting" << std::endl;

  data_size = nsv + nv+1 + max_iter;
  actual_vector_len = data_size;
  if (info.is_square())
     actual_vector_len = 2*data_size;

  assert(pengine);
  pengine->map_reduce_vertices<graphlab::empty>(Axb::init_lanczos);
   // g->vertex_data(i).pvec = zeros(actual_vector_len);

  logstream(LOG_INFO)<<"Allocated a total of: " << ((double)actual_vector_len * g->num_vertices() * sizeof(double)/ 1e6) << " MB for storing vectors." << std::endl;
}

vec lanczos(bipartite_graph_descriptor & info, timer & mytimer, vec & errest, 
            const std::string & vecfile){
   

   int nconv = 0;
   int its = 1;
   int mpd = 24;
   DistMat A(info);
   DistSlicedMat U(info.is_square() ? data_size : 0, info.is_square() ? 2*data_size : data_size, true, info, "U");
   DistSlicedMat V(0, data_size, false, info, "V");
   vec alpha, beta, b;
   vec sigma = zeros(data_size);
   errest = zeros(nv);
   DistVec v_0(info, 0, false, "v_0");
   if (vecfile.size() == 0)
     v_0 = randu(size(A,2));
   PRINT_VEC2("svd->V", v_0);
/* Example Usage:
  DECLARE_TRACER(classname_someevent)
  INITIALIZE_TRACER(classname_someevent, "hello world");
  Then later on...
  BEGIN_TRACEPOINT(classname_someevent)
  ...
  END_TRACEPOINT(classname_someevent)
 */
   DistDouble vnorm = norm(v_0);
   v_0=v_0/vnorm;
   PRINT_INT(nv);

   while(nconv < nsv && its < max_iter){
     logstream(LOG_INFO)<<"Starting iteration: " << its << " at time: " << mytimer.current_time() << std::endl;
     int k = nconv;
     int n = nv;
     PRINT_INT(k);
     PRINT_INT(n);

     alpha = zeros(n);
     beta = zeros(n);

     U[k] = V[k]*A._transpose();
     orthogonalize_vs_all(U, k, alpha(0));
     PRINT_VEC3("alpha", alpha, 0);

     for (int i=k+1; i<n; i++){
       logstream(LOG_INFO) <<"Starting step: " << i << " at time: " << mytimer.current_time() << std::endl;
       PRINT_INT(i);

       V[i]=U[i-1]*A;
       orthogonalize_vs_all(V, i, beta(i-k-1));
      
       PRINT_VEC3("beta", beta, i-k-1); 
      
       U[i] = V[i]*A._transpose();
       orthogonalize_vs_all(U, i, alpha(i-k));

       PRINT_VEC3("alpha", alpha, i-k);
     }

     V[n]= U[n-1]*A;
     orthogonalize_vs_all(V, n, beta(n-k-1));
     PRINT_VEC3("beta", beta, n-k-1);

  //compute svd of bidiagonal matrix
  PRINT_INT(nv);
  PRINT_NAMED_INT("svd->nconv", nconv);
  PRINT_NAMED_INT("svd->mpd", mpd);
  n = nv - nconv;
  PRINT_INT(n);
  alpha.conservativeResize(n);
  beta.conservativeResize(n);

  PRINT_MAT2("Q",eye(n));
  PRINT_MAT2("PT",eye(n));
  PRINT_VEC2("alpha",alpha);
  PRINT_VEC2("beta",beta);
 
  mat T=diag(alpha);
  for (int i=0; i<n-1; i++)
    set_val(T, i, i+1, beta(i));
  PRINT_MAT2("T", T);
  mat a,PT;
  svd(T, a, PT, b);
  PRINT_MAT2("Q", a);
  alpha=b.transpose();
  PRINT_MAT2("alpha", alpha);
  for (int t=0; t< n-1; t++)
     beta(t) = 0;
  PRINT_VEC2("beta",beta);
  PRINT_MAT2("PT", PT.transpose());

  //estiamte the error
  int kk = 0;
  for (int i=nconv; i < nv; i++){
    int j = i-nconv;
    PRINT_INT(j);
    sigma(i) = alpha(j);
    PRINT_NAMED_DBL("svd->sigma[i]", sigma(i));
    PRINT_NAMED_DBL("Q[j*n+n-1]",a(n-1,j));
    PRINT_NAMED_DBL("beta[n-1]",beta(n-1));
    errest(i) = abs(a(n-1,j)*beta(n-1));
    PRINT_NAMED_DBL("svd->errest[i]", errest(i));
    if (alpha(j) >  tol){
      errest(i) = errest(i) / alpha(j);
      PRINT_NAMED_DBL("svd->errest[i]", errest(i));
    }
    if (errest(i) < tol){
      kk = kk+1;
      PRINT_NAMED_INT("k",kk);
    }


    if (nconv +kk >= nsv){
      printf("set status to tol\n");
      finished = true;
    }
  }//end for
  PRINT_NAMED_INT("k",kk);


  vec v;
  if (!finished){
    vec swork=get_col(PT,kk); 
    PRINT_MAT2("swork", swork);
    v = zeros(size(A,1));
    for (int ttt=nconv; ttt < nconv+n; ttt++){
      v = v+swork(ttt-nconv)*(V[ttt].to_vec());
    }
    PRINT_VEC2("svd->V",V[nconv]);
    PRINT_VEC2("v[0]",v); 
  }


INITIALIZE_TRACER(matproduct, "computing ritz eigenvectors");
   //compute the ritz eigenvectors of the converged singular triplets
  if (kk > 0){
    PRINT_VEC2("svd->V", V[nconv]);
BEGIN_TRACEPOINT(matproduct);
    mat tmp= V.get_cols(nconv,nconv+n)*PT;
    V.set_cols(nconv, nconv+kk, get_cols(tmp, 0, kk));
    PRINT_VEC2("svd->V", V[nconv]);
    PRINT_VEC2("svd->U", U[nconv]);
    tmp= U.get_cols(nconv, nconv+n)*a;
END_TRACEPOINT(matproduct);
    U.set_cols(nconv, nconv+kk,get_cols(tmp,0,kk));
    PRINT_VEC2("svd->U", U[nconv]);
  }

  nconv=nconv+kk;
  if (finished)
    break;

  V[nconv]=v;
  PRINT_VEC2("svd->V", V[nconv]);
  PRINT_NAMED_INT("svd->nconv", nconv);

  its++;
  PRINT_NAMED_INT("svd->its", its);
  PRINT_NAMED_INT("svd->nconv", nconv);
  PRINT_NAMED_INT("nv",nv);

} // end(while)

  printf(" Number of computed signular values %d",nconv);
  printf("\n");
  DistVec normret(info, nconv, false, "normret");
  DistVec normret_tranpose(info, nconv, true, "normret_tranpose");
  for (int i=0; i < nconv; i++){
    normret = V[i]*A._transpose() -U[i]*sigma(i);
    double n1 = norm(normret).toDouble();
    PRINT_DBL(n1);
    normret_tranpose = U[i]*A -V[i]*sigma(i);
    double n2 = norm(normret_tranpose).toDouble();
    PRINT_DBL(n2);
    double err=sqrt(n1*n1+n2*n2);
    PRINT_DBL(err);
    PRINT_DBL(tol);
    if (sigma(i)>tol){
      err = err/sigma(i);
    }
    PRINT_DBL(err);
    PRINT_DBL(sigma(i));
    printf("Singular value %d \t%13.6g\tError estimate: %13.6g\n", i, sigma(i),err);
  }

  if (save_vectors){
     if (nconv == 0)
       logstream(LOG_FATAL)<<"No converged vectors. Aborting the save operation" << std::endl;
 
    //TOOD
    std::cout << "Saving predictions" << std::endl;
    /*const bool gzip_output = false;
    const bool save_vertices = false;
    const bool save_edges = true;
    const size_t threads_per_machine = 1;*/
    //save the linear model
    //graph.save(predictions + ".U", linear_model_saver_U(),
		//gzip_output, save_edges, save_vertices, threads_per_machine);
    //graph.save(predictions + ".V", linear_model_saver_V(),
//		gzip_output, save_edges, save_vertices, threads_per_machine);
     
    //write_output_vector(datafile + ".singular_values", format, singular_values,false, "%GraphLab SVD Solver library. This file contains the singular values.");

    //for (int i=0; i< nconv; i++){
        //TODO
        //write_output_vector(datafile + ".U." + boost::lexical_cast<std::string>(i), format, U[i].to_vec(), false, "GraphLab v2 SVD output. This file contains eigenvector number " + boost::lexical_cast<std::string>(i) + " of the matrix U");
        //write_output_vector(datafile + ".V." + boost::lexical_cast<std::string>(i), format, V[i].to_vec(), false, "GraphLab v2 SVD output. This file contains eigenvector number " + boost::lexical_cast<std::string>(i) + " of the matrix V'");
     }
  //}
  return sigma;
}

void start_engine(){
  pengine->signal_all();
  pengine->start();
}


int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the gklanczos factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir, output_dir;
  std::string predictions;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("initial_vector", vecfile,"optional initial vector");
  clopts.attach_option("debug", debug, "Display debug output.");
  clopts.attach_option("unittest", unittest,  
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("max_iter", max_iter, "max iterations");
  clopts.attach_option("ortho_repeats", ortho_repeats, "orthogonalization iterations. 1 = low accuracy but fast, 2 = medium accuracy, 3 = high accuracy but slow.");
  clopts.attach_option("nv", nv, "Number of vectors in each iteration");
  clopts.attach_option("nsv", nsv, "Number of requested singular values to comptue"); 
  clopts.attach_option("regularization", regularization, "regularization");
  clopts.attach_option("tol", tol, "convergence threshold");
  clopts.attach_option("save_vectors", save_vectors, "save output matrices U and V.");
  clopts.attach_option("nodes", nodes, "number of rows/cols in square matrix (optional)");
  clopts.attach_option("no_edge_data", no_edge_data, "matrix is binary (optional)");

  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (nv < nsv){
    logstream(LOG_FATAL)<<"Please set the number of vectors --nv=XX, to be at least the number of support vectors --nsv=XX or larger" << std::endl;
  }

  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  
  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; 
  graph_type graph(dc, clopts);  
  graph.load(input_dir, graph_loader); 
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


  // Signal all vertices on the vertices on the left (libersvd) 
  //TODO engine.map_reduce_vertices<graphlab::empty>(svd_vertex_program::signal_left);
 

  // Run the PageRank ---------------------------------------------------------
  dc.cout() << "Running SVD (gklanczos)" << std::endl;
  dc.cout() << "(C) Code by Danny Bickson, CMU " << std::endl;
  dc.cout() << "Please send bug reports to danny.bickson@gmail.com" << std::endl;
  dc.cout() << "Time   Training    Validation" <<std::endl;
  dc.cout() << "       RMSE        RMSE " <<std::endl;
  timer.start();

  init_lanczos(&graph, info);
  init_math(&graph, info, ortho_repeats, update_function);
  if (vecfile.size() > 0){
    std::cout << "Load inital vector from file" << vecfile << std::endl;
    //TODO load_vector(vecfile, format, info, graph, 0, true, false);
  }  
 
  vec errest;
 
  vec singular_values = lanczos( info, timer, errest, vecfile);
 
  //TODO std::cout << "\t Updates: " << core.last_update_count() << " per node: " 
     //<< core.last_update_count() / core.graph().num_vertices() << std::endl;

  //vec ret = fill_output(&core.graph(), bipartite_graph_descriptor, JACOBI_X);

  //TODO write_output_vector(datafile + ".singular_values", format, singular_values,false, "%GraphLab SVD Solver library. This file contains the singular values.");


  //engine.start();  

  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
            << std::endl
            << "Final Runtime (seconds):   " << runtime 
            << std::endl
            << "Updates executed: " << engine.num_updates() << std::endl
            << "Update Rate (updates/second): " 
            << engine.num_updates() / runtime << std::endl;

  // Compute the final training error -----------------------------------------
  dc.cout() << "Final error: " << std::endl;
  engine.aggregate_now("error");

  if (unittest == 1){
    assert(errest.size() == 3);
    for (int i=0; i< errest.size(); i++)
      assert(errest[i] < 1e-30);
  }
  else if (unittest == 2){
     assert(errest.size() == 10);
    for (int i=0; i< errest.size(); i++)
      assert(errest[i] < 1e-15);
  }



  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



