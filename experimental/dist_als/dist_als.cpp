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


#include <Eigen/Dense>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

/**
 * Define global constants for debugging.
 */
double TOLERANCE = 1e-2;
size_t NLATENT = 20;
double LAMBDA = 0.065;

void save(graphlab::oarchive& arc, const Eigen::VectorXd& vec);
void load(graphlab::iarchive& arc, Eigen::VectorXd& vec);
void save(graphlab::oarchive& arc, const Eigen::MatrixXd& mat);
void load(graphlab::iarchive& arc, Eigen::MatrixXd& mat);

/** Vertex and edge data types **/
struct vertex_data {
  uint32_t nupdates; //! the number of times the vertex was updated
  float residual; 
  Eigen::VectorXd latent; //! vector of learned values 
  //constructor
  vertex_data() : 
    residual(std::numeric_limits<float>::max()), 
    nupdates(0) { latent.resize(NLATENT); latent.setZero(); }
  void save(graphlab::oarchive& arc) const { 
    arc << residual << nupdates << latent;        
  }
  void load(graphlab::iarchive& arc) { 
    arc >> residual >> nupdates >> latent;
  }
}; // end of vertex data

/**
 * The edge data is just an observation float
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  bool is_test; float rating;
  edge_data(const float& rating = 0, const bool& is_test = false) : 
    rating(rating), is_test(is_test) { }
}; // end of edge data

/**
 * The graph type;
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



/***
 * UPDATE FUNCTION
 */
class als_update : 
  public graphlab::iupdate_functor<graph_type, als_update> {
  double error;
  Eigen::MatrixXd XtX;
  Eigen::VectorXd Xty;
public:
  als_update(double error = 0) : error(error) { }
  double priority() const { return error; }
  void operator+=(const als_update& other) { error += other.error; }
  consistency_model gather_consistency() { return graphlab::EDGE_CONSISTENCY; }
  consistency_model scatter_consistency() { return graphlab::EDGE_CONSISTENCY; }
  edge_set gather_edges() const { return graphlab::ALL_EDGES; }
  edge_set scatter_edges() const { return graphlab::ALL_EDGES; }
  bool is_factorizable() const { return true; }

  // Reset the accumulator before running the gather
  void init_gather(iglobal_context_type& context) {    
    XtX.resize(NLATENT, NLATENT); XtX.setZero();
    Xty.resize(NLATENT); Xty.setZero();
  } // end of init gather

  // Gather the XtX and Xty matrices from the neighborhood of this vertex
  void gather(icontext_type& context, const edge_type& edge) {
    // Get the edge data and skip the edge if it is for testing only
    const edge_data& edata = context.const_edge_data(edge);
    if(edata.is_test) return;
    const vertex_id_type neighbor_id = context.vertex_id() == edge.target()?
      edge.source() : edge.target();
    const vertex_data& neighbor = context.const_vertex_data(neighbor_id);
    for(int i = 0; i < XtX.rows(); ++i) {
      Xty(i) += neighbor.latent(i) * edata.rating;
      // Compute the upper triangular component of XtX
      for(int j = i; j < XtX.rows(); ++j) 
        XtX(j,i) += neighbor.latent(i) * neighbor.latent(j);
    }
  } // end of gather

  // Update the center vertex
  void apply(icontext_type& context) {
    // Get and reset the vertex data
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    // Determine the number of neighbors.  Each vertex has only in or
    // out edges depending on which side of the graph it is located
    const size_t nneighbors = context.num_in_edges() + context.num_out_edges();
    if(nneighbors == 0) return;
    // Fill in the lower triangular components of XtX
    for(size_t i = 0; i < NLATENT; ++i) 
      for(size_t j = i+1; j < NLATENT; ++j) XtX(i,j) = XtX(j,i);
    // Add regularization
    for(int i = 0; i < XtX.rows(); ++i) XtX(i,i) += (LAMBDA)*nneighbors;
    // Solve the least squares problem using eigen ----------------------------
    const vec old_latent = vdata.latent;
    vdata.latent = XtX.ldlt().solve(Xty);
    // Compute the residual change in the latent factor -----------------------
    vdata.residual = 0;
    for(int i = 0; i < NLATENT; ++i)
      vdata.residual += std::fabs(old_latent(i) - vdata.latent(i));
    vdata.residual /= XtX.rows();
  } // end of apply

  void scatter(icontext_type& context, const edge_type& edge) {
    const vertex_data& vdata = context.const_vertex_data();
    const vertex_id_type neighbor_id = context.vertex_id() == edge.target()?
      edge.source() : edge.target();
    const vertex_data& neighbor = context.const_vertex_data(neighbor_id);
    const edge_data& edata = context.const_edge_data(edge);
    // Compute the prediction on the edge
    const double pred = vdata.latent.dot(neighbor.latent);
    const double error = std::fabs(edata.rating - pred);
    // Reschedule neighbors ------------------------------------------------
    if( error > TOLERANCE && vdata.residual > TOLERANCE) 
      context.schedule(neighbor_id, als_update(error * vdata.residual));
  } // end of scatter
  void save(graphlab::oarchive& arc) const { arc << XtX << Xty << error; }
  void load(graphlab::iarchive& arc) { arc >> XtX >> Xty >> error; }  
}; // end of class ALS update


#ifdef FSCOPE
typedef graphlab::distributed_fscope_engine<graph_type, als_update> engine_type;
#else
typedef graphlab::distributed_engine<graph_type, als_update> engine_type;
#endif



//! Graph loading code 
void load_graph_dir(graph_type& graph, const std::string& matrix_dir);
void load_graph_file(graph_type& graph, const std::string& fname);





int main(int argc, char** argv) {
  //global_logger().set_log_level(LOG_DEBUG);
  //global_logger().set_log_to_console(true);

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
  
  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  clopts.use_distributed_options();
  std::string matrix_dir; 
  bool savebin = false;
  std::string binpath = "./";
  std::string binprefix = "als_graph";
  bool output = false;
  clopts.attach_option("matrix", &matrix_dir, matrix_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
   clopts.attach_option("D",
                       &NLATENT, NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("LAMBDA", &LAMBDA, LAMBDA, "ALS regularization weight"); 
  clopts.attach_option("tol",
                       &TOLERANCE, TOLERANCE,
                       "residual termination threshold");
  clopts.attach_option("output", &output, output,
                       "Output results");
  clopts.attach_option("savebin", &savebin, savebin,
                       "Option to save the graph as binary\n");
  clopts.attach_option("binpath", &binpath, binpath,
                       "The path for save binary file\n");
  clopts.attach_option("binprefix", &binprefix, binprefix,
                       "The prefix for load/save binary file\n");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }


  std::cout << dc.procid() << ": Loading graph." << std::endl;
  graphlab::timer timer; timer.start();
  graph_type graph(dc, clopts);
  load_graph_dir(graph, matrix_dir);
  std::cout << dc.procid() << ": Finalizing graph." << std::endl;
  graph.finalize();
  std::cout << dc.procid() << ": Finished in " << timer.current_time() << std::endl;


  if(dc.procid() == 0){
    std::cout
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
  }

  if (savebin) graph.save(binpath, binprefix);
 
 
  std::cout << dc.procid() << ": Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts.get_ncpus());
  std::cout << dc.procid() << ": Intializing engine" << std::endl;
  engine.set_options(clopts);
  engine.initialize();
  std::cout << dc.procid() << ": Scheduling all" << std::endl;
  engine.schedule_all(als_update(10000));
  dc.full_barrier();
  
  // Run the PageRank ---------------------------------------------------------
  std::cout << "Running ALS" << std::endl;
  timer.start();
  engine.start();  

  const double runtime = timer.current_time();
  std::cout << "Graphlab finished, runtime: " << runtime << " seconds." 
            << std::endl
            << "Updates executed: " << engine.last_update_count() << std::endl
            << "Update Rate (updates/second): " 
            << engine.last_update_count() / runtime << std::endl;


  if (output) {
    std::string fname = "results_";
    fname = fname + graphlab::tostr(size_t(dc.procid()));
    std::ofstream fout(fname.c_str());
    for (size_t i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
      if (graph.l_get_vertex_record(i).owner == dc.procid()) {
        fout << graph.l_get_vertex_record(i).gvid << "\t" 
             << graph.l_get_vertex_record(i).num_in_edges + 
          graph.l_get_vertex_record(i).num_out_edges << "\t" 
             << graph.get_local_graph().vertex_data(i).residual << "\t"
             << graph.get_local_graph().vertex_data(i).nupdates << "\n";
      }
    }
  }
  
  if (output) {
    std::string fname = "results_local_";
    fname = fname + graphlab::tostr((size_t)dc.procid());
    std::ofstream fout(fname.c_str());
    for (size_t i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
      fout << graph.l_get_vertex_record(i).gvid << "\t" 
           << graph.l_get_vertex_record(i).num_in_edges + 
        graph.l_get_vertex_record(i).num_out_edges << "\t" 
           << graph.get_local_graph().vertex_data(i).residual << "\t"
           << graph.get_local_graph().vertex_data(i).nupdates << "\n";
    }
  } 

  if (output) {
    std::string fname = "adj_";
    fname = fname + graphlab::tostr((size_t)dc.procid());
    std::ofstream fout(fname.c_str());
    typedef graph_type::local_graph_type::edge_type etype;
    for (size_t i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
      foreach(graph_type::edge_type e, graph.l_in_edges(i)) {
        fout << e.source() << "\t" << e.target() << "\n";
      }
    }
  } 

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main








void save(graphlab::oarchive& arc, const Eigen::VectorXd& vec) {
  typedef Eigen::VectorXd::Index index_type;
  const index_type size = vec.size();
  arc << size;
  graphlab::archive_detail::serialize(arc, vec.data(), 
                                      size * sizeof(double));
} // end of save vector


void load(graphlab::iarchive& arc, Eigen::VectorXd& vec) {
  typedef Eigen::VectorXd::Index index_type;
  index_type size = 0;
  arc >> size;
  vec.resize(size);
  graphlab::archive_detail::deserialize(arc, vec.data(), 
                                        size * sizeof(double));
} // end of save vector


void save(graphlab::oarchive& arc, const Eigen::MatrixXd& mat) {
  typedef Eigen::MatrixXd::Index index_type;
  const index_type rows = mat.rows();
  const index_type cols = mat.cols();
  arc << rows << cols;
  graphlab::archive_detail::serialize(arc, mat.data(), 
                                      rows*cols*sizeof(double));
} // end of save vector


void load(graphlab::iarchive& arc,  Eigen::MatrixXd& mat) {
  typedef Eigen::MatrixXd::Index index_type; 
  typedef Eigen::MatrixXd::Scalar scalar_type;
  index_type rows=0, cols=0;
  arc >> rows >> cols;
  mat.resize(rows,cols);
  graphlab::archive_detail::deserialize(arc, mat.data(), 
                                        rows*cols*sizeof(scalar_type));
} // end of save vector


template<typename Stream>
void load_graph_stream(graph_type& graph, Stream& fin) {
  // Loop over the contents
  while(fin.good()) {
    // Load a vertex
    graph_type::vertex_id_type source = 0, target = 0;
    float value = 0;
    try { fin >> source >> target >> value; } catch ( ... ) { 
      logstream(LOG_WARNING) 
        << "Error reading source." << std::endl; return false;
    }
    if(!fin.good()) break;
    graph.add_edge(source, target, edge_data(value));
  } // end of loop over file
} // end of load graph from stream;

void load_graph_file(graph_type& graph, const std::string& fname) {
  const bool gzip = boost::ends_with(fname, ".gz");
  // test to see if the graph_dir is an hadoop path
  if(boost::starts_with(fname, "hdfs://")) {
    graphlab::hdfs hdfs;
    graphlab::hdfs::fstream in_file(hdfs, fname);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    fin.set_auto_close(false);
    if(gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good()) {
      logstream(LOG_FATAL) << "Error opening file:" << fname << std::endl;
    }
    load_graph_from_stream(graph, fin);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
    return success;
  } else {
    std::ifstream in_file(fname.c_str(), 
                          std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    if (gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good()) {
      logstream(LOG_FATAL) << "Error opening file:" << fname << std::endl;
    }
    const bool success = load_structure_from_stream(fin, format, graph);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
} // end of load graph from file


void load_graph_dir(graph_type& graph, const std::string& matrix_dir) {
  std::vector<std::string> graph_files;
  if(boost::starts_with(matrix_dir, "hdfs://")) {
    graphlab::hdfs hdfs;
    graph_files = hdfs.list_files(matrix_dir);
  } else {
    graphlab::fs_util::list_files_with_prefix(matrix_dir, "", graph_files);
    for(size_t i = 0; i < graph_files.size(); ++i)
      graph_files[i] = matrix_dir + graph_files[i];
  }
  std::sort(graph_files.begin(), graph_files.end());
  
  for(size_t i = 0; i < graph_files.size(); ++i) {
    if (i % dc.numprocs() == dc.procid()) {
      std::cout << "Loading graph from structure file: " 
                << graph_files[i] << std::endl;
      load_graph_from_file(graph, fname);
      ASSERT_TRUE(success);
    }
  }
} // end of load graph from directory
