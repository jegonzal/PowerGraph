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

#include <distributed_graphlab.hpp>
#include <graphlab/util/stl_util.hpp>






#include <graphlab/macros_def.hpp>

/**
 * Define global constants for debugging.
 */
double TOLERANCE = 1e-2;
size_t NLATENT = 20;
double LAMBDA = 0.01;
size_t max_wordid = 0;
size_t max_vid = 0;
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::VectorXd& vec);
graphlab::iarchive& operator>>(graphlab::iarchive& arc, Eigen::VectorXd& vec);
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::MatrixXd& mat);
graphlab::iarchive& operator>>(graphlab::iarchive& arc, Eigen::MatrixXd& mat);

/** Vertex and edge data types **/
struct vertex_data {
  float residual; 
  uint32_t nupdates; //! the number of times the vertex was updated
  Eigen::VectorXd latent; //! vector of learned values 
  //constructor
  vertex_data() : 
    residual(std::numeric_limits<float>::max()), 
    nupdates(0) { 
    latent.resize(NLATENT);
    // Initialize the latent variables
    for(size_t i = 0; i < NLATENT; ++i) 
      latent(i) = graphlab::random::gaussian();
  }
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
  float rating, error;
  edge_data(const float& rating = 0) :
    rating(rating), error(std::numeric_limits<float>::max()) { }
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
    //if(edata.is_test) return;
    const vertex_id_type neighbor_id = context.vertex_id() == edge.target()?
      edge.source() : edge.target();
    const vertex_data& neighbor = context.const_vertex_data(neighbor_id);
    ASSERT_EQ(neighbor.latent.size(), NLATENT);
    for(int i = 0; i < XtX.rows(); ++i) {
      Xty(i) += neighbor.latent(i) * edata.rating;
      // Compute the upper triangular component of XtX
      for(int j = i; j < XtX.rows(); ++j) 
        XtX(j,i) += neighbor.latent(i) * neighbor.latent(j);
    }
  } // end of gather

  // Merge two updates
  void merge(const als_update& other) {
    ASSERT_EQ(XtX.rows(), NLATENT);
    ASSERT_EQ(XtX.cols(), NLATENT);
    ASSERT_EQ(other.XtX.rows(), NLATENT);
    ASSERT_EQ(other.XtX.cols(), NLATENT);
    ASSERT_EQ(Xty.size(), NLATENT);
    ASSERT_EQ(other.Xty.size(), NLATENT);
    error += other.error; XtX += other.XtX; Xty += other.Xty;
  } // end of merge


  // Update the center vertex
  void apply(icontext_type& context) {
    // Get and reset the vertex data
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    ASSERT_EQ(vdata.latent.size(), NLATENT);
    // Determine the number of neighbors.  Each vertex has only in or
    // out edges depending on which side of the graph it is located
    const size_t nneighbors = context.num_in_edges() + context.num_out_edges();
    if(nneighbors == 0) return;
    // Fill in the lower triangular components of XtX
    for(size_t i = 0; i < NLATENT; ++i) 
      for(size_t j = i+1; j < NLATENT; ++j) XtX(i,j) = XtX(j,i);
    // Add regularization
    for(int i = 0; i < XtX.rows(); ++i) XtX(i,i) += LAMBDA*nneighbors;
    // Solve the least squares problem using eigen ----------------------------
    const Eigen::VectorXd old_latent = vdata.latent;
    vdata.latent = XtX.ldlt().solve(Xty);
    // Compute the residual change in the latent factor -----------------------
    vdata.residual = 0;
    for(int i = 0; i < XtX.rows(); ++i)
      vdata.residual += std::fabs(old_latent(i) - vdata.latent(i));
    vdata.residual /= XtX.rows();
    // if(context.vertex_id() % 1000 == 0) {
    //   std::cout << context.vertex_id() << ": r = " 
    //             << vdata.residual << "; [";
    //   for(int i = 0; i < XtX.rows(); ++i)
    //     std::cout << vdata.latent(i) << " ";
    //   std::cout << "]" << std::endl;
    // } // display output
  } // end of apply

  void scatter(icontext_type& context, const edge_type& edge) {
    const vertex_data& vdata = context.const_vertex_data();
    const vertex_id_type neighbor_id = context.vertex_id() == edge.target()?
      edge.source() : edge.target();
    const vertex_data& neighbor = context.const_vertex_data(neighbor_id);
    edge_data& edata = context.edge_data(edge);
    // Compute the prediction on the edge
    const double pred = vdata.latent.dot(neighbor.latent);
    const double error = std::fabs(edata.rating - pred);
    edata.error = error;
    // std::cout << edata.rating << '\t' << pred 
    //           << '\t' << vdata.residual << std::endl;
    // Reschedule neighbors ------------------------------------------------
    if( (error * vdata.residual) > TOLERANCE ) 
      context.schedule(neighbor_id, als_update(error * vdata.residual));
  } // end of scatter
  void save(graphlab::oarchive& arc) const { arc << XtX << Xty << error; }
  void load(graphlab::iarchive& arc) { arc >> XtX >> Xty >> error; }  
}; // end of class ALS update


graphlab::timer timer;


class aggregator :
  public graphlab::iaggregator<graph_type, als_update, aggregator>, 
  public graphlab::IS_POD_TYPE {
private:
  float rmse, max_error, residual, max_residual;
  size_t nedges, min_updates, max_updates, total_updates;
public:
  aggregator() : 
    rmse(0), max_error(0), residual(0), max_residual(0),
    nedges(0), min_updates(-1), max_updates(0), total_updates(0) { }
  inline bool is_factorizable() const { return true; }

  void gather(icontext_type& context, const edge_type& edge) {
    const edge_data& edata = context.const_edge_data(edge);
    rmse += (edata.error * edata.error);
    ++nedges;
    max_error = std::max(max_error, edata.error);
  }

  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.const_vertex_data();
    residual += vdata.residual;
    max_residual = std::max(max_residual, vdata.residual);
    min_updates = std::min(min_updates, size_t(vdata.nupdates));
    max_updates = std::max(max_updates, size_t(vdata.nupdates));
    total_updates += vdata.nupdates;
  }

  void operator+=(const aggregator& other) { 
    rmse += other.rmse; 
    nedges += other.nedges;
    max_error = std::max(max_error, other.max_error);
    residual += other.residual;
    max_residual = std::max(max_residual, other.residual);
    min_updates = std::min(min_updates, other.min_updates);
    max_updates = std::max(max_updates, other.max_updates);
    total_updates += other.total_updates;
  }

  void finalize(iglobal_context_type& context) {
    std::cout 
      << "results:\t" 
      << timer.current_time() << '\t'
      << std::setw(10) << sqrt( rmse / nedges ) << '\t'
      << std::setw(10) << max_error << '\t'
      << std::setw(10) << max_residual << '\t'
      << std::setw(10) << max_updates << '\t'
      << std::setw(10) << min_updates << '\t'
      << std::setw(10) << total_updates << '\t'
      << std::setw(10) << (double(total_updates) / context.num_vertices()) 
      << std::endl;
  }
}; // end of  aggregator





#ifdef FSCOPE
typedef graphlab::distributed_fscope_engine<graph_type, als_update> engine_type;
#elif SYNC
typedef graphlab::distributed_synchronous_engine<graph_type, als_update> engine_type;
#else
typedef graphlab::distributed_engine<graph_type, als_update> engine_type;
#endif



//! Graph loading code 
void load_graph_dir(graphlab::distributed_control& dc,
                    graph_type& graph, const std::string& matrix_dir,
                    const size_t procid, const size_t numprocs);
void load_graph_file(graph_type& graph, const std::string& fname);





int main(int argc, char** argv) {
  //global_logger().set_log_level(LOG_DEBUG);
  //global_logger().set_log_to_console(true);

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

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
  


  std::cout << dc.procid() << ": Loading graph." << std::endl;
  timer.start();
  graph_type graph(dc, clopts);
  load_graph_dir(dc, graph, matrix_dir, dc.procid(), dc.numprocs());
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
  std::cout << dc.procid() << ": Initializign vertex data. " 
            << timer.current_time() << std::endl;




  if (savebin) graph.save(binpath, binprefix);
 
 
  std::cout << dc.procid() << ": Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts.get_ncpus());
  std::cout << dc.procid() << ": Intializing engine" << std::endl;
  engine.set_options(clopts);
  engine.initialize();
  engine.add_aggregator("error", aggregator(), 10); 
  std::cout << "Precomputing the error" << std::endl;
  // engine.aggregate_now("error");

  std::cout << dc.procid() << ": Scheduling all" << std::endl;
  for(size_t lvid = 0; lvid < graph.num_local_vertices(); ++lvid) {
    if(graph.l_is_master(lvid) && graph.global_vid(lvid) < max_wordid)  
      engine.schedule_local(lvid, als_update(10000));
  }
  //  engine.schedule_all(als_update(10000));
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








graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::VectorXd& vec) {
  typedef Eigen::VectorXd::Index index_type;
  const index_type size = vec.size();
  arc << size;
  graphlab::serialize(arc, vec.data(), size * sizeof(double));
  return arc;
} // end of save vector

graphlab::iarchive& operator>>(graphlab::iarchive& arc, Eigen::VectorXd& vec) {
  typedef Eigen::VectorXd::Index index_type;
  index_type size = 0;
  arc >> size;
  vec.resize(size);
  graphlab::deserialize(arc, vec.data(), size * sizeof(double));
  return arc;
} // end of save vector


graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::MatrixXd& mat) {
  typedef Eigen::MatrixXd::Index index_type;
  const index_type rows = mat.rows();
  const index_type cols = mat.cols();
  arc << rows << cols;
  graphlab::serialize(arc, mat.data(), rows*cols*sizeof(double));
  return arc;
} // end of save vector

graphlab::iarchive& operator>>(graphlab::iarchive& arc,  Eigen::MatrixXd& mat) {
  typedef Eigen::MatrixXd::Index index_type; 
  typedef Eigen::MatrixXd::Scalar scalar_type;
  index_type rows=0, cols=0;
  arc >> rows >> cols;
  mat.resize(rows,cols);
  graphlab::deserialize(arc, mat.data(), rows*cols*sizeof(scalar_type));
  return arc;
} // end of save vector


template<typename Stream>
void load_graph_stream(graph_type& graph, Stream& fin) {
  // Loop over the contents
  while(fin.good()) {
    // Load a vertex
    size_t source = 0, target = 0;
    float value = 0;
    fin >> source >> target >> value; 
    if(!fin.good()) break;
    ASSERT_LT(target, source);
    max_wordid = std::max(target, max_wordid);
    max_vid = std::max(max_vid, std::max(target, source));
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
    load_graph_stream(graph, fin);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } else {
    std::ifstream in_file(fname.c_str(), 
                          std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    if (gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good()) {
      logstream(LOG_FATAL) << "Error opening file:" << fname << std::endl;
    }
    load_graph_stream(graph, fin);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
} // end of load graph from file


void load_graph_dir(graphlab::distributed_control& dc,
                    graph_type& graph, const std::string& matrix_dir,
                    const size_t procid, const size_t numprocs) {
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
    if (i % numprocs == procid) {
      std::cout << "Loading graph from structure file: " 
                << graph_files[i] << std::endl;
      load_graph_file(graph, graph_files[i]);
    }
  }
  {
    std::cout << "Determining max wordid: ";
    std::vector<size_t> max_wordids(dc.numprocs());
    max_wordids[dc.procid()] = max_wordid;
    dc.all_gather(max_wordids);
    for(size_t i = 0; i < max_wordids.size(); ++i)
      max_wordid = std::max(max_wordid, max_wordids[i]);
    std::cout << max_wordid << std::endl;
  }
  {
    std::cout << "Determining max vid: ";
    std::vector<size_t> max_vids(dc.numprocs());
    max_vids[dc.procid()] = max_vid;
    dc.all_gather(max_vids);
    for(size_t i = 0; i < max_vids.size(); ++i)
      max_vid = std::max(max_vid, max_vids[i]);
    std::cout << max_vid << std::endl;
  }
  std::cout << "Adding vertex data" << std::endl;
  for(size_t vid = dc.procid(); vid <= max_vid; vid += dc.numprocs()) 
    graph.add_vertex(vid, vertex_data());
 

} // end of load graph from directory
