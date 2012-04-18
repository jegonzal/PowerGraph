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

size_t MAX_VID = 0;
size_t NWORDS = 0;
size_t NDOCS = 0;


graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::VectorXd& vec);
graphlab::iarchive& operator>>(graphlab::iarchive& arc, Eigen::VectorXd& vec);
graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::MatrixXd& mat);
graphlab::iarchive& operator>>(graphlab::iarchive& arc, Eigen::MatrixXd& mat);

/** Vertex and edge data types **/
struct vertex_data {
  float residual; 
  float neighborhood_total; //! sum of values on edges
  uint32_t nupdates; //! the number of times the vertex was updated

  Eigen::VectorXd latent; //! vector of learned values 
  //constructor
  vertex_data() : 
    residual(std::numeric_limits<float>::max()), 
    neighborhood_total(0), nupdates(0) { }
  void randomize() {
    latent.resize(NLATENT);
    // Initialize the latent variables
    for(size_t i = 0; i < NLATENT; ++i) 
      latent(i) = graphlab::random::gaussian();
  }
  void save(graphlab::oarchive& arc) const { 
    arc << residual << neighborhood_total << nupdates << latent;        
  }
  void load(graphlab::iarchive& arc) { 
    arc >> residual >> neighborhood_total >> nupdates >> latent;
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
  void init_gather(icontext_type& context) {    
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
    Xty += edata.rating * neighbor.latent;
    XtX.triangularView<Eigen::Upper>() += 
      (neighbor.latent * neighbor.latent.transpose());
  } // end of gather

  // Merge two updates
  void merge(const als_update& other) {
    ASSERT_EQ(XtX.rows(), NLATENT);
    ASSERT_EQ(XtX.cols(), NLATENT);
    ASSERT_EQ(other.XtX.rows(), NLATENT);
    ASSERT_EQ(other.XtX.cols(), NLATENT);
    ASSERT_EQ(Xty.size(), NLATENT);
    ASSERT_EQ(other.Xty.size(), NLATENT);
    error += other.error; 
    XtX.triangularView<Eigen::Upper>() += other.XtX;
    Xty += other.Xty;
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
    // Add regularization
    for(size_t i = 0; i < NLATENT; ++i) XtX(i,i) += LAMBDA*nneighbors;
    // Solve the least squares problem using eigen ----------------------------
    const Eigen::VectorXd old_latent = vdata.latent;
    vdata.latent = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xty);
    // Compute the residual change in the latent factor -----------------------
    vdata.residual = 0;
    for(int i = 0; i < XtX.rows(); ++i)
      vdata.residual += std::fabs(old_latent(i) - vdata.latent(i));
    vdata.residual /= XtX.rows();
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
  float max_priority, rmse, max_error, residual, max_residual;
  size_t min_updates, max_updates, 
    total_vupdates, total_eupdates;
public:
  aggregator() : 
    max_priority(0), rmse(0), max_error(0), residual(0), max_residual(0),
    min_updates(-1), max_updates(0), 
    total_vupdates(0), total_eupdates(0) { }
  inline bool is_factorizable() const { return true; }
  inline edge_set gather_edges() const { return graphlab::IN_EDGES; }

  void gather(icontext_type& context, const edge_type& edge) {
    const edge_data& edata = context.const_edge_data(edge);
    const vertex_data& source_vdata = context.const_vertex_data(edge.source());
    const vertex_data& target_vdata = context.const_vertex_data(edge.target());
    rmse += (edata.error * edata.error);
    max_priority =  std::max(std::max(edata.error * source_vdata.residual, 
                                      edata.error * target_vdata.residual),
                             max_priority);
    max_error = std::max(max_error, edata.error);
  }

  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.const_vertex_data();
    residual += vdata.residual;
    max_residual = std::max(max_residual, vdata.residual);
    min_updates = std::min(min_updates, size_t(vdata.nupdates));
    max_updates = std::max(max_updates, size_t(vdata.nupdates));
    total_vupdates += vdata.nupdates;
    total_eupdates += vdata.nupdates *
      (context.num_in_edges() + context.num_out_edges());
  }

  void operator+=(const aggregator& other) { 
    max_priority = std::max(max_priority, other.max_priority);
    rmse += other.rmse; 
    max_error = std::max(max_error, other.max_error);
    residual += other.residual;
    max_residual = std::max(max_residual, other.residual);
    min_updates = std::min(min_updates, other.min_updates);
    max_updates = std::max(max_updates, other.max_updates);
    total_vupdates += other.total_vupdates;
    total_eupdates += other.total_eupdates;
  }

  void finalize(iglobal_context_type& context) {
    std::cout 
      << "results:\t" 
      << timer.current_time() << '\t'
      << std::setw(10) << max_priority << '\t'
      << std::setw(10) << sqrt( rmse / context.num_edges() ) << '\t'
      << std::setw(10) << max_error << '\t'
      << std::setw(10) << max_residual << '\t'
      << std::setw(10) << max_updates << '\t'
      << std::setw(10) << min_updates << '\t'
      << std::setw(10) << total_vupdates << '\t'
      << std::setw(10) << (double(total_vupdates) / context.num_vertices()) 
      << '\t'
      << std::setw(10) << total_eupdates << '\t'
      << std::setw(10) << (double(total_eupdates) / context.num_edges())
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

void initialize_vertex_data(graphlab::distributed_control& dc, 
                            graph_type& graph);
void make_tfidf(graph_type& graph);

void save_graph_info(const size_t procid,
                     const graph_type& graph);




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
  size_t interval = 10;
  clopts.attach_option("matrix", &matrix_dir, matrix_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D",
                       &NLATENT, NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("lambda", &LAMBDA, LAMBDA, "ALS regularization weight"); 
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
  clopts.attach_option("interval", &interval, interval,
                       "statistics reporting interval");
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
  initialize_vertex_data(dc, graph);
  //make_tfidf(graph);
  std::cout << dc.procid() << ": Finished initializign vertex data. " 
            << timer.current_time() << std::endl;



  if (savebin) graph.save(binpath, binprefix);
 
 
  std::cout << dc.procid() << ": Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts.get_ncpus());
  std::cout << dc.procid() << ": Intializing engine" << std::endl;
  engine.set_options(clopts);
  engine.initialize();
  engine.add_aggregator("error", aggregator(), interval); 

  std::cout << dc.procid() << ": Scheduling all" << std::endl;
  for(size_t lvid = 0; lvid < graph.num_local_vertices(); ++lvid) {
    // Schedule only "left side" vertices
    if(graph.l_is_master(lvid) && 
       graph.l_get_vertex_record(lvid).num_in_edges == 0) {
      engine.schedule_local(lvid, als_update(10000));
    }
  }
  dc.full_barrier();
  

  // Run the PageRank ---------------------------------------------------------
  std::cout << "Running ALS" << std::endl;
  timer.start();
  engine.start();  

  const double runtime = timer.current_time();
  std::cout << "Runtime: " << runtime << " seconds." 
            << std::endl
            << "Updates executed: " << engine.last_update_count() << std::endl
            << "Update Rate (updates/second): " 
            << engine.last_update_count() / runtime << std::endl;


  if (output) save_graph_info(dc.procid(), graph);

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
    MAX_VID = std::max(MAX_VID, std::max(target, source));
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
    std::cout << "Determining max vid: ";
    std::vector<size_t> max_vids(dc.numprocs());
    max_vids[dc.procid()] = MAX_VID;
    dc.all_gather(max_vids);
    for(size_t i = 0; i < max_vids.size(); ++i)
      MAX_VID = std::max(MAX_VID, max_vids[i]);
    std::cout << MAX_VID << std::endl;
  }
  std::cout << "Adding vertex data" << std::endl;
  for(size_t vid = dc.procid(); vid <= MAX_VID; vid += dc.numprocs()) 
    graph.add_vertex(vid, vertex_data());
} // end of load graph from directory



void initialize_vertex_data(graphlab::distributed_control& dc, 
                            graph_type& graph) {
  typedef graph_type::vertex_id_type vertex_id_type;
  typedef graph_type::edge_type edge_type;
  typedef std::pair<vertex_id_type, float> pair_type;
  std::cout << "Computing neighborhood totals" << std::endl;
  graphlab::buffered_exchange<pair_type> vertex_exchange(dc, 100000);
  graphlab::buffered_exchange<pair_type>::buffer_type recv_buffer;
  graphlab::procid_t sending_proc;
  for(size_t lvid = 0; lvid < graph.num_local_vertices(); ++lvid) {
    if(graph.l_is_master(lvid) ) {
      if(graph.l_get_vertex_record(lvid).num_out_edges > 0) ++NDOCS;
      else ++NWORDS;
      graph.get_local_graph().vertex_data(lvid).randomize();
    }
    float sum = 0;
    // Compute the sum of the neighborhood (either in or out edges)
    foreach(const edge_type& edge, graph.l_in_edges(lvid))
      sum += graph.edge_data(edge).rating;
    foreach(const edge_type& edge, graph.l_out_edges(lvid))
      sum += graph.edge_data(edge).rating;

    const pair_type rec(graph.global_vid(lvid), sum);
    const graphlab::procid_t owner = graph.l_get_vertex_record(lvid).owner;
    vertex_exchange.send(owner, rec);
    // recv any buffers if necessary
    if(lvid + 1 == graph.num_local_vertices()) vertex_exchange.flush();
    while(vertex_exchange.recv(sending_proc, recv_buffer)) {
      foreach(const pair_type& pair, recv_buffer) 
        graph.vertex_data(pair.first).neighborhood_total += pair.second;
      recv_buffer.clear();
    }    
  }
  ASSERT_TRUE(vertex_exchange.empty());
  std::cout << "Synchronizing neighborhood totals" << std::endl;
  graph.synchronize();
  {
    std::cout << "Determining NWORDS: ";
    std::vector<size_t> counts(dc.numprocs());
    counts[dc.procid()] = NWORDS;
    dc.all_gather(counts);
    NWORDS = 0;
    for(size_t i = 0; i < counts.size(); ++i) NWORDS += counts[i];
     std::cout << NWORDS << std::endl;
  }
  {
    std::cout << "Determining NDOCS: ";
    std::vector<size_t> counts(dc.numprocs());
    counts[dc.procid()] = NDOCS;
    dc.all_gather(counts);
    NDOCS = 0;
    for(size_t i = 0; i < counts.size(); ++i) NDOCS += counts[i];
     std::cout << NDOCS << std::endl;
  }
} // end of initialize_vertex_data



void make_tfidf(graph_type& graph) {
  typedef graph_type::edge_type edge_type;
  std::cout << "Updating edge values with TFIDF information." << std::endl;
  for(size_t lvid = 0; lvid < graph.num_local_vertices(); ++lvid) {
    foreach(const edge_type& edge, graph.l_in_edges(lvid)) {    
      const float words_in_doc = 
        graph.vertex_data(edge.source()).neighborhood_total;
      ASSERT_GT(words_in_doc, 0);
      const float doc_freq = 
        1 + graph.num_in_edges(edge.target());
      edge_data& edata = graph.edge_data(edge);
      edata.rating = 
        log((edata.rating / words_in_doc) * log( NDOCS / doc_freq ));
    }
  }
} // end of make tfidf




void save_graph_info(const size_t procid,
                     const graph_type& graph) {
  typedef graph_type::vertex_id_type vertex_id_type;
  typedef graph_type::edge_type edge_type;

  std::stringstream sstrm;
  sstrm << "vinfo_" << procid;
  const std::string vinfo_fname = sstrm.str();
  sstrm.str("");
  sstrm << "einfo_" << procid;
  const std::string einfo_fname = sstrm.str();
  std::cout << "[" << vinfo_fname << ", " << einfo_fname << "]" << std::endl;

  std::ofstream vfout(vinfo_fname.c_str());
  std::ofstream efout(einfo_fname.c_str());
  for(size_t lvid = 0; lvid < graph.num_local_vertices(); ++lvid) {
    const size_t gvid = graph.global_vid(lvid);
    const vertex_data& vdata = graph.vertex_data(gvid);
    vfout << gvid << '\t'
          << graph.l_is_master(lvid) << '\t'
          << vdata.residual << '\t'
          << vdata.neighborhood_total << '\t'
          << vdata.nupdates << '\t';
    ASSERT_EQ(vdata.latent.size(), NLATENT);
    for(size_t i = 0; i < NLATENT; ++i)
      vfout << vdata.latent(i) << ( (i+1) < NLATENT ? '\t' : '\n');
   
    // save edge information
    foreach(const edge_type& edge, graph.l_in_edges(lvid)) {
      const edge_data& edata = graph.edge_data(edge);
      efout << edge.source() << '\t' << edge.target() << '\t'
            << edata.rating << '\t' << edata.error << '\n';
    }
  }
  vfout.close();
  efout.close();

} // end of save graph info
