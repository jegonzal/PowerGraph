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


#include <vector>

#include <distributed_graphlab.hpp>
#include <graphlab/util/stl_util.hpp>


#include <graphlab/macros_def.hpp>


typedef uint32_t doc_id_type;
typedef uint16_t topic_id_type;
typedef int      count_type;

const topic_id_type NULL_TOPIC(-1);
typedef std::vector<count_type> factor_type;
typedef std::vector<topic_id_type> token_asg_type;


// Global variables
size_t NTOPICS = 100;
size_t NWORDS = 0;
size_t NDOCS = 0;
double ALPHA = 0.01;
double BETA = 1;
std::vector< std::string > dictionary;
std::vector< graphlab::atomic<count_type> > global_factor;






void operator+=(factor_type& a, const factor_type& b) {
  ASSERT_EQ(a.size(), b.size());
  for(size_t i = 0; i < a.size(); ++i) a[i] += b[i];
}

/**
 * The vertex data type
 */
struct vertex_data {
  factor_type factor;
  vertex_data() : factor(NTOPICS)  { }
  void save(graphlab::oarchive& arc) const { arc << factor; }
  void load(graphlab::iarchive& arc) { arc >> factor; }
};


/**
 * The edge data type
 */
struct edge_data {
  token_asg_type asgs;
  edge_data(count_type count = 0) : asgs(count, NULL_TOPIC) { }  
  void save(graphlab::oarchive& arc) const { arc << asgs; }
  void load(graphlab::iarchive& arc) { arc >> asgs; }
};


/**
 * The graph type;
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;






class lda_cgs :
  public graphlab::iupdate_functor<graph_type, lda_cgs> {
private:
  factor_type factor;
  typedef std::pair<token_asg_type, factor_type> edge_pair_type;
  typedef std::map<vertex_id_type, edge_pair_type> nbr_info_map_type;
  nbr_info_map_type nbr_info;

public:
  lda_cgs() { }
  void save(graphlab::oarchive& arc) const { arc << factor << nbr_info; }
  void load(graphlab::iarchive& arc) { arc >> factor >> nbr_info; }  
  edge_set gather_edges() const { return graphlab::ALL_EDGES; }
  edge_set scatter_edges() const { return graphlab::ALL_EDGES; }
  bool is_factorizable() const { return true; }

  // Reset the accumulator before running the gather
  void init_gather(icontext_type& context) {  } // end of init gather

  void gather(icontext_type& context, const edge_type& edge) {
    // Assume: Document ----> word
    const bool is_word = context.vertex_id() == edge.target();
    const vertex_id_type neighbor_id = context.vertex_id() == edge.target()?
      edge.source() : edge.target();    
    const vertex_data& neighbor_vdata = context.const_vertex_data(neighbor_id); 
    const edge_data& edata = context.const_edge_data(edge);
    if(is_word) { // accumulating the new count
      if(factor.empty()) factor.resize(NTOPICS);
      foreach(const topic_id_type asg, edata.asgs) 
        if(asg != NULL_TOPIC) ++factor[asg];
    } else { // This is a document vertex
      // Collect the information about the neighboring word and document
      nbr_info[neighbor_id] = 
        edge_pair_type(edata.asgs, neighbor_vdata.factor);
    }
  } // end of gather

  void merge(const lda_cgs& other) {    
    if(!other.factor.empty()) {
      if(factor.empty()) factor.resize(NTOPICS);
      factor += other.factor;
    }
    nbr_info.insert(other.nbr_info.begin(), other.nbr_info.end());
  } // end of merge

  void apply(icontext_type& context) {
    vertex_data& vdata = context.vertex_data();
    const size_t nneighbors = context.num_in_edges() + context.num_out_edges();
    ASSERT_TRUE(nneighbors > 0);
    // Assume: Document ----> word    
    const bool is_word = context.num_in_edges() > 0;
    if(is_word) {
      vdata.factor = factor;
      factor_type().swap(factor); // Clear the factor
      ASSERT_EQ(nbr_info.size(), 0);
    } else {
      factor_type& doc_factor = vdata.factor;
      if(doc_factor.size() != NTOPICS) doc_factor.resize(NTOPICS);
      // run the actual gibbs sampling routine
      typedef nbr_info_map_type::value_type pair_type;
      foreach(pair_type& pair, nbr_info) {
        token_asg_type& asgs = pair.second.first;
        factor_type& word_factor = pair.second.second;
        if(word_factor.size() != NTOPICS) word_factor.resize(NTOPICS);
        // Resample the topics
        foreach(topic_id_type& asg, asgs) {
          if(asg != NULL_TOPIC) { // construct the cavity
            --doc_factor[asg];
            --word_factor[asg];
            --global_factor[asg];
          }
          std::vector<double> prob(NTOPICS);
          // double normalizer = 0;
          for(size_t t = 0; t < NTOPICS; ++t) {
            const double n_dt = std::max(doc_factor[t], count_type(0));
            const double n_wt = std::max(word_factor[t], count_type(0));
            const double n_t  = std::max(count_type(global_factor[t]), count_type(0));
            prob[t] = (ALPHA + n_dt) * (BETA + n_wt) / (BETA * NWORDS + n_t);
            // normalizer += prob[t];
          }
          // ASSERT_GT(normalizer, 0);
          // for(size_t t = 0; t < ntopics; ++t) prob[t] /= normalizer;
          asg = graphlab::random::multinomial(prob);
          doc_factor[asg]++;
          word_factor[asg]++;                    
          global_factor[asg]++;
        } // End of loop over each token
      } // end of loop over neighbors      
    }
  } // end of apply

  void scatter(icontext_type& context, const edge_type& edge) {
    // Assume: Document ----> word
    const bool is_doc = context.vertex_id() == edge.source();
    const vertex_id_type neighbor_id = context.vertex_id() == edge.target()?
      edge.source() : edge.target();    
    if(is_doc) {
      edge_data& edata = context.edge_data(edge);
      // update the new assignment
      edata.asgs.swap(nbr_info[neighbor_id].first);      
    }
    // reschedule the neighbor to be run again in the future
    context.schedule(neighbor_id, lda_cgs());
  } // end of gather



}; // end of lad




#ifdef FSCOPE
typedef graphlab::distributed_fscope_engine<graph_type, lda_cgs> engine_type;
#elif SYNC
typedef graphlab::distributed_synchronous_engine<graph_type, lda_cgs> engine_type;
#else
typedef graphlab::distributed_engine<graph_type, lda_cgs> engine_type;
#endif



//! Graph loading code 

/** populate the global dictionary */
void load_dictionary(const std::string& fname);

void load_graph_dir(graphlab::distributed_control& dc,
                    graph_type& graph, const std::string& matrix_dir,
                    const size_t procid, const size_t numprocs);
void load_graph_file(graph_type& graph, const std::string& fname);





int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Run the asynchronous collapsed Gibbs Sampler.";
  graphlab::command_line_options clopts(description);
  clopts.use_distributed_options();
  std::string matrix_dir; 
  std::string dictionary_file;
  size_t interval = 10;
  clopts.attach_option("matrix", &matrix_dir, matrix_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("dictionary", &dictionary, dictionary,
                       "The file containing the list of unique words");
  clopts.add_positional("dictionary");

  clopts.attach_option("ntopics", &NTOPICS, NTOPICS,
                       "Number of topics to use.");
  clopts.attach_option("alpha", &ALPHA, ALPHA,
                       "The document hyper-prior");
  clopts.attach_option("beta", &ALPHA, ALPHA,
                       "The word hyper-prior");
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
  
  std::cout << dc.procid() << ": Loading the dictionary." << std::endl;
  load_dictionary(dictionary_file);
  dc.barrier();

  std::cout << dc.procid() << ": Loading graph." << std::endl;
  graphlab::timer timer;
  timer.start();
  graph_type graph(dc, clopts);
  std::cout << "Loading text matrix file" << std::endl;
  load_graph_dir(dc, graph, matrix_dir, dc.procid(), dc.numprocs());
  std::cout << dc.procid() << ": Finalizing graph." << std::endl;
  graph.finalize();
  std::cout << dc.procid() << ": Finished in " 
            << timer.current_time() << std::endl;  
  
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
  dc.barrier();


 
  std::cout << dc.procid() << ": Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts.get_ncpus());
  std::cout << dc.procid() << ": Intializing engine" << std::endl;
  engine.set_options(clopts);
  engine.initialize();



  std::cout << dc.procid() << ": Scheduling all documents" << std::endl;
  std::vector<graph_type::vertex_id_type> vtxs;
  vtxs.reserve(graph.num_local_vertices());
  for(graph_type::lvid_type lvid = 0; lvid < graph.num_local_vertices(); 
      ++lvid) {
    if(graph.l_is_master(lvid) && 
       graph.l_get_vertex_record(lvid).num_in_edges == 0)
      vtxs.push_back(lvid);             
  }
  graphlab::random::shuffle(vtxs.begin(), vtxs.end());
  foreach(graph_type::lvid_type lvid, vtxs)
    engine.schedule_local(lvid, lda_cgs());  
  dc.full_barrier();
  

  std::cout << "Running The Collapsed Gibbs Sampler" << std::endl;
  timer.start();
  engine.start();  

  const double runtime = timer.current_time();
  if(dc.procid() == 0) {
    std::cout << "----------------------------------------------------------"
              << std::endl;
    std::cout << "Final Runtime (seconds):   " << runtime 
              << std::endl
              << "Updates executed: " << engine.last_update_count() << std::endl
              << "Update Rate (updates/second): " 
              << engine.last_update_count() / runtime << std::endl;
  }

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main



























/** populate the global dictionary */
void load_dictionary(const std::string& fname) {
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
      logstream(LOG_FATAL) << "Error loading dictionary: "
                           << fname << std::endl;
    }
    std::string term;
    while(std::getline(fin,term) && fin.good()) dictionary.push_back(term);
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
    std::string term;
    while(std::getline(fin,term) && fin.good()) dictionary.push_back(term);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
  NWORDS = dictionary.size();
  std::cout << "Dictionary size: " << NWORDS << std::endl;
} // end of load dictionary




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
    graph.add_edge(source + NWORDS, target, edge_data(value));
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
} // end of load graph from directory




