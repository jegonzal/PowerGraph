/*
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
#include <string>
#include <fstream>

#include <graphlab.hpp>
// #include <graphlab/macros_def.hpp>

// Global random reset probability
double RESET_PROB = 0.15;

double TOLERANCE = 1.0E-2;

size_t ITERATIONS = 0;

bool USE_DELTA_CACHE = false;

// The vertex data is just the pagerank value (a double)
//typedef double vertex_data_type;



//struct vertex_t : public graphlab::IS_POD_TYPE {
struct vertex_t {
  int count;
  bool is_doc; //
  //std::string str;
  std::vector<std::string> str_vec;


  vertex_t(std::vector<std::string>& vec) :
    count(0){str_vec=vec;is_doc=true;}

  vertex_t() :
    count(0){is_doc=false;}

  void save(graphlab::oarchive& arc) const { 
    arc << count <<str_vec <<is_doc;         
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) { 
    arc >> count >> str_vec >>is_doc;
  }

}; // end of edge data

typedef vertex_t vertex_data_type;

// There is no edge data in the pagerank application
typedef graphlab::empty edge_data_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;


inline graph_type::vertex_type
get_other_vertex(graph_type::edge_type& edge, 
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex

/*
 * A simple function used by graph.transform_vertices(init_vertex);
 * to initialize the vertes data.
 */



void init_vertex(graph_type::vertex_type& vertex) { vertex.data().count = 0; }


inline bool graph_loader(graph_type& graph,
                         const std::string& filename,
                         const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;

  graph_type::vertex_id_type source_id(-1), target_id(-1);
  std::vector<std::string> str_vec;

  if(boost::ends_with(filename,".edge")){
    const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(source_id) = qi::_1] >> -qi::char_(',') >>
      qi::ulong_[phoenix::ref(target_id) = qi::_1] 
      )
     ,
     //  End grammar
     ascii::space); 
      
    if(!success) return false;
    graph.add_edge(source_id, target_id);
    return true; // successful load
  } 
 

  const bool success = qi::parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::omit[qi::ulong_[phoenix::ref(source_id) = qi::_1] ]>> qi::omit[+qi::space] >> //-qi::char_(',') >>
       +qi::alnum % qi::omit[+qi::space]
       ),
        str_vec);
     
      
  if(!success) return false;
  vertex_data_type v_data(str_vec);
//  v_data.count=v_data.str_vec.size();
//  graph.add_edge(source_id, target_id, edge_data(obs, role));
  graph.add_vertex(source_id, v_data); 
  return true; // successful load
} // end of graph_loader


/*
 * The factorized page rank update function extends ivertex_program
 * specifying the:
 *
 *   1) graph_type
 *   2) gather_type: double (returned by the gather function). Note
 *      that the gather type is not strictly needed here since it is
 *      assumed to be the same as the vertex_data_type unless
 *      otherwise specified
 *
 * In addition ivertex program also takes a message type which is
 * assumed to be empty. Since we do not need messages no message type
 * is provided.
 *
 * pagerank also extends graphlab::IS_POD_TYPE (is plain old data type)
 * which tells graphlab that the pagerank program can be serialized
 * (converted to a byte stream) by directly reading its in memory
 * representation.  If a vertex program does not exted
 * graphlab::IS_POD_TYPE it must implement load and save functions.
 */
//graphlab::ivertex_program< Graph, GatherType, MessageType >

class pagerank :

  public graphlab::ivertex_program<graph_type, int> {

  double last_change;
public:

  /**
   * Gather only in edges.
   */
/*
  void init(icontext_type& context, const vertex_type& vertex,const message_type& message) {
    
    //if(!vertex.data().is_doc) 
      context.signal(vertex);
  } 
*/
  edge_dir_type gather_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    if(vertex.data().is_doc)
      return graphlab::NO_EDGES;
    else
      return graphlab::ALL_EDGES;
  } // end of Gather edges


  /* Gather the weighted rank of the adjacent page   */
  int gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    if(vertex.data().is_doc)
      return 0;
    else{
      int count=0;
      vertex_type doc=get_other_vertex(edge,vertex);
      for(int i=0;i<doc.data().str_vec.size();i++){
        if(doc.data().str_vec[i]=="google"){
          count ++;
        }
      }
      return count;
    }
      
  }

  /* Use the total rank of adjacent pages to update this page */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& total) {
    vertex.data().count=total;
  }

  /* The scatter edges depend on whether the pagerank has converged */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
    
  }

  /* The scatter function just signal adjacent pages */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    
  }

  void save(graphlab::oarchive& oarc) const {
    // If we are using iterations as a counter then we do not need to
    // move the last change in the vertex program along with the
    // vertex data.
    if (ITERATIONS == 0) oarc << last_change;
  }

  void load(graphlab::iarchive& iarc) {
    if (ITERATIONS == 0) iarc >> last_change;
  }

}; // end of factorized_pagerank update functor




int map_rank(const graph_type::vertex_type& v) { return v.data().count; }

edge_data_type signal_target(graphlab::omni_engine<pagerank>::icontext_type& context,
                                    graph_type::edge_type edge) {
  context.signal(edge.target());
  return graphlab::empty();
}



int pagerank_sum(graph_type::vertex_type v) {
  return v.data().count;
}

int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_dir;
  std::string format = "adj";
  std::string exec_type = "synchronous";
  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("engine", exec_type,
                       "The engine type synchronous or asynchronous");
  clopts.attach_option("tol", TOLERANCE,
                       "The permissible change at convergence.");
  clopts.attach_option("format", format,
                       "The graph file format");
  size_t powerlaw = 0;
  clopts.attach_option("powerlaw", powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
  clopts.attach_option("iterations", ITERATIONS,
                       "If set, will force the use of the synchronous engine"
                       "overriding any engine option set by the --engine parameter. "
                       "Runs complete (non-dynamic) PageRank for a fixed "
                       "number of iterations. Also overrides the iterations "
                       "option in the engine");
  clopts.attach_option("use_delta", USE_DELTA_CACHE,
                       "Use the delta cache to reduce time in gather.");
  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }


  // Enable gather caching in the engine
  clopts.get_engine_args().set_option("use_cache", USE_DELTA_CACHE);

  if (ITERATIONS) {
    // make sure this is the synchronous engine
    dc.cout() << "--iterations set. Forcing Synchronous engine, and running "
              << "for " << ITERATIONS << " iterations." << std::endl;
    //clopts.get_engine_args().set_option("type", "synchronous");
    clopts.get_engine_args().set_option("max_iterations", ITERATIONS);
    clopts.get_engine_args().set_option("sched_allv", true);
  }

  // Build the graph ----------------------------------------------------------
  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; 
  graph_type graph(dc, clopts);
  if(powerlaw > 0) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, 2.1, 100000000);
  }
  else if (graph_dir.length() > 0) { // Load the graph from a file

    dc.cout() << "Loading graph in format: "<< format << std::endl;
    graph.load(graph_dir, graph_loader); 
//    graph.load_format(graph_dir, format);
  }
  else {
    dc.cout() << "graph or powerlaw option must be specified" << std::endl;
    clopts.print_description();
    return 0;
  }
  dc.cout() << "Loading graph. Finished in " 
    << timer.current_time() << std::endl;

  size_t bytes_sent=dc.bytes_sent();
  size_t calls_sent=dc.calls_sent();
  size_t network_bytes_sent=dc.network_bytes_sent();
  size_t bytes_received=dc.bytes_received();
  size_t calls_received=dc.calls_received();
  dc.cout() << "load_Bytes_Sent: "     << bytes_sent        << std::endl;
  dc.cout() << "load_Calls_Sent: "     << calls_sent        << std::endl;
  dc.cout() << "load_Network_Sent: "   << network_bytes_sent<< std::endl;
  dc.cout() << "load_Bytes_Received: " << bytes_received    << std::endl;
  dc.cout() << "load_Calls_Received: " << calls_received    << std::endl;

  // must call finalize before querying the graph
  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " 
    << timer.current_time() << std::endl;

  bytes_sent=dc.bytes_sent()-bytes_sent;
  calls_sent=dc.calls_sent()-calls_sent;
  network_bytes_sent=dc.network_bytes_sent()-network_bytes_sent;
  bytes_received=dc.bytes_received()-bytes_received;
  calls_received=dc.calls_received()-calls_received;
  dc.cout() << "finalize_Bytes_Sent: "     << bytes_sent        << std::endl;
  dc.cout() << "finalize_Calls_Sent: "     << calls_sent        << std::endl;
  dc.cout() << "finalize_Network_Sent: "   << network_bytes_sent<< std::endl;
  dc.cout() << "finalize_Bytes_Received: " << bytes_received    << std::endl;
  dc.cout() << "finalize_Calls_Received: " << calls_received    << std::endl;


  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;

  // Initialize the vertex data
  //graph.transform_vertices(init_vertex);

  // Running The Engine -------------------------------------------------------
  graphlab::omni_engine<pagerank> engine(dc, graph, exec_type, clopts);

  
  engine.map_reduce_edges<graphlab::empty>(signal_target);

  bytes_sent=dc.bytes_sent()-bytes_sent;
  calls_sent=dc.calls_sent()-calls_sent;
  network_bytes_sent=dc.network_bytes_sent()-network_bytes_sent;
  bytes_received=dc.bytes_received()-bytes_received;
  calls_received=dc.calls_received()-calls_received;
  dc.cout() << "before_start_Bytes_Sent: "     << bytes_sent        << std::endl;
  dc.cout() << "before_start_Calls_Sent: "     << calls_sent        << std::endl;
  dc.cout() << "before_start_Network_Sent: "   << network_bytes_sent<< std::endl;
  dc.cout() << "before_start_Bytes_Received: " << bytes_received    << std::endl;
  dc.cout() << "before_start_Calls_Received: " << calls_received    << std::endl;

  //engine.signal_all();
  timer.start();
  engine.start();
  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
            << std::endl
            << "Final Runtime (seconds):   " << runtime 
            << std::endl
            << "Updates executed: " << engine.num_updates() << std::endl
            << "Update Rate (updates/second): " 
            << engine.num_updates() / runtime << std::endl;

  const int total_rank = graph.map_reduce_vertices<int>(map_rank);
  std::cout << "Total rank: " << total_rank << std::endl;

  // Save the final graph -----------------------------------------------------

  bytes_sent=dc.bytes_sent()-bytes_sent;
  calls_sent=dc.calls_sent()-calls_sent;
  network_bytes_sent=dc.network_bytes_sent()-network_bytes_sent;
  bytes_received=dc.bytes_received()-bytes_received;
  calls_received=dc.calls_received()-calls_received;
  dc.cout() << "compute_Bytes_Sent: "     << bytes_sent        << std::endl;
  dc.cout() << "compute_Calls_Sent: "     << calls_sent        << std::endl;
  dc.cout() << "compute_Network_Sent: "   << network_bytes_sent<< std::endl;
  dc.cout() << "compute_Bytes_Received: " << bytes_received    << std::endl;
  dc.cout() << "compute_Calls_Received: " << calls_received    << std::endl;

  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


// We render this entire program in the documentation


