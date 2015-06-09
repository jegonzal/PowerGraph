/*  
 * Copyright (c) 2013 Shanghai Jiao Tong University. 
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
 *      http://ipads.se.sjtu.edu.cn/projects/powerlyra.html
 *
 *
 * 2014.02  implement word search application for testing bipartite-aware partitiong
 *              with affinity
 *
 */



#include <vector>
#include <string>
#include <fstream>

#include <graphlab.hpp>

bool USE_DELTA_CACHE = false;

struct vertex_t {
  bool is_doc;
  int count;
  std::vector<std::string> str_vec; // only for doc

  vertex_t(std::vector<std::string>& vec) : count(0)
    { str_vec = vec; is_doc = true; }

  vertex_t() : count(0)
    { is_doc = false; }

  void save(graphlab::oarchive& arc) const { 
    arc << is_doc << count << str_vec ;         
  }

  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) { 
    arc >> is_doc >> count >> str_vec ;
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
  graph.add_vertex(source_id, v_data); 
  return true; // successful load
} // end of graph_loader


class wsearch :

  public graphlab::ivertex_program<graph_type, int> {

public:

  edge_dir_type gather_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    if(vertex.data().is_doc) return graphlab::NO_EDGES;
    else return graphlab::ALL_EDGES;
  } // end of Gather edges


  int gather(icontext_type& context,
              const vertex_type& vertex, edge_type& edge) const {
    if(vertex.data().is_doc) return 0;
    else {
      int count = 0;
      vertex_type doc = get_other_vertex(edge, vertex);
      for(int i = 0; i < doc.data().str_vec.size(); i++)
        if(doc.data().str_vec[i] == "google") count ++;
      return count;
    }
  }

  void apply(icontext_type& context, vertex_type& vertex,
              const gather_type& total) {
    vertex.data().count = total;
  }

  edge_dir_type scatter_edges(icontext_type& context,
                                const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }

  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const { }

  void save(graphlab::oarchive& oarc) const { }

  void load(graphlab::iarchive& iarc) { }

};

int map_count(const graph_type::vertex_type& v) { return v.data().count; }

edge_data_type 
signal_target(graphlab::omni_engine<wsearch>::icontext_type& context,
               graph_type::edge_type edge) {
  context.signal(edge.target());
  return graphlab::empty();
}

int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Word serach application for data affinity.");
  std::string graph_dir;
  std::string exec_type = "synchronous";
  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("engine", exec_type,
                       "The engine type synchronous or asynchronous");
  clopts.attach_option("use_delta", USE_DELTA_CACHE,
                       "Use the delta cache to reduce time in gather.");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }


  // Enable gather caching in the engine
  clopts.get_engine_args().set_option("use_cache", USE_DELTA_CACHE);

  // Build the graph ----------------------------------------------------------
  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; 
  graph_type graph(dc, clopts);
  if (graph_dir.length() > 0) { // Load the graph from a file
    graph.load(graph_dir, graph_loader); 
  } else {
    clopts.print_description();
    return 0;
  }
  dc.cout() << "Loading graph. Finished in " 
    << timer.current_time() << std::endl;

  size_t bytes_sent = dc.bytes_sent();
  size_t calls_sent = dc.calls_sent();
  size_t network_bytes_sent = dc.network_bytes_sent();
  size_t bytes_received = dc.bytes_received();
  size_t calls_received = dc.calls_received();
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

  bytes_sent = dc.bytes_sent() - bytes_sent;
  calls_sent = dc.calls_sent() - calls_sent;
  network_bytes_sent = dc.network_bytes_sent() - network_bytes_sent;
  bytes_received = dc.bytes_received() - bytes_received;
  calls_received = dc.calls_received() - calls_received;
  dc.cout() << "finalize_Bytes_Sent: "     << bytes_sent        << std::endl;
  dc.cout() << "finalize_Calls_Sent: "     << calls_sent        << std::endl;
  dc.cout() << "finalize_Network_Sent: "   << network_bytes_sent<< std::endl;
  dc.cout() << "finalize_Bytes_Received: " << bytes_received    << std::endl;
  dc.cout() << "finalize_Calls_Received: " << calls_received    << std::endl;


  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;


  // Running The Engine -------------------------------------------------------
  graphlab::omni_engine<wsearch> engine(dc, graph, exec_type, clopts);

  // Initialize the vertex data  
  engine.map_reduce_edges<graphlab::empty>(signal_target);

  bytes_sent = dc.bytes_sent() - bytes_sent;
  calls_sent = dc.calls_sent() - calls_sent;
  network_bytes_sent = dc.network_bytes_sent() - network_bytes_sent;
  bytes_received = dc.bytes_received() - bytes_received;
  calls_received = dc.calls_received() - calls_received;
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

  const int total_count = graph.map_reduce_vertices<int>(map_count);
  std::cout << "Total count: " << total_count << std::endl;

  bytes_sent = dc.bytes_sent() - bytes_sent;
  calls_sent = dc.calls_sent() - calls_sent;
  network_bytes_sent = dc.network_bytes_sent() - network_bytes_sent;
  bytes_received = dc.bytes_received() - bytes_received;
  calls_received = dc.calls_received() - calls_received;
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


