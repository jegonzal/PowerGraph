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
#include <string>
#include <fstream>


#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/unordered_set.hpp>

#include <graphlab.hpp>

#include <graphlab/util/stl_util.hpp>


#include <graphlab/macros_def.hpp>


/**
 * \brief The type used to measure distances in the graph.
 */
typedef float distance_type;

/**
 * \brief The current distance of the vertex.
 */
typedef distance_type vertex_data;


/**
 * \brief The distance associated with the edge.
 */
typedef distance_type edge_data;


/**
 * \brief The graph type encodes the distances between vertices and
 * edges
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



/**
 * \brief Initialize all vertices to have distance infinity except for
 * the source vertices which should have distance 0.
 *
 */
void initialize_vertices(graph_type::vertex_type& vertex, 
                         boost::unordered_set<graphlab::vertex_id_type>* sources_ptr) {
  if(sources_ptr->count(vertex.id()) > 0) vertex.data() = 0;
  else vertex.data() = std::numeric_limits<distance_type>::max();
}


/**
 * \brief The graph loader is used by graph.load to parse lines of the
 * text data file.
 *
 * We use the relativley fast boost::spirit parser to parse each line.
 */
bool graph_loader(graph_type& graph, const std::string& fname,
                  const std::string& line) {
  ASSERT_FALSE(line.empty());
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  graphlab::vertex_id_type source_id(-1), target_id(-1);
  float weight = 1;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(source_id) = qi::_1] >> -qi::char_(',') >>
      qi::ulong_[phoenix::ref(target_id) = qi::_1] >> 
      -(-qi::char_(',') >> qi::float_[phoenix::ref(weight) = qi::_1])
      )
     ,
     //  End grammar
     ascii::space); 
  if(!success) return false;
  if(source_id == target_id) {
    logstream(LOG_ERROR) 
      << "Self edge to vertex " << source_id << " is not permitted." << std::endl;
  }
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, weight);
  return true; // successful load
}; // end of graph loader



/**
 * \brief Get the other vertex in the edge.
 */
inline graph_type::vertex_type
get_other_vertex(const graph_type::edge_type& edge,
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}


/**
 * \brief Use directed or undireced edges.
 */
bool DIRECTED_SSSP = false;


/**
 * \brief This class is used as the gather type.
 */
struct min_distance_type : graphlab::IS_POD_TYPE {
  distance_type dist;
  min_distance_type(distance_type dist = 
                    std::numeric_limits<distance_type>::max()) : dist(dist) { }
  min_distance_type& operator+=(const min_distance_type& other) {
    dist = std::min(dist, other.dist);
    return *this;
  }
};


/**
 * \brief The single source shortest path vertex program.
 */
class sssp :
  public graphlab::ivertex_program<graph_type, min_distance_type>,
  public graphlab::IS_POD_TYPE {
  bool changed;
public:

  /**
   * \brief Determine if SSSP should run on all edges or just in edges
   */
  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const { 
    return DIRECTED_SSSP? graphlab::IN_EDGES : graphlab::ALL_EDGES; 
  }; // end of gather_edges 


  /** 
   * \brief Collect the distance to the neighbor
   */
  min_distance_type gather(icontext_type& context, const vertex_type& vertex, 
                           edge_type& edge) const {
    return min_distance_type(edge.data() + 
                             get_other_vertex(edge, vertex).data());
  } // end of gather function


  /**
   * \brief If the distance is smaller then update
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const min_distance_type& total) {
    changed = false;
    if(vertex.data() > total.dist) {
      vertex.data() = total.dist;
      changed = true;
    }
  }

  /**
   * \brief Determine if SSSP should run on all edges or just in edges
   */
  edge_dir_type scatter_edges(icontext_type& context, 
                             const vertex_type& vertex) const {
    if(changed) return DIRECTED_SSSP? graphlab::OUT_EDGES : graphlab::ALL_EDGES; 
    else return graphlab::NO_EDGES;
  }; // end of scatter_edges

  /**
   * \brief The scatter function just signal adjacent pages 
   */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    context.signal(get_other_vertex(edge, vertex));
  } // end of scatter

}; // end of shortest path vertex program




/**
 * \brief We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", pagerank_writer()) to save the graph.
 */
struct shortest_path_writer {
  std::string save_vertex(const graph_type::vertex_type& vtx) {
    std::stringstream strm;
    strm << vtx.id() << "\t" << vtx.data() << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of shortest_path_writer





int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options 
    clopts("Single Source Shortest Path Algorithm.");
  std::string graph_dir;
  std::string format = "adj";
  std::string exec_type = "synchronous";
  size_t powerlaw = 0;
  std::vector<graphlab::vertex_id_type> sources;
  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");

  clopts.attach_option("source", sources,
                       "The source vertices");
  clopts.add_positional("source");

  clopts.attach_option("directed", DIRECTED_SSSP,
                       "Treat edges as directed.");

  clopts.attach_option("engine", exec_type, 
                       "The engine type synchronous or asynchronous");
 
  
  clopts.attach_option("powerlaw", powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  if(sources.empty()) {
    dc.cout()
      << "No source vertex provided.  Using vertex 0 as source."
      << std::endl;
    sources.push_back(0);
  }



  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  if(powerlaw > 0) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, 2.1, 100000000);
  } else { // Load the graph from a file
    dc.cout() << "Loading graph in format: "<< format << std::endl;
    graph.load(graph_dir, graph_loader);
  }
  // must call finalize before querying the graph
  graph.finalize();
  dc.cout() << "#vertices:  " << graph.num_vertices() << std::endl
            << "#edges:     " << graph.num_edges() << std::endl;


  // Initialize the vertex data
  boost::unordered_set<graphlab::vertex_id_type> source_set;
  for(size_t i = 0; i < sources.size(); ++i) 
    source_set.insert(sources[i]);
  graph.transform_vertices(boost::bind(initialize_vertices, _1, &source_set));

  // Running The Engine -------------------------------------------------------
  graphlab::omni_engine<sssp> engine(dc, graph, exec_type, clopts);
  
  


  engine.start();
  const float runtime = engine.elapsed_seconds();
  dc.cout() << "Finished Running engine in " << runtime
            << " seconds." << std::endl;


  // Save the final graph -----------------------------------------------------
  if (saveprefix != "") {
    graph.save(saveprefix, shortest_path_writer(),
               false,    // do not gzip
               true,     // save vertices
               false);   // do not save edges
  }

  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


// We render this entire program in the documentation



