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
#include <set>
#include <iterator>

#include <graphlab.hpp>

/**
 * \brief The type used to measure distances in the graph.
 */
typedef float distance_type;

/**
 * \brief The type used to collect paths in the graph.
 */
typedef std::set<double> path_type;

/**
 * \brief The current distance of the vertex.
 */
struct vertex_data {
    distance_type dist;
    distance_type parent_node;
    path_type paths;

    vertex_data(distance_type dist = std::numeric_limits<distance_type>::max()) :
            dist(dist), parent_node(dist) {
    }

    void save(graphlab::oarchive& oarc) const {
        oarc << paths << dist << parent_node;
    }

    void load(graphlab::iarchive& iarc) {
        iarc >> paths >> dist >> parent_node;
    }
}; // end of vertex data



/**
 * \brief The distance associated with the edge.
 */
struct edge_data : graphlab::IS_POD_TYPE {
  distance_type dist;
  edge_data(distance_type dist = 1) : dist(dist) { }
}; // end of edge data


/**
 * \brief The graph type encodes the distances between vertices and
 * edges
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/**
 * \brief Get the other vertex in the edge.
 */
inline graph_type::vertex_type
get_other_vertex(const graph_type::edge_type& edge,
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}

/**
 * \brief Get the parent node of the vertex ID.
 */
inline graph_type::vertex_type get_parent_node(const graph_type::edge_type& edge,
        const graph_type::vertex_type& vertex) {
    return vertex.id() == edge.source().id() ? edge.source() : edge.target();
}

/**
 * \brief Use directed or undireced edges.
 */
bool DIRECTED_SSSP = false;


/**
 * \brief This class is used as the gather type.
 */
struct min_distance_type {
  distance_type dist;
  distance_type parent_node;
  path_type alternate_paths;
  bool isAlternatePath;
  min_distance_type(distance_type dist = 
                    std::numeric_limits<distance_type>::max(),
                    distance_type parent_node =
                    std::numeric_limits<distance_type>::max()) : dist(dist), parent_node(parent_node), isAlternatePath(false) { }
  min_distance_type& operator+=(const min_distance_type& other) {
    if (dist == other.dist) {
            alternate_paths.insert(other.parent_node);
            isAlternatePath = true;
        } else {
            isAlternatePath = false;
        }
    dist = std::min(dist, other.dist);
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
        oarc << alternate_paths << dist << parent_node << isAlternatePath;
    }

    void load(graphlab::iarchive& iarc) {
        iarc >> alternate_paths >> dist >> parent_node >> isAlternatePath;
    }
};


/**
 * \brief The single source shortest path vertex program.
 */
class sssp :
  public graphlab::ivertex_program<graph_type, 
                                   graphlab::empty,
                                   min_distance_type>
  {
  distance_type min_dist;
  bool changed;
  distance_type parent_node;
  path_type alternate_paths;
  bool isAlternatePath;
public:


  void init(icontext_type& context, const vertex_type& vertex,
            const min_distance_type& msg) {
    min_dist = msg.dist;
    parent_node = msg.parent_node;
    isAlternatePath = msg.isAlternatePath;
    alternate_paths = msg.alternate_paths;
  } 

  /**
   * \brief We use the messaging model to compute the SSSP update
   */
  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const { 
    return graphlab::NO_EDGES;
  }; // end of gather_edges 


  // /** 
  //  * \brief Collect the distance to the neighbor
  //  */
  // min_distance_type gather(icontext_type& context, const vertex_type& vertex, 
  //                          edge_type& edge) const {
  //   return min_distance_type(edge.data() + 
  //                            get_other_vertex(edge, vertex).data());
  // } // end of gather function


  /**
   * \brief If the distance is smaller then update
   */
  void apply(icontext_type& context, vertex_type& vertex,
             const graphlab::empty& empty) {
    changed = false;
    if(vertex.data().dist > min_dist) {
      changed = true;
      vertex.data().dist = min_dist;
    }
    if (isAlternatePath == false) {
            vertex.data().paths.insert(parent_node);
        } else {
            vertex.data().paths = alternate_paths;
            vertex.data().paths.insert(parent_node);
        }
  }

  /**
   * \brief Determine if SSSP should run on all edges or just in edges
   */
  edge_dir_type scatter_edges(icontext_type& context, 
                             const vertex_type& vertex) const {
    if(changed)
      return DIRECTED_SSSP? graphlab::OUT_EDGES : graphlab::ALL_EDGES; 
    else return graphlab::NO_EDGES;
  }; // end of scatter_edges

  /**
   * \brief The scatter function just signal adjacent pages 
   */
  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    const vertex_type other = get_other_vertex(edge, vertex);
    const vertex_type path = get_parent_node(edge, vertex);
    distance_type newd = vertex.data().dist + edge.data().dist;
    if (other.data().dist > newd) {
      const min_distance_type msg(newd,path.id());
      context.signal(other, msg);
    }
  } // end of scatter

  void save(graphlab::oarchive& oarc) const {
        oarc << alternate_paths << parent_node << min_dist << isAlternatePath
                << changed;
    }

    void load(graphlab::iarchive& iarc) {
        iarc >> alternate_paths >> parent_node >> min_dist >> isAlternatePath
                >> changed;
    }

}; // end of shortest path vertex program




/**
 * \brief We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", pagerank_writer()) to save the graph.
 */
struct shortest_path_writer {
  std::string save_vertex(const graph_type::vertex_type& vtx) {
    std::stringstream strm;
    if (vtx.data().dist == 0) {
            strm << vtx.id() << "\t" << vtx.data().dist << "\t" << "'"
                    << vtx.id() << "'" << std::endl;
        } else {
            strm << vtx.id() << "\t" << vtx.data().dist << "\t";
            std::set<double>::iterator it;

            for (it = vtx.data().paths.begin(); it != vtx.data().paths.end();
                    it++) {
                strm << "'" << *it << "'";
            }
            strm << std::endl;
        }
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of shortest_path_writer



struct max_deg_vertex_reducer: public graphlab::IS_POD_TYPE {
  size_t degree;
  graphlab::vertex_id_type vid;
  max_deg_vertex_reducer& operator+=(const max_deg_vertex_reducer& other) {
    if (degree < other.degree) {
      (*this) = other;
    }
    return (*this);
  }
};

max_deg_vertex_reducer find_max_deg_vertex(const graph_type::vertex_type vtx) {
  max_deg_vertex_reducer red;
  red.degree = vtx.num_in_edges() + vtx.num_out_edges();
  red.vid = vtx.id();
  return red;
}

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
  bool max_degree_source = false;
  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "graph format");
  clopts.attach_option("source", sources,
                       "The source vertices");
  clopts.attach_option("max_degree_source", max_degree_source,
                       "Add the vertex with maximum degree as a source");

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


  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  if(powerlaw > 0) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, 2, 100000000);
  } else if (graph_dir.length() > 0) { // Load the graph from a file
    dc.cout() << "Loading graph in format: "<< format << std::endl;
    graph.load_format(graph_dir, format);
  } else {
    dc.cout() << "graph or powerlaw option must be specified" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  // must call finalize before querying the graph
  graph.finalize();
  dc.cout() << "#vertices:  " << graph.num_vertices() << std::endl
            << "#edges:     " << graph.num_edges() << std::endl;



  if(sources.empty()) {
    if (max_degree_source == false) {
      dc.cout()
        << "No source vertex provided. Adding vertex 0 as source" 
        << std::endl;
      sources.push_back(0);
    }
  }

  if (max_degree_source) {
    max_deg_vertex_reducer v = graph.map_reduce_vertices<max_deg_vertex_reducer>(find_max_deg_vertex);
    dc.cout()
      << "No source vertex provided.  Using highest degree vertex " << v.vid << " as source."
      << std::endl;
    sources.push_back(v.vid);
  }



  // Running The Engine -------------------------------------------------------
  graphlab::omni_engine<sssp> engine(dc, graph, exec_type, clopts);


  
  // Signal all the vertices in the source set
  for(size_t i = 0; i < sources.size(); ++i) {
    engine.signal(sources[i], min_distance_type(0));
  }

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



