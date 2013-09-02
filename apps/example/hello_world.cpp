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


#include <graphlab.hpp>

// The vertex data is just the pagerank value (a float)
typedef float vertex_data_type;

// There is no edge data in the pagerank application
typedef float edge_data_type;
typedef float message_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

/*
 * A simple function used by graph.transform_vertices(init_vertex);
 * to initialize the vertex data.
 */
void init_vertex(graph_type::vertex_type& vertex) { vertex.data() = vertex.id(); }

struct min_combiner {
    graphlab::vertex_id_type v;
    min_combiner& operator+=(const min_combiner& other) { 
        v = std::min(v, other.v);  
        return *this; 
    }
};

class concomp :
    public graphlab::ivertex_program<graph_type, float>,
    public graphlab::IS_POD_TYPE {
    bool changed;
public:
    float gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
       float edge_data = edge.source().data();
       float vertex_data = vertex.data();

       std::cout << "current vertex id: " << vertex.id() << " data: " << vertex.data() << "\n";


       std::cout << "\tedge vertex id: " << edge.source().id() << " data: " << edge_data << "\n";
       if (edge_data < vertex_data) {
           std::cout << "returning edge data: " << edge_data << "\n";
           return edge_data;
       } else {
           std::cout << "returning vertex data: " << vertex_data << "\n";
           return vertex_data;
       }
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& smallest) {
        std::cout << "vertex id: " << vertex.id() << " data: " << vertex.data() << "\n";
        std::cout << "smallest: " << smallest << "\n";

        if (smallest < vertex.data()) {
            vertex.data() = smallest;
            changed = true;
        } else {
            changed = false;
        }
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
        if (changed) {
            return graphlab::ALL_EDGES;
        } else {
            return graphlab::NO_EDGES;
        }
    }

    /* The scatter function just signal adjacent pages */
    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
        context.signal(edge.target());
    }
};

/* We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", concomp_writer()) to save the graph.
 */
struct concomp_writer {
    std::string save_vertex(graph_type::vertex_type v) {
        std::stringstream strm;
        strm << v.id() << "\t" << v.data() << "\n";
        return strm.str();
    }
    std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of concomp writer
            

int main(int argc, char** argv) {
   // Initialize control plain using mpi
   graphlab::mpi_tools::init(argc, argv);
   graphlab::distributed_control dc;
   global_logger().set_log_level(LOG_INFO);
    
   // Parse command line options -----------------------------------------------
   graphlab::command_line_options clopts("PageRank algorithm.");
   std::string graph_dir;
   std::string format = "snap";
   clopts.attach_option("graph", graph_dir, "The graph file. Required ");
   clopts.add_positional("graph");
   clopts.attach_option("format", format, "The graph file format");
   if(!clopts.parse(argc, argv)) {
       dc.cout() << "Error in parsing command line arguments." << std::endl;
       return EXIT_FAILURE;
   }
   if (graph_dir == "") {
       dc.cout() << "Graph not specified. Cannot continue";
       return EXIT_FAILURE;
   }
   
   // Build the graph ----------------------------------------------------------
   graph_type graph(dc, clopts);
   dc.cout() << "Loading graph in format: "<< format << std::endl;
   graph.load_format(graph_dir, format);
   // must call finalize before querying the graph
   graph.finalize();

   graph.transform_vertices(init_vertex);

   dc.cout() << "#vertices: " << graph.num_vertices() << " #edges:" << graph.num_edges() << std::endl;

   graphlab::omni_engine<concomp> engine(dc, graph, "synchronous", clopts);
   engine.signal_all();

   engine.start();

   const float runtime = engine.elapsed_seconds();
   dc.cout() << "Finished Running engine in " << runtime << " seconds." << std::endl;

   graph.save("output" + graph_dir + ".txt", concomp_writer(),
               false,    // do not gzip
               true,     // save vertices
               false);   // do not save edges
}
