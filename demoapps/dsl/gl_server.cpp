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
#include <graphlab/engine/gl3engine.hpp>

#include <dlfcn.h>

using namespace graphlab;

#include "graph_typedefs.gen"

#define GET_EDGE_DATA 100
#define GET_VERTEX_DATA 101

typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;
typedef gl3engine<graph_type> engine_type;


//here we keep one global copy of the user defined function pointers
void (*user_program)(user_funs*) = NULL;
void (*vertex_reduce)(vertex_data_type&,const vertex_data_type&) = NULL;
void (*edge_reduce)(edge_data_type&,const edge_data_type&) = NULL;

//here are definitions of functions passed to dynamically linked user code
graph_type::vertex_data_type get_vertex_data(graph_type::vertex_type* v) {
    return v->data();
}
void set_vertex_data(graph_type::vertex_data_type& d, graph_type::vertex_type* v) {
    v->data() = d;
}

// graph_type::edge_data_type edge_get_map(graph_type::edge_type& e) {
//     return e.data();
// }
graph_type::vertex_data_type vertex_get_map(const graph_type::vertex_type& v) {
    return v.data();
}

// void edge_get_reduce(edge_data_type& ev, const edge_data_type ed) {
//     edge_reduce(ev,ed);
// }

void vertex_get_reduce(vertex_data_type& ev, const vertex_data_type ed) {
    vertex_reduce(ev,ed);
}

//wonder why this isn't working. hm...
// void edge_get_reduce(std::vector<graph_type::edge_data_type>& ev, const graph_type::edge_data_type ed) {
//     ev.push_back(ed);
// }

// void vertex_get_reduce(std::vector<graph_type::vertex_data_type>& vv, const graph_type::edge_data_type vd) {
//     vv.push_back(vd);
// }

// std::vector<edge_data_type> get_neighboring_edges() {
//     return current_context->map_reduce<std::vector<edge_data_type> >(GET_EDGE_DATA, ALL_EDGES);
// }
// std::vector<vertex_data_type> get_neighboring_vertices() {
//     return current_context->map_reduce<std::vector<vertex_data_type> >(GET_VERTEX_DATA, ALL_EDGES);
// }

vertex_data_type reduce_neighbors(edge_dir_type d, engine_type::context_type* ctx) {
    return ctx->map_reduce<vertex_data_type>(GET_VERTEX_DATA, d);
}

void signal_neighbors(edge_dir_type d, engine_type::context_type* ctx) {
    ctx->broadcast_signal(d);
}

// edge_data_type reduce_edges() {
//     return current_context->map_reduce<edge_data_type>(GET_EDGE_DATA, ALL_EDGES);
// }

void server_program(engine_type::context_type& context,
		      graph_type::vertex_type& vertex,
		      const engine_type::message_type& unused) {



        //capture functions for user to call
    user_funs f;
    f.vtx = &vertex;
    f.ctx = &context;
    f._get_vertex_data = (vertex_data_type (*)(void*))get_vertex_data;
    f._set_vertex_data = (void (*)(vertex_data_type&,void*))set_vertex_data;
    f._reduce_neighbors = (vertex_data_type (*)(edge_dir_type,void*))reduce_neighbors;
    f._signal_neighbors = (void (*)(edge_dir_type,void*))signal_neighbors;
    //f.reduce_edges = reduce_edges;

    //call dynamically linked user function
    user_program(&f);
    
    // float prev = vertex.data();
    // // map reduce over neighbors
    // vertex.data() = 0.15 + 0.85 *
    // 	context.map_reduce<float>(PAGERANK_MAP_REDUCE, IN_EDGES);

    // float last_change = std::fabs((vertex.data()- prev) / vertex.num_out_edges());
    // if (last_change > TOLERANCE) {
    // 	// signals out neighbors if I change substantially
    // 	context.broadcast_signal(OUT_EDGES);
    // }
}



void init_vertex(graph_type::vertex_type& vertex) { vertex.data() = 1; }

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
                       "The graph file. Required ");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "The graph file format");
  clopts.attach_option("engine", exec_type, 
                       "The engine type synchronous or asynchronous");
  //clopts.attach_option("tol", TOLERANCE,
  //                     "The permissible change at convergence.");
  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");
  //set single threaded
  //clopts.set_ncpus(1);
  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  if (graph_dir == "") {
    dc.cout() << "Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }


  // Dynamically link the functions
  dc.cout() << "Starting dynamic linking" << std::endl;
  void* handle = dlopen("/home/emullen/graphlab/graphlab2.2/graphlabapi/debug/demoapps/dsl/libgen_impl.so", RTLD_LAZY);
  if (handle == NULL) {
      dc.cout() << dlerror() << std::endl;
      assert(handle != NULL);
  }

  // pagerank_map = (float (*)(const graphlab::distributed_graph<float, graphlab::empty>::vertex_type&))dlsym(pagerank_handle, "pagerank_map");
  // assert(pagerank_map != NULL);
  // pagerank_combine = (void (*)(float&, const float&))dlsym(pagerank_handle, "pagerank_combine");
  // assert(pagerank_combine != NULL);
  user_program = (void (*)(user_funs*))dlsym(handle, "user_program");
  assert(user_program != NULL);
  vertex_reduce = (void (*)(vertex_data_type&,const vertex_data_type&))dlsym(handle, "vertex_reduce");
  assert(vertex_reduce != NULL);
  // edge_reduce = (void (*)(edge_data_type&,const edge_data_type&))dlsym(handle, "edge_reduce");
  // assert(edge_reduce != NULL);
  dc.cout() << "Finished dynamic linking" << std::endl;

  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  dc.cout() << "Loading graph in format: "<< format << std::endl;
  graph.load_format(graph_dir, format);
  // must call finalize before querying the graph
  graph.finalize();
  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;

  // Initialize the vertex data
  graph.transform_vertices(init_vertex);

    // Running The Engine -------------------------------------------------------
  engine_type engine(dc, graph, clopts);

  // register the map reduce operation before usage
  // Each task registration must have a distinct ID ranging fro 0 to 223
  //engine.register_map_reduce(PAGERANK_MAP_REDUCE,
  //                           pagerank_map,
  //                           pagerank_combine);
  // engine.register_map_reduce(GET_EDGE_DATA,
  // 			     edge_get_map,
  // 			     edge_get_reduce);
  engine.register_map_reduce(GET_VERTEX_DATA,
			     vertex_get_map,
			     vertex_get_reduce);

  engine.set_vertex_program(server_program);
  engine.signal_all();
  engine.wait();
  
  // Save the final graph -----------------------------------------------------
  // if (saveprefix != "") {
  //   graph.save(saveprefix, pagerank_writer(),
  //              false,    // do not gzip
  //              true,     // save vertices
  //              false);   // do not save edges
  // }

  dlclose(handle);

  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main
