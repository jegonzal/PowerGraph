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


/**
   This demo provides a synthetic application which intentionally uses
   most of the GraphLab V2 features.  Most GraphLab applications will
   not use all of these features.

   Picture of grid model:

      x----x----x
      |    |    |
      x----o----o
      |    |    |
      o----x----x

   Given a DIMxDIM undirected grid graph where each vertex is assigned
   a random color, either black or red.  Each vertex then makes a
   local random decision:

   - become red with probabilty proportionate to the number of red
   neighbors

   - become black with probabilty proportionate to the number of black
   neighbors

   Clearly, the two stable outcomes are where all vertices are black,
   or where all vertices are red. We are interested in knowing how
   many flips each vertex took on average.

   Also, since GraphLab only has directed edges, we will build the
   graph by duplicating every edge in both directions.

*/


// standard C++ headers
#include <iostream>
#include <algorithm>
// includes the entire graphlab framework
#include <distributed_graphlab.hpp>

// A simple enum describing the two possible colors
enum color_type {RED, BLACK};

/**
   First we will design the graph data. For each vertex, we will need
   to know its current color, and a counter to count the number of
   flips it took.

   GraphLab provides facilities to directly save/load graphs from
   disk. However, to do so, GraphLab must be able to understand your
   datastructures.  If you are not interested in saving/loading
   graphs, you can simply inherit from
   unsupported_serialize. Otherwise, you will need write a save/load
   function pair for your struct.

   To write a save/load function see the commented region in the
   struct.  The serialization mechanism is simple to use and it
   understands all basic datatypes as well as standard STL
   containers. If the STL container contains non-basic datatypes (such
   as a struct), save/load functions must be written for the datatype.
   */
struct vertex_data {
  size_t     num_flips;
  color_type color; 
  void save(graphlab::oarchive &oarc) const {
    oarc << num_flips << color;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> num_flips >> color;
  }
};


/** 
 * In this example, we do not need edge data. However GraphLab
 * currently does not have a mechanism to completely disable the use
 * of edge data.  Therefore, we will just put an empty struct as the
 * edge type.  
 */
struct edge_data { 
  // No edge data required
  void save(graphlab::oarchive &oarc) const { }
  void load(graphlab::iarchive &iarc) { }
};


/**
   The GraphLab graph is templatized over the vertex data as well as the
   edge data.  Here we define the type of the graph for convenience.  */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::graph<vertex_data, edge_data> memory_graph_type;



/**
   Now we can begin to write the update function class. GraphLab V2
   supports only one update function but since the update function is
   now a class and can hold state, it is easy to simulate multiple
   different update functions.

   \param context 

   The context provides access to a local neighborhood of a vertex
   (context.vertex_id()) in the graph.  The context includes all
   adjacent edges and vertices.  All vertex and edge ids are
   identified by the vertex_id_type and edge_id_type.


   GraphLab guarantees that all vertices are sequentially numbered
   from 0 (so the largest vertex id is |num_vertices| - 1), and
   similarly for edges.  All edges are directed.

   \param scheduler
   There are two basic types of schedulers.
   The synchronous / round_robin scheduler takes a single fixed set of tasks
   and repeatedly executes them until some termination condition is achieved.
   Using these schedulers generally means that the update function will not use
   this parameter.

   The task schedulers, which include fifo, multiqueue_fifo, priority,
   clustered_priority, all operate on the idea that executing an
   update_function on a vertex can be thought of as a task. Update functions
   can therefore inject new jobs into the task scheduler through this parameter.
   Since the task scheduler is slightly more complex to use, in this example,
   we will demonstrate task schedulers. Each update will decide whether to 
   schedule its neighbors, and the algorithm terminates when there are no
   tasks remaining. (See the comments in update_function for details).

   There are other methods for terminating execution, such as registering a
   termination evaluator with the engine, but we are not going to describe
   that here.

*/
struct update_functor : 
  public graphlab::iupdate_functor<graph_type, update_functor> {
  void operator()(icontext_type& context) {
    //context.vertex_data allows me to grab a reference to the vertex
    // data on the graph
    vertex_data& curvdata = context.vertex_data();
    // the in_edge_ids() function provide a vector of the edge ids of
    // the edges entering the current vertex
    const graph_type::edge_list_type in_edges = context.in_edge_ids();
    // a counter for the number of red neighbors
    size_t num_red_neighbors = 0;  
    for (size_t i = 0; i < in_edges.size(); ++i) {
      // eid is the current edge id
      graph_type::edge_id_type eid = in_edges[i];    
      // the target(eid) function allows to get the vertex at the destination
      // of the edge 'eid'. The source(eid) function provides me with the
      // source vertex.. Since I am looking at in_edges, the source vertex
      // will be my adjacent vertices
      graph_type::vertex_id_type sourcev = context.source(eid);
      // the neighbor_vertex_data() function allow me to read the
      // vertex data of a vertex adjacent to the current vertex.
      // since I am not going to change this data, I can just grab a
      // const reference. You should always try to use const
      // references whenever you know that you will definitely not be
      // changing the data, since GraphLab could make use of this
      // knowledge to perform other optimizations
      const vertex_data& nbrvertex = context.vertex_data(sourcev);
      // if red, add to our counter
      if (nbrvertex.color == RED) ++num_red_neighbors;
    }
    // get the total number of neighbors we have
    const double num_neighbors = in_edges.size();

    // Determine the new color by drawing a random number.  There are 2
    // functions. rand01() provides a random floating point number
    // between 0 and 1. rand_int(max) provides a random integer between
    // 0 and max inclusive
    color_type new_color = 
      graphlab::random::bernoulli(num_red_neighbors / num_neighbors)?
      RED : BLACK;
  
    // Determine if the coin was deterministic probability 1 or 0 of
    // landing red
    const bool is_deterministic =
      num_neighbors == num_red_neighbors || num_red_neighbors == 0;

    // see if I flip and update the current vertex data.
    const bool color_changed = (new_color != curvdata.color);
    if (color_changed) ++curvdata.num_flips;

    // Assign the new color
    curvdata.color = new_color;

    // If I flipped, all my neighbors could be affected, loop through
    // all my neighboring vertices and add them as tasks.
    if (color_changed) {
      for (size_t i = 0; i < in_edges.size(); ++i) {
        const graph_type::vertex_id_type sourcev = context.source(in_edges[i]);
        // add the task the gl::update_task object takes a vertex id,
        // and the update function to execute on. add_task also takes
        // another argument, which is the priority of this task. This
        // value should be strictly > 0.  The priority parameter of
        // course, is only used by the priority schedulers. In this
        // demo app, we don't really care about the priority, so we
        // will just set it to 1.0
        context.schedule(sourcev, update_functor());
      }
    }
    // now if I flipped myself based on a random number. This means
    // that if I update myself again, I could switch colors. Therefore
    // I should add myself as a task
    if (is_deterministic == false) {
      const update_functor fun;
      context.schedule(context.vertex_id(), fun);
    }
  }
  
  void save(graphlab::oarchive &oarc) const { }
  void load(graphlab::iarchive &iarc) { }
}; // end of update_functor




// A sync is defined by a pair of functions, a reducer, and an apply
// The reducer is exactly a fold over all the vertices, and the apply
// takes the final value at the end of the reduce, and performs whatever
// transformation it needs, before writing it into the Shared Data table.
//
// for instance, an L2 sum can be computed by having the reducer add the
// squares of the values at each vertex, then the apply function performs
// the square root.
//
// We will use this to implement the RED_PROPORTION sync. The way we will
// implement this is to use the reducer to count the number of Red
// vertices. The apply function will then divide the result by the value in
// the NUM_VERTICES table entry

/**
   This is the reducer for the RED_PROPORTION sync.
   We just count the number of red verices.

   \param context The context on the vertex we are currently accessing

   \param accumulator The input and output of the fold/reduce operation.
*/       
class accumulator :
  public graphlab::iaccumulator<graph_type, update_functor, accumulator> {
private:
  size_t red_count, flips_count;
public:
  accumulator() : red_count(0), flips_count(0) { }
  void operator()(icontext_type& context) {
    red_count += (context.vertex_data().color == RED)? 1 : 0;
    flips_count += context.vertex_data().num_flips;
  }
  void operator+=(const accumulator& other) { 
    red_count += other.red_count; 
    flips_count += other.flips_count; 
  }
  void finalize(iglobal_context_type& context) {
    const size_t numvertices = context.num_vertices();
    const double proportion = double(red_count) / numvertices;
    // here we can output something as a progress monitor
    std::cout << "Red Proportion: " << proportion << std::endl
              << "Num Flips: " << flips_count << std::endl;
    // write the final result into the shared data table
    context.set_global("RED_PROPORTION", proportion);
    context.set_global("NUM_FLIPS", flips_count);
  }
}; // end of  accumulator


/**
   In this function, we construct the grid graph
*/
void init_graph(memory_graph_type& g,
                size_t dim) {
  // here we create dim * dim vertices.  The graph
  // add_vertex(vertexdata) function takes the vertex data as input
  // and returns the vertex id of the new vertex.  The ids are
  // guaranteed to be sequentially numbered
  for (size_t i = 0;i < dim * dim; ++i) {
    // create the vertex data, randomizing the color
    vertex_data vdata;
    vdata.num_flips = 0;
    vdata.color = (graphlab::random::bernoulli())? RED : BLACK;
    // create the vertex
    g.add_vertex(vdata);
  }

  // create the edges. The add_edge(i,j,edgedata) function creates
  // an edge from i->j. with the edgedata attached. It then returns the id
  // of the new edge. The ids are guaranteed to be sequentially numbered.
  // GraphLab does NOT support duplicated edges, and currently has no facilities
  // for checking for accidental duplicated edge insertions at the
  // graph construction stage. (It is quite costly to do so)
  //
  // Any duplicated edges will result in an assertion failure at the later
  // 'finalize' stage.
  edge_data edata;
  for (size_t i = 0;i < dim; ++i) {
    for (size_t j = 0;j < dim - 1; ++j) {
      // add the horizontal edges in both directions
      g.add_edge(dim * i + j, dim * i + j + 1, edata);
      g.add_edge(dim * i + j + 1, dim * i + j, edata);

      // add the vertical edges in both directions
      g.add_edge(dim * j + i, dim * (j + 1) + i, edata);
      g.add_edge(dim * (j + 1) + i, dim * j + i, edata);
    }
  }

  // the graph is now constructed
  // we need to call finalize. 
  g.finalize();
}



int main(int argc,  char *argv[]) {

  // sets the logging level of graphlab
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  
  // Seed the random number generator
  graphlab::random::nondet_seed();

  // Parse the command line using the command line options tool
  // and context type on the command line
  graphlab::command_line_options opts;
  
  size_t dimensions = 20;
  bool makegraph = false;
  opts.use_distributed_options();
  opts.attach_option("makegraph", 
         &makegraph, makegraph, 
         "Makes Graph");
  opts.attach_option("dim", 
		     &dimensions, size_t(20), 
		     "the dimension of the grid");
  // parse the command line
  bool success = opts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }
  
  
  if (makegraph) {
    // call init_graph to create the graph
    memory_graph_type g;
    init_graph(g, dimensions);
    std::vector<graphlab::graph_partitioner::part_id_type> parts;
    graphlab::graph_partitioner::graph_partitioner::metis_partition(g, 16, parts);
    graphlab::disk_graph<vertex_data, edge_data> dg("demograph", 16,  
                            graphlab::disk_graph_atom_type::WRITE_ONLY_ATOM);
    dg.create_from_graph(g, parts);
    return 0;
  }

  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param param;
  ASSERT_TRUE(graphlab::init_param_from_mpi(param));
  graphlab::distributed_control dc(param);

  // Display the values
  opts.print();

  // create a graphlab core which contains the graph, shared data, and
  // engine
  graphlab::distributed_core<graph_type, update_functor> glcore(dc, "demograph.idx");

  // Initialize the core with the command line arguments
  glcore.set_options(opts);

  glcore.build_engine();
  // Initialize the shared data.
  accumulator initial_accum;
  size_t sync_interval = 100;  
  glcore.add_sync("sync", initial_accum, sync_interval);
  glcore.add_global("NUM_FLIPS", size_t(0));
  glcore.add_global("RED_PROPORTION", double(0));


  // since we are using a task scheduler, we need to
  // to create tasks. otherwise the engine will just terminate immediately
  // there are DIM * DIM vertices
  const update_functor functor;
  glcore.schedule_all( functor );
  
  // Run the graphlab engine 
  double runtime = glcore.start();

  // output the runtime
  std::cout << "Completed in " << runtime << " seconds" << std::endl;

  // since it is possible for the engine to terminate in between syncs
  // if we want to get a correct value for the syncs we should run them again
  // we can do his with
  glcore.sync_now("sync");

  // now we can look the values using the get() function
  size_t numberofflips = glcore.get_global<size_t>("NUM_FLIPS");
  double redprop = glcore.get_global<double>("RED_PROPORTION");

  // output some interesting statistics
  std::cout << "Number of flips: " <<  numberofflips << std::endl;
  std::cout << "Red prop: " << redprop << std::endl;

  // output the graph
  // note that here we take advantage of the fact that vertex insertion
  // gives sequential numberings
  size_t ctr = 0;
  for (size_t i = 0;i < dimensions; ++i) {
    for (size_t j = 0;j < dimensions; ++j) {
      std::cout << size_t(glcore.graph().get_vertex_data(ctr).color) << " ";
      ++ctr;
    }
    std::cout << std::endl;
  }

}


/*
  As a final comment. Since the update function only requires reading of
  neighboring vertex data, the edge_consistency model is guaranteed to have
  sequential consistency, and the algorithm is therefore guaranteed to be
  correct if executed with --scope=edge or --scope=full.

  Sequential consistency is not guaranteed under --scope=vertex, though it could
  be quite difficult in practice to construct a race.
*/

