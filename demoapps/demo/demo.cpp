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

   This demo provides a synthetic application which makes use a good
   number of GraphLab concepts. Note that this demo app is
   intentionally built to use as many of graphlab concepts as
   possible. This may not be typical for most GraphLab applications.


   Picture of grid model:

   x----x----x
   |    |    |
   x----o----o
   |    |    |
   o----x----x


   Given a DIMxDIM undirected grid graph where each vertex is assigned
   a random color, either black or white (represented as a boolean
   value) Each vertex then makes a local random decision:

   - become white with probabilty proportionate to the number of white
   neighbors

   - become black with probabilty proportionate to the number of black
   neighbors

   Clearly, the two stable outcomes are where all vertices are black,
   or where all vertices are white. We are interested in knowing how
   many flips each vertex took on average.

   Also, since GraphLab only has directed edges, we will build the
   graph by duplicating every edge in both directions.

*/


// standard C++ headers
#include <iostream>

// includes the entire graphlab framework
#include <graphlab.hpp>


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
  size_t numflips;
  bool color;     // black == FALSE, red == TRUE,
};


/** In this example, we do not need edge data. However GraphLab
    currently does not have a mechanism to completely disable the use of
    edge data.  Therefore, we will just put an arbitrary small
    placeholder type on the edges.

    Note that we do not need to write a save/load function here since
    GraphLab's serializer already understands basic datatypes. */
typedef double edge_data;


/**
   The GraphLab graph is templatized over the vertex data as well as the
   edge data.  Here we define the type of the graph for convenience.  */
typedef graphlab::graph<vertex_data, edge_data> graph_type;

/**
   In order to take advantage of the types we must first predeclare the
   update functor. */
class update_functor;

/**
   Since graphlab is heavily templatized and can be inconvenient to use
   in its standard form, the graphlab::types structure provides
   convenient typedefed "shortcuts" to figure out the other graphlab
   types easily.  */
typedef graphlab::types<graph_type, update_functor> gl;


/**

   Now we can begin to write the update function. This is the standard
   form of an update function. You may specify more than one update
   function, but we only need one for this application.

   \param scope
   The scope provides access to a local neighborhood of a graph.
   The scope is centered on a particular vertex, ( scope.vertex() ), and includes
   all adjacent edges and vertices. \
   \ 
   All vertices are identified by an unsigned integer type
   vertex_id_type, and all edges are similarly identified by an unsigned
   integer type edge_id_type.  GraphLab guarantees that all vertices are
   sequentially numbered from 0 (so the largest vertex id is
   |num_vertices| - 1), and similarly for edges.  All edges are directed.

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
struct update_functor : public gl::iupdate_functor {
  void operator()(gl::iscope& scope, gl::icallback& callback) {
    //scope.vertex_data allows me to grab a reference to the vertex data
    // on the graph
    vertex_data& curvdata = scope.vertex_data();

    // the in_edge_ids() function provide a vector of the edge ids of the edges
    // entering the current vertex
    graph_type::edge_list_type in_edges = scope.in_edge_ids();
    // a counter for the number of red neighbors
    size_t num_red_neighbors = 0;  
    for (size_t i = 0; i < in_edges.size(); ++i) {
      // eid is the current edge id
      size_t eid = in_edges[i];    
      // the target(eid) function allows to get the vertex at the destination
      // of the edge 'eid'. The source(eid) function provides me with the
      // source vertex.. Since I am looking at in_edges, the source vertex
      // will be my adjacent vertices
      size_t sourcev = scope.source(eid);
      // the neighbor_vertex_data() function allow me to read the vertex data
      // of a vertex adjacent to the current vertex.
      // since I am not going to change this data, I can just grab a const
      // reference. You should always try to use const references whenever
      // you know that you will definitely not be changing the data, since
      // GraphLab could make use of this knowledge to perform other optimizations
      const vertex_data& nbrvertex = scope.neighbor_vertex_data(sourcev);
      // if red, add to our counter
      if (nbrvertex.color) ++num_red_neighbors;
    }
    // get the total number of neighbors we have
    size_t num_neighbors = in_edges.size();

    // Determine the new color by drawing a random number.  There are 2
    // functions. rand01() provides a random floating point number
    // between 0 and 1. rand_int(max) provides a random integer between
    // 0 and max inclusive
    bool new_color =
      graphlab::random::rand01() < (double(num_red_neighbors) / num_neighbors);
  
    // Determine if the coin was deterministic probability 1 or 0 of
    // landing red
    bool is_deterministic =
      num_neighbors == num_red_neighbors || num_red_neighbors == 0;

    // see if I flip and update the current vertex data.
    bool color_changed = new_color != curvdata.color;
    if (color_changed) ++curvdata.numflips;

    // Assign the new color
    curvdata.color = new_color;

    // If I flipped, all my neighbors could be affected, loop through
    // all my neighboring vertices and add them as tasks.
    if (color_changed) {
      for (size_t i = 0; i < in_edges.size(); ++i) {
        const gl::vertex_id sourcev = scope.source(in_edges[i]);
        // add the task the gl::update_task object takes a vertex id,
        // and the update function to execute on. add_task also takes
        // another argument, which is the priority of this task. This
        // value should be strictly > 0.  The priority parameter of
        // course, is only used by the priority schedulers. In this
        // demo app, we don't really care about the priority, so we
        // will just set it to 1.0
        callback.schedule(sourcev, update_functor());
      }
    }
    // now if I flipped myself based on a random number. This means
    // that if I update myself again, I could switch colors. Therefore
    // I should add myself as a task
    if (is_deterministic == false) {
      callback.schedule(scope.vertex(), update_functor());
    }
  }
}; // end of update_functor


/**
   In this function, we construct the grid graph
*/
void init_graph(graph_type& g,
                size_t dim) {
  // here we create dim * dim vertices.  The graph
  // add_vertex(vertexdata) function takes the vertex data as input
  // and returns the vertex id of the new vertex.  The ids are
  // guaranteed to be sequentially numbered
  for (size_t i = 0;i < dim * dim; ++i) {
    // create the vertex data, randomizing the color
    vertex_data vdata;
    vdata.numflips = 0;
    if (graphlab::random::bernoulli())  vdata.color = true;
    else vdata.color = false;
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

/*
  Say if we are interested in having an incremental counter which provides
  the total number of flips executed so far, as well as a ratio of the total
  number of red vertices vs black vertices.
  we can do this via the shared data manager's Sync mechanism.
 
  The Sync mechanism allows you to build a 'Fold / Reduce' operation
  across all the vertices in the graph, and store the results in the
  Shared Data object. The Shared Data table is essentially a big table
  mapping integer ids -> arbitrary data types

  First we need to define the entries of the table. The data we need are:
  - the total number of vertices (constant)
  - red vertex proportion      (synced)
  - the total number of flips    (synced)
  
  We will therefore define 3 entries in the Shared Data table.
*/


gl::glshared_const<size_t> NUM_VERTICES;
gl::glshared<double> RED_PROPORTION;
gl::glshared<size_t> NUM_FLIPS;




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

   \param scope The scope on the vertex we are currently accessing

   \param accumulator The input and output of the fold/reduce operation.
*/       
void reduce_red_proportion(gl::iscope& scope, double& acc) {
  // each entry in the shared_data table is a special data type called
  // graphlab::any (which is derived and modified from boost::any).
  // This allows you to store arbitrary datatypes into the shared data table,
  // with the minor caveat that the user must know EXACTLY what is the data
  // type stored at each entry. In this case, we will simply
  // store doubles.
  if (scope.vertex_data().color) acc += 1.0;
}

/**
   This is the apply for the RED_PROPORTION sync.
   We divide the accumulated value by the number of vertices

   \param current_data The current (old) value in the shared data table entry.
   Overwriting this will update the shared data table entry

   \param new_data The result of the reduce operation on all the vertices.
*/       
void apply_red_proportion(double& current_data, 
                          double& accum) {
  // get the number of vertices from the constant section of the shared data
  const size_t numvertices = NUM_VERTICES.get_val();
  // compute the proportion
  const double proportion = accum / numvertices;
  // here we can output something as a progress monitor
  std::cout << "Red Proportion: " << proportion << std::endl;
  // write the final result into the shared data table
  current_data = accum;
}


/**
   This is the merge function for the RED_PROPORTION sync
   Since it is just a sum, intermediate results simply add`
*/
void merge_red_proportion(double& target, const double& source) {
  target += source;
}



/**
   GraphLab provides a number of predefined syncing operations which allow
   simple reductions / applies to be implemented very quickly. 
   We will implement the NUM_FLIPS entry using one of these predefined
   operations. The predefined operations typically require the user to
   provide a simple function which extracts the information of interest
   from the vertex data. In this case, the numflips field.
*/
size_t get_flip(const vertex_data &v) {
  return v.numflips;
}

/**
   Here we create the shared data values
*/
void init_shared_data(gl::core &core, size_t dim) {
  // the number of vertices is a constant and is just dim * dim
  // since this is a constant we will just use the "constant" part of the table
  // using the function add_constant(index, value)
  //
  // Since the 'any' allows you to store any datatype, it is therefore good
  // practice to explicitly state the data type of the value you are storing
  // (size_t here). You will see this theme alot in all uses of the shared
  // data table
  NUM_VERTICES.set(dim*dim);
 
  // create the sync for the red_proportion entriy

  core.engine().
    set_sync(RED_PROPORTION,       // The value we are syncing
             double(0),             // the initial value for the fold/reduce
             128,                   // syncing frequency.in #updates             
             reduce_red_proportion, // the reduce function
             apply_red_proportion,  // the apply function           
             merge_red_proportion); // merge function



  // for the number of flips counter, we will demonstrate
  // the use of GraphLab's predefined reduce and apply operations
  // we will use set_sync as usual, but using something different for the
  // reduce and apply functions

  // glshared_sync_ops::sum<size_t, get_flip> is a predefined reduce operation
  // which sums over all the result of running get_flip() on all the vertex
  // data. The first template field (size_t) is the type of the accumulator.

  // glshared_apply_ops::identity<size_t> is the identity apply which directly
  // writes the result of the reduction into the shared data table entry.
  // The template field (size_t) is the type of the entry.
  
  // glshared_merge_ops::sum<size_t> simply returns the sum of intermediate results

  core.set_sync(NUM_FLIPS,  
                gl::glshared_sync_ops::sum<size_t, get_flip>,
                gl::glshared_apply_ops::identity<size_t>,
                size_t(0),
                128,
                gl::glshared_merge_ops::sum<size_t>);

}





int main(int argc,  char *argv[]) {

  // sets the logging level of graphlab
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse the command line using the command line options tool
  // and scope type on the command line
  graphlab::command_line_options opts;
  
  size_t dimensions = 20;
  opts.attach_option("dim", 
		     &dimensions, size_t(20), 
		     "the dimension of the grid");
  // parse the command line
  bool success = opts.parse(argc, argv);
  if(!success) {
    return EXIT_FAILURE;
  }
  // Display the values
  opts.print();

  // create a graphlab core which contains the graph, shared data, and
  // engine
  gl::core glcore;

  // Initialize the core with the command line arguments
  glcore.set_options(opts);
  
  // call init_graph to create the graph
  init_graph(glcore.graph(), dimensions);
  // call create shared_data to create the shared data
  init_shared_data(glcore, dimensions);

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
  glcore.sync_now(NUM_FLIPS);
  //  glcore.shared_data().sync(RED_PROPORTION_KEY);

  // now we can look the values using the get() function
  size_t numberofflips = NUM_FLIPS.get_val();
  double redprop = RED_PROPORTION.get_val();

  // output some interesting statistics
  std::cout << "Number of flips: " <<  numberofflips << std::endl;
  std::cout << "Red prop: " << redprop << std::endl;

  // output the graph
  // note that here we take advantage of the fact that vertex insertion
  // gives sequential numberings
  size_t ctr = 0;
  for (size_t i = 0;i < dimensions; ++i) {
    for (size_t j = 0;j < dimensions; ++j) {
      std::cout << size_t(glcore.graph().vertex_data(ctr).color) << " ";
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

