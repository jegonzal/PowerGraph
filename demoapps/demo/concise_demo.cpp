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


#include <iostream>
#include <graphlab.hpp>

struct vertex_data {
  size_t numflips;
  bool color;     // black == FALSE, red == TRUE,
};

typedef char edge_data;

typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl;



void update_function(gl::iscope& scope,
                     gl::icallback& scheduler) {
  vertex_data& curvdata = scope.vertex_data();


  gl::edge_list in_edges = scope.in_edge_ids();
  size_t num_red_neighbors = 0;  
  for (size_t i = 0; i < in_edges.size(); ++i) {
    // eid is the current edge id
    size_t eid = in_edges[i];    
    size_t sourcev = scope.source(eid);
    const vertex_data& nbrvertex = scope.neighbor_vertex_data(sourcev);
    if (nbrvertex.color) ++num_red_neighbors;
  }
  // get the total number of neighbors we have
  size_t num_neighbors = in_edges.size();

  bool new_color =
    gl::random::rand01() < (double(num_red_neighbors) / num_neighbors);
  
  bool is_deterministic =
    num_neighbors == num_red_neighbors || num_red_neighbors == 0;

  bool color_changed = new_color != curvdata.color;
  if (color_changed) ++curvdata.numflips;

  curvdata.color = new_color;

  if (color_changed) {
    for (size_t i = 0; i < in_edges.size(); ++i) {
      size_t sourcev = scope.source(in_edges[i]);
      scheduler.add_task(gl::update_task(sourcev, update_function),
                        1.0);
    }
  }
  if (is_deterministic == false) {
    scheduler.add_task(gl::update_task(scope.vertex(), update_function),
                        1.0);
  }
}

/**
  In this function, we construct the grid graph
*/
void init_graph(graph_type& g,
                  size_t dim) {
  for (size_t i = 0;i < dim * dim; ++i) {
    // create the vertex data, randomizing the color
    vertex_data vdata;
    vdata.numflips = 0;
    if (gl::random::rand_int(1) == 1)  vdata.color = true;
    else vdata.color = false;
    // create the vertex
    g.add_vertex(vdata);
  }

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


gl::glshared_const<size_t> NUM_VERTICES;
gl::glshared<double> RED_PROPORTION;
gl::glshared<size_t> NUM_FLIPS;


   
void reduce_red_proportion(gl::iscope& scope,
                           graphlab::any& accumulator) {
  if (scope.vertex_data().color) accumulator.as<double>() += 1.0;
}

void apply_red_proportion(graphlab::any& current_data, 
                          const graphlab::any& new_data) {
  size_t numvertices = NUM_VERTICES.get_val();

  double numred = new_data.as<double>();
  // compute the proportion
  double proportion = numred / numvertices;
  // here we can output something as a progress monitor
  std::cout << "Red Proportion: " << proportion << std::endl;
  // write the final result into the shared data table
  current_data.as<double>() = proportion;
}

void merge_red_proportion(graphlab::any& target, 
                          const graphlab::any& source) {
  target.as<double>() += source.as<double>();
}

size_t get_flip(const vertex_data &v) {
  return v.numflips;
}

void init_shared_data(gl::core &core, size_t dim) {
  NUM_VERTICES.set(dim*dim);
 
  // create the sync for the red_proportion entriy

  core.set_sync(RED_PROPORTION,       // The value we are syncing
                reduce_red_proportion, // the reduce function
                apply_red_proportion,  // the apply function
                double(0),             // the initial value for the fold/reduce
                128,                   // syncing frequency.in #updates
                merge_red_proportion); // merge function


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
  glcore.set_engine_options(opts);
  
  // call init_graph to create the graph
  init_graph(glcore.graph(), dimensions);
  // call create shared_data to create the shared data
  init_shared_data(glcore, dimensions);

  // since we are using a task scheduler, we need to
  // to create tasks. otherwise the engine will just terminate immediately
  // there are DIM * DIM vertices
  for (size_t i = 0;i < dimensions * dimensions; ++i) {
    glcore.add_task(gl::update_task(i, update_function), 1.0);
  }
  
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

