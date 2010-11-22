// standard C++ headers
#include <iostream>

// includes the entire graphlab framework
#include <graphlab.hpp>

struct vertex_data  {
  size_t numflips;
  bool color; // black == FALSE, red == TRUE,
  void save(oarchive& archive) const {
    archive << numflips << color;
  }
  void load(iarchive& archive) {
    archive >> numflips >> color;
  }
};

typedef char edge_data;
typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl;


void update_function(gl::iscope& scope,
                     gl::icallback& scheduler,
                     gl::ishared_data* shared_data) {
  vertex_data& curvdata = scope.vertex_data();
  const std::vector<gl::edge_id_t>& in_edges = scope.in_edge_ids();
  size_t num_red_neighbors = 0;  
  for (size_t i = 0; i < in_edges.size(); ++i) {
    size_t eid = in_edges[i];    
    size_t sourcev = scope.source(eid);
    const vertex_data& nbrvertex =
      scope.const_neighbor_vertex_data(sourcev);
    if (nbrvertex.color) ++num_red_neighbors;
  }
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



void init_graph(graph_type& g,
                  size_t dim) {
  for (size_t i = 0;i < dim * dim; ++i) {
    vertex_data vdata;
    vdata.numflips = 0;
    if (gl::random::rand_int(1) == 1)  vdata.color = true;
    else vdata.color = false;
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
  g.finalize();
}




enum shared_data_keys {
  NUM_VERTICES_KEY,
  RED_PROPORTION_KEY,
  NUM_FLIPS_KEY
};




void reduce_red_proportion(size_t index,
                           const gl::ishared_data& shared_data,
                           gl::iscope& scope,
                           graphlab::any& accumulator) {
  if (scope.vertex_data().color) accumulator.as<double>() += 1.0;
}



static void apply_red_proportion(size_t index,
                                  const gl::ishared_data& shared_data,
                                  graphlab::any& current_data,
                                  const graphlab::any& new_data) {

  size_t numvertices = shared_data.get_constant(NUM_VERTICES_KEY).as<size_t>();
  double numred = new_data.as<double>();
  double proportion = numred / numvertices;
  std::cout << "Red Proportion: " << proportion << std::endl;
  current_data = (double)proportion;
}


size_t get_flip(const vertex_data &v) {
  return v.numflips;
}




void init_shared_data(gl::thread_shared_data& sdm,
                      size_t dim) {
  sdm.set_constant(NUM_VERTICES_KEY, (size_t)(dim * dim));
  sdm.set_sync(RED_PROPORTION_KEY,
               reduce_red_proportion, 
               apply_red_proportion,  
               double(0),             
               100);                  
  sdm.set_sync(NUM_FLIPS_KEY,
               gl::sync_ops::sum<size_t, get_flip>,
               gl::apply_ops::identity<size_t>,
               size_t(0),
               100);
}





int main(int argc,  char *argv[]) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options
  graphlab::command_line_options opts;
  size_t dimensions = 20;
  opts.attach_option("dim",
                     &dimensions, dimensions,
                     "the dimension of the grid");
  opts.scheduler_type = "fifo";
  opts.scope_type = "edge";
  if(!opts.parse(argc, argv)) return EXIT_FAILURE;
  opts.print();

  gl::core glcore;
  glcore.set_engine_options(opts);

  // Initialize the the data structures
  init_graph(glcore.graph(), dimensions);
  init_shared_data(glcore.shared_data(), dimensions);
  // Add all starting tasks
  glcore.add_task_to_all(update_function, 1.0);
  
  // Run the graphlab engine 
  double runtime = glcore.start();
  
  std::cout << "Completed in " << runtime << " seconds" << std::endl;
  glcore.shared_data().sync(NUM_FLIPS_KEY);
  glcore.shared_data().sync(RED_PROPORTION_KEY);

  // now we can look the values using the get() function
  size_t numberofflips =
    glcore.shared_data().get(NUM_FLIPS_KEY).as<size_t>();
  double redprop =
    glcore.shared_data().get(RED_PROPORTION_KEY).as<double>();
  std::cout << "Number of flips: " <<  numberofflips << std::endl;
  std::cout << "Red prop: " << redprop << std::endl;

  // output the graph
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

