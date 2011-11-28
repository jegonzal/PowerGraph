// standard C++ headers
#include <iostream>

// includes the entire graphlab framework
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

struct vertex_data {
  size_t num_flips;
};

struct edge_data { 
  int from; 
  int to;
  edge_data (int f, int t) : from(f), to(t) {}
};

typedef graphlab::graph<vertex_data, edge_data> graph_type2;
typedef graph_type2::edge_list edge_list;
typedef graph_type2::edge_wrapper_type edge_wrapper_type;
typedef graph_type2::vertex_id_type vertex_id_type;

struct update_functor2 : public graphlab::iupdate_functor<graph_type2, update_functor2> {
};

void sparseGraphtest (graph_type2& g) {
  size_t num_v = 10;
  size_t num_e = 6;

  for (size_t i = 0; i < num_v; ++i)
  {
    vertex_data vdata;
    g.add_vertex(vdata);
  }

  g.add_edge(1,3,edge_data(1,3));
  g.add_edge(2,3,edge_data(2,3));
  g.add_edge(4,3,edge_data(4,3));
  g.add_edge(5,3,edge_data(5,3));

  g.add_edge(3,2, edge_data(3,2));
  g.add_edge(3,5, edge_data(3,5));

  g.finalize();

  ASSERT_EQ(g.num_vertices(), num_v);
  ASSERT_EQ(g.num_edges(), num_e);

  for (vertex_id_type i = 0; i < 6; ++i) {
    std::cout << i << std::endl;
    edge_list inedges = g.get_in_edges(i);
    edge_list outedges = g.get_out_edges(i);
    size_t arr_insize[] = {0,0,1,4,0,1};
    size_t arr_outsize[] = {0,1,1,2,1,1};
    if (i != 3) {
      ASSERT_EQ(inedges.size(), arr_insize[i]);
      ASSERT_EQ(outedges.size(), arr_outsize[i]);
      if (outedges.size() > 0)
      {
        ASSERT_EQ(outedges[0].src, i);
        ASSERT_EQ(outedges[0].target, 3);

        edge_data data = outedges[0].get_edge_data();
        ASSERT_EQ(data.from, i);
        ASSERT_EQ(data.to, 3);
      }
    } else {

      ASSERT_EQ(outedges.size(), 2);
      size_t arr_out[] = {2,5};
      for (size_t j = 0; j < 2; ++j) {
         edge_data data = outedges[j].get_edge_data();
         ASSERT_EQ(data.from, 3);
         ASSERT_EQ(data.to, arr_out[j]);
      }

      size_t arr_in[] = {1,2,4,5};
      ASSERT_EQ(inedges.size(), 4);
      for (size_t j = 0; j < 4; ++j) {
         edge_data data = inedges[j].get_edge_data();
         ASSERT_EQ(data.from, arr_in[j]);
         ASSERT_EQ(data.to, 3);
      }
    }
  }

  for (vertex_id_type i = 6; i < num_v; ++i) {
     edge_list inedges = g.get_in_edges(i);
     edge_list outedges = g.get_out_edges(i);
     ASSERT_EQ(0, inedges.size());
     ASSERT_EQ(0, outedges.size());
  }
}

/**
   In this function, we construct the 3 by 3 grid graph.
*/
void grid_graph_test(graph_type2& g) {
  g.resetMem();
  std::cout << "-----------Begin Grid Test--------------------" << std::endl;
  size_t dim = 3;
  size_t num_vertices = 0;
  size_t num_edge = 0;
  typedef uint32_t vertex_id_type;
  typedef uint32_t edge_id_type;

  // here we create dim * dim vertices.
  for (size_t i = 0;i < dim * dim; ++i) {
    // create the vertex data, randomizing the color
    vertex_data vdata;
    vdata.num_flips = 0;
    // create the vertex
    g.add_vertex(vdata);
    ++num_vertices;
  }

  // create the edges. The add_edge(i,j,edgedata) function creates
  // an edge from i->j. with the edgedata attached.   edge_data edata;

  for (size_t i = 0;i < dim; ++i) {
    for (size_t j = 0;j < dim - 1; ++j) {
      // add the horizontal edges in both directions
      //
      g.add_edge(dim * i + j, dim * i + j + 1, edge_data(dim*i+j, dim*i+j+1));
      g.add_edge(dim * i + j + 1, dim * i + j, edge_data(dim*i+j+1, dim*i+j));

      // add the vertical edges in both directions
      g.add_edge(dim * j + i, dim * (j + 1) + i, edge_data(dim*j+i, dim*(j+1)+i));
      g.add_edge(dim * (j + 1) + i, dim * j + i, edge_data(dim*(j+1)+i, dim*j+i));
      num_edge += 4;
    }
  }

  // the graph is now constructed
  // we need to call finalize. 
  g.finalize();

  printf("Test num_vertices()...\n");
  ASSERT_EQ(g.num_vertices(), num_vertices);
  printf("+ Pass test: num_vertices :)\n\n");

  printf("Test num_edges()...\n");
  ASSERT_EQ(g.num_edges(), num_edge);
  printf("+ Pass test: num_edges :)\n\n");

  // Symmetric graph: #inneighbor == outneighbor
  printf("Test num_in_neighbors() == num_out_neighbors() ...\n");
  for (size_t i = 0; i < num_vertices; ++i)
  {
    ASSERT_EQ(g.num_in_neighbors(i), g.num_out_neighbors(i));
  }
  ASSERT_EQ(g.num_in_neighbors(4), 4);
  ASSERT_EQ(g.num_in_neighbors(0), 2);
  printf("+ Pass test: #in = #out...\n\n");


  printf("Test iterate over in/out_edges and get edge data: \n");
  for (vertex_id_type i = 0; i < num_vertices; ++i) {
    const edge_list& out_edges = g.get_out_edges(i);
    const edge_list& in_edges = g.get_in_edges(i);

    printf("Test v: %u\n", i);
    printf("In edge ids: ");
    foreach(edge_wrapper_type ewrapper, in_edges) std::cout << "(" << g.edge_data(ewrapper.src, ewrapper.target).from << ","<< g.edge_data(ewrapper.src, ewrapper.target).to << ") ";
    std::cout <<std::endl;

    printf("Out edge ids: ");
    foreach(edge_wrapper_type ewrapper, out_edges) std::cout << "(" << g.edge_data(ewrapper.src, ewrapper.target).from << "," << g.edge_data(ewrapper.src, ewrapper.target).to << ") ";
    std::cout <<std::endl;

    foreach(edge_wrapper_type ewrapper, out_edges) {
      edge_data edata = ewrapper.get_edge_data();
      ASSERT_EQ(ewrapper.src, i);
      ASSERT_EQ(edata.from, ewrapper.src);
      ASSERT_EQ(edata.to, ewrapper.target);
    }

    foreach(edge_wrapper_type ewrapper, in_edges) {
      edge_data edata = ewrapper.get_edge_data();
      ASSERT_EQ(ewrapper.target, i);
      ASSERT_EQ(edata.from, ewrapper.src);
      ASSERT_EQ(edata.to, ewrapper.target);
    }
  }
  printf("+ Pass test: iterate edgelist and get data. :) \n");
  std::cout << "-----------End Grid Test--------------------" << std::endl;
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
  
  size_t dimensions = 4;
  opts.attach_option("dim", 
		     &dimensions, size_t(dimensions), 
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
  //
  graphlab::core<graph_type2, update_functor2> glcore2;

  // Initialize the core with the command line arguments
  glcore2.set_options(opts);
  sparseGraphtest(glcore2.graph());
  grid_graph_test(glcore2.graph());
}
#include <graphlab/macros_undef.hpp>
