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

struct edge_data_empty {
};

typedef graphlab::graph2<vertex_data, edge_data> graph_type;
typedef graph_type::edge_list_type edge_list_type;
typedef graph_type::edge_type edge_type;

typedef graphlab::graph2<vertex_data, edge_data_empty> graph_type_e;
typedef graph_type_e::edge_list_type edge_list_type_e;
typedef graph_type_e::edge_type edge_type_e;


typedef uint32_t vertex_id_type;


struct update_functor : public graphlab::iupdate_functor<graph_type, update_functor> {};
struct update_functor_e : public graphlab::iupdate_functor<graph_type_e, update_functor_e> {};

void sparseGraphtest (graph_type& g) {
  size_t num_v = 10;
  size_t num_e = 6;

  for (size_t i = 0; i < num_v; ++i) {
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
    edge_list_type inedges = g.in_edges(i);
    edge_list_type outedges = g.out_edges(i);
    size_t arr_insize[] = {0,0,1,4,0,1};
    size_t arr_outsize[] = {0,1,1,2,1,1};
    if (i != 3) {
      ASSERT_EQ(inedges.size(), arr_insize[i]);
      ASSERT_EQ(outedges.size(), arr_outsize[i]);
      if (outedges.size() > 0)
        {
          ASSERT_EQ(outedges[0].source(), i);
          ASSERT_EQ(outedges[0].target(), 3);

          edge_data data = g.edge_data(outedges[0]);
          ASSERT_EQ(data.from, i);
          ASSERT_EQ(data.to, 3);
        }
    } else {

      ASSERT_EQ(outedges.size(), 2);
      size_t arr_out[] = {2,5};
      for (size_t j = 0; j < 2; ++j) {
        edge_data data = g.edge_data(outedges[j]);
        ASSERT_EQ(data.from, 3);
        ASSERT_EQ(data.to, arr_out[j]);
      }

      size_t arr_in[] = {1,2,4,5};
      ASSERT_EQ(inedges.size(), 4);
      for (size_t j = 0; j < 4; ++j) {
        edge_data data = g.edge_data(inedges[j]);
        ASSERT_EQ(data.from, arr_in[j]);
        ASSERT_EQ(data.to, 3);
      }
    }
  }

  for (vertex_id_type i = 6; i < num_v; ++i) {
    edge_list_type inedges = g.in_edges(i);
    edge_list_type outedges = g.out_edges(i);
    ASSERT_EQ(0, inedges.size());
    ASSERT_EQ(0, outedges.size());
  }
}

/**
   In this function, we construct the 3 by 3 grid graph.
*/
void grid_graph_test(graph_type& g) {
  g.clear_reserve();
  std::cout << "-----------Begin Grid Test--------------------" << std::endl;
  size_t dim = 3;
  size_t num_vertices = 0;
  size_t num_edge = 0;
  typedef uint32_t vertex_id_type;


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
  for (size_t i = 0; i < num_vertices; ++i) {
    ASSERT_EQ(g.in_edges(i).size(), g.out_edges(i).size());
  }
  ASSERT_EQ(g.in_edges(4).size(), 4);
  ASSERT_EQ(g.in_edges(0).size(), 2);
  printf("+ Pass test: #in = #out...\n\n");


  printf("Test iterate over in/out_edges and get edge data: \n");
  for (vertex_id_type i = 0; i < num_vertices; ++i) {
    const edge_list_type& out_edges = g.out_edges(i);
    const edge_list_type& in_edges = g.in_edges(i);

    printf("Test v: %u\n", i);
    printf("In edge ids: ");
    foreach(edge_type edge, in_edges) 
      std::cout << "(" << g.edge_data(edge.source(), edge.target()).from << ","
                << g.edge_data(edge.source(), edge.target()).to << ") ";
    std::cout <<std::endl;

    printf("Out edge ids: ");
    foreach(edge_type edge, out_edges) 
      std::cout << "(" << g.edge_data(edge.source(), edge.target()).from << "," 
                << g.edge_data(edge.source(), edge.target()).to << ") ";
    std::cout <<std::endl;

    foreach(edge_type edge, out_edges) {
      edge_data edata = g.edge_data(edge);
      ASSERT_EQ(edge.source(), i);
      ASSERT_EQ(edata.from, edge.source());
      ASSERT_EQ(edata.to, edge.target());
    }

    foreach(edge_type edge, in_edges) {
      edge_data edata = g.edge_data(edge);
      ASSERT_EQ(edge.target(), i);
      ASSERT_EQ(edata.from, edge.source());
      ASSERT_EQ(edata.to, edge.target());
    }
  }
  printf("+ Pass test: iterate edgelist and get data. :) \n");
  std::cout << "-----------End Grid Test--------------------" << std::endl;
}

/*
 * Test undirected graph
 * */
void undirected_graph_test(graph_type& g) {
  g.clear_reserve();
  std::cout << "-----------Begin Undirected Grid Test--------------------" << std::endl;
  size_t dim = 3;
  size_t num_vertices = 0;
  size_t num_edge = 0;
  typedef uint32_t vertex_id_type;

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
      // add the horizontal edges in one direction
      g.add_edge(dim * i + j, dim * i + j + 1, edge_data(dim*i+j, dim*i+j+1));
      // add the vertical edges in one direction
      g.add_edge(dim * j + i, dim * (j + 1) + i, edge_data(dim*j+i, dim*(j+1)+i));
      num_edge += 2;
    }
  }
  // the graph is now constructed
  // we need to call finalize. 
  g.set_is_directed(false);
  g.finalize();

  printf("Test num_vertices()...\n");
  ASSERT_EQ(g.num_vertices(), num_vertices);
  printf("+ Pass test: num_vertices :)\n\n");

  printf("Test num_edges()...\n");
  ASSERT_EQ(g.num_edges(), num_edge);
  printf("+ Pass test: num_edges :)\n\n");

  printf("Test iterate over in/out_edges and get edge data: \n");
  for (vertex_id_type i = 0; i < num_vertices; ++i) {
    const edge_list_type& out_edges = g.out_edges(i);
    const edge_list_type& in_edges = g.in_edges(i);

    printf("Test v: %u\n", i);
    printf("In edge ids: ");
    foreach(edge_type edge, in_edges) 
      std::cout << "(" << g.edge_data(edge.source(), edge.target()).from << ","
                << g.edge_data(edge.source(), edge.target()).to << ") ";
    std::cout <<std::endl;

    printf("Out edge ids: ");
    foreach(edge_type edge, out_edges) 
      std::cout << "(" << g.edge_data(edge.source(), edge.target()).from << "," 
                << g.edge_data(edge.source(), edge.target()).to << ") ";
    std::cout <<std::endl;

    // foreach(edge_type edge, out_edges) {
    //   edge_data edata = g.edge_data(edge);
    //   ASSERT_EQ(edge.source(), i);
    //   ASSERT_EQ(edata.from, edge.source());
    //   ASSERT_EQ(edata.to, edge.target());
    // }

    // foreach(edge_type edge, in_edges) {
    //   edge_data edata = g.edge_data(edge);
    //   ASSERT_EQ(edge.target(), i);
    //   ASSERT_EQ(edata.from, edge.source());
    //   ASSERT_EQ(edata.to, edge.target());
    // }
  }
  std::cout << "-----------End Undirected Grid Test--------------------" << std::endl;
}

/*
 * Test empty_edge graph
 */
void empty_edge_graph_test(graph_type_e& g) {
  std::cout << "-----------Begin Empty Edge Grid Test--------------------" << std::endl;
  g.clear_reserve();
  size_t dim = 3;
  size_t num_vertices = 0;
  size_t num_edge = 0;
  typedef uint32_t vertex_id_type;

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
      // add the horizontal edges in one direction
      g.add_edge(dim * i + j, dim * i + j + 1);
      // add the vertical edges in one direction
      g.add_edge(dim * j + i, dim * (j + 1) + i);
      num_edge += 2;
    }
  }
  // the graph is now constructed
  // we need to call finalize. 
  g.set_is_directed(false);
  g.finalize();

  printf("Test num_vertices()...\n");
  ASSERT_EQ(g.num_vertices(), num_vertices);
  printf("+ Pass test: num_vertices :)\n\n");

  printf("Test num_edges()...\n");
  ASSERT_EQ(g.num_edges(), num_edge);
  printf("+ Pass test: num_edges :)\n\n");

  printf("Test iterate over in/out_edges and get edge data: \n");
  for (vertex_id_type i = 0; i < num_vertices; ++i) {
    const edge_list_type_e& out_edges = g.out_edges(i);
    const edge_list_type_e& in_edges = g.in_edges(i);

    printf("Test v: %u\n", i);
    printf("In edge ids: ");
    foreach(edge_type_e edge, in_edges) 
      std::cout << "(" << edge.source() << ", " << edge.target() << ") ";
    std::cout << std::endl;

    printf("Out edge ids: ");
    foreach(edge_type_e edge, out_edges) 
      std::cout << "(" << edge.source() << ", " << edge.target() << ") ";
    std::cout <<std::endl;
  }
  std::cout << "-----------End Empty Edge Grid Test--------------------" << std::endl;
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
  graphlab::core<graph_type, update_functor> glcore;


  // Initialize the core with the command line arguments
  glcore.set_options(opts);
  sparseGraphtest(glcore.graph());
  grid_graph_test(glcore.graph());
  undirected_graph_test(glcore.graph());

  graphlab::core<graph_type_e, update_functor_e> glcore_e;
  empty_edge_graph_test(glcore_e.graph());
}
#include <graphlab/macros_undef.hpp>
