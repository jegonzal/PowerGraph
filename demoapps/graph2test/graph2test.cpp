// standard C++ headers
#include <iostream>

// includes the entire graphlab framework
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

struct vertex_data {
  size_t num_flips;
};

struct edge_data { 
  float random;
  edge_data (float f) : random(f){}
};

typedef graphlab::graph2<vertex_data, edge_data> graph_type2;
typedef graphlab::graph<vertex_data, edge_data> graph_type1;
typedef graph_type2::edge_list edge_list;

struct update_functor2 : public graphlab::iupdate_functor<graph_type2, update_functor2> {
};
struct update_functor1 : public graphlab::iupdate_functor<graph_type1, update_functor1> {
};


void size_test(graph_type1& g1, graph_type2& g2, unsigned int N = 1000) {
  g2.resetMem();
  g1.resetMem();

  size_t num_vertices = 0;
  size_t num_edge = 0;

  std::cout <<  "--------------Begin size test: building N-clique, Compare size of G1, G2 ---------------" << " N = " << N << std::endl;
  std::cout << "Initial graph1 size: " << g1.get_graph_size() << std::endl;
  std::cout << "Initial graph2 size: " << g2.get_graph_size() << std::endl;

  std::cout << "Adding vertices..." << std::endl;
  for (size_t i = 0; i < N; ++i) {
    vertex_data vdata;
    // create the vertex
    g2.add_vertex(vdata);
    g1.add_vertex(vdata);
    ++num_vertices;
  }

  std::cout << "Adding edges..." << std::endl;
  edge_data edata(0.0);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < i; ++j) {
      g2.add_edge(i,j, edata);
      g2.add_edge(j,i, edata);

      g1.add_edge(i,j,edata);
      g1.add_edge(j,i,edata);
      ++num_edge;
    }
  }

  std::cout << "Intermediate graph1 size: " << g1.get_graph_size() << std::endl;
  std::cout << "Intermediate graph2 size: " << g2.get_graph_size() << std::endl;

  std::cout << "Finalize graph1..." << std::endl;
  g1.finalize();

  std::cout << "Finalize graph2..." << std::endl;
  g2.finalize();

  std::cout << "Size of graph1 after finalized: " << g1.get_graph_size() << std::endl;
  std::cout << "Size of graph2 after finalized: " << g2.get_graph_size() << std::endl;
  std::cout << "Saving: " << double(g1.get_graph_size() - g2.get_graph_size())/(1024*1024) << "MB" << std::endl;
  std::cout << "----------------End size test --------------------" << std::endl;
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
  edge_data edata(0.0);

  for (size_t i = 0;i < dim; ++i) {
    for (size_t j = 0;j < dim - 1; ++j) {
      // add the horizontal edges in both directions
      g.add_edge(dim * i + j, dim * i + j + 1, edata);
      g.add_edge(dim * i + j + 1, dim * i + j, edata);

      // add the vertical edges in both directions
      g.add_edge(dim * j + i, dim * (j + 1) + i, edata);
      g.add_edge(dim * (j + 1) + i, dim * j + i, edata);
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

  // Test find
  printf("Test find...\n");
  typedef std::pair<bool, edge_id_type> edge_find_type;
  edge_find_type p0;
  for (size_t i = 0;i < dim; ++i) {
    for (size_t j = 0;j < dim - 1; ++j) {
      // add the horizontal edges in both directions
      p0 = g.find(dim * i + j + 1, dim * i + j);
      ASSERT_TRUE(p0.first);
      p0 = g.find(dim * i + j, dim * i + j + 1);
      ASSERT_TRUE(p0.first);
      p0 = g.find(dim * j + i, dim * (j + 1) + i);
      ASSERT_TRUE(p0.first);
      p0 = g.find(dim * (j + 1) + i, dim * j + i);
      ASSERT_TRUE(p0.first);
    }
  }
  edge_find_type pno = g.find(0, 5); ASSERT_FALSE(pno.first);
  pno = g.find(0, 4); ASSERT_FALSE(pno.first);
  printf("+ Pass test: find(src, dest) :)\n\n");

  // Test in vertices and out vertices
  printf("Test in_vertices(), out_vertices()... \n");
  std::vector<vertex_id_type> ins = g.in_vertices(4);
  std::vector<vertex_id_type> outs = g.out_vertices(4);
  ASSERT_EQ(ins.size(), 4);
  printf("+ Pass test: in/out vertices\n\n");

  // Test in source(), target(), rev_edge_id
  printf("Test source(), target(), rev_edge_id()...\n");
  for (size_t i = 0; i < num_edge; ++i) {
    // Test edge src/target
    vertex_id_type src = g.source(i);
    vertex_id_type dst = g.target(i);
    std::pair<bool, edge_id_type> eid_pair = g.find(src, dst);
    ASSERT_TRUE(eid_pair.first);
    ASSERT_EQ(i, eid_pair.second);

    // Test rev_edge_id
    edge_id_type reid = g.rev_edge_id(i);
    ASSERT_EQ(g.source(reid), dst);
    ASSERT_EQ(g.target(reid), src);
  }
  printf("+ Pass test: edge source, target\n\n");


  // Test edge_list
  edge_list in_edges = g.in_edge_ids(4);
  printf("Test in_edges of V4: \n");
  ASSERT_EQ(in_edges.size(), 4);
  foreach(edge_id_type eid, in_edges) {
    ASSERT_EQ(g.target(eid), 4);
  }

  printf("Test out_edges of V4: \n");
  edge_list out_edges = g.out_edge_ids(4);
  ASSERT_EQ(out_edges.size(), 4);
  foreach(edge_id_type eid, out_edges) {
    ASSERT_EQ(g.source(eid), 4);
  }
  printf("+ Pass test: in_edge_ids, out_edge_ids\n\n");

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
  graphlab::core<graph_type1, update_functor1> glcore1;
  graphlab::core<graph_type2, update_functor2> glcore2;

  // Initialize the core with the command line arguments
  glcore1.set_options(opts);
  glcore2.set_options(opts);

  size_test(glcore1.graph(), glcore2.graph(), 1000);

  grid_graph_test(glcore2.graph());
}
#include <graphlab/macros_undef.hpp>
