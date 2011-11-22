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

typedef graphlab::graph2<vertex_data, edge_data> graph_type2;
typedef graphlab::graph<vertex_data, edge_data> graph_type1;
typedef graph_type2::edge_list edge_list;

struct update_functor2 : public graphlab::iupdate_functor<graph_type2, update_functor2> { };
struct update_functor1 : public graphlab::iupdate_functor<graph_type1, update_functor1> { };

void sparseGraphtest (graph_type2& g) {
  size_t num_v = 5000000;
  for (size_t i = 0; i < num_v; ++i)
  {
    vertex_data vdata;
    g.add_vertex(vdata);
  }
  g.add_edge(1,3,edge_data(1,3));
  g.add_edge(2,5,edge_data(2,5));
  g.add_edge(5,8,edge_data(5,8));
  g.add_edge(8,9,edge_data(8,9));
  g.finalize();
  for (size_t i = 0; i < 4; ++i)
    std::cout << g.source(i) << " ";
  std::cout << std::endl;
}


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
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < i; ++j) {

      edge_data iedata(i,j);
      edge_data oedata(j,i);
      g2.add_edge(i,j, iedata);
      g2.add_edge(j,i, oedata);

      g1.add_edge(i,j,iedata);
      g1.add_edge(j,i,oedata);
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


  printf("Test iterate over in/out_edges and get edge data: \n");

  for (vertex_id_type i = 0; i < num_vertices; ++i) {
    const edge_list& out_eids = g.out_edge_ids(i);
    const edge_list& in_eids = g.in_edge_ids(i);
    printf("Test v: %u\n", i);
    printf("In edge ids: ");
    foreach(edge_id_type in, in_eids) std::cout << in << "(" << g.edge_data(in).from << ","<< g.edge_data(in).to << ") ";
    std::cout <<std::endl;
    foreach(edge_id_type out, out_eids) std::cout << out << "(" << g.edge_data(out).from << "," << g.edge_data(out).to << ") ";
    std::cout <<std::endl;

    foreach(edge_id_type eid, out_eids) {
      vertex_id_type o = g.target(eid);
      edge_data edata = g.edge_data(eid);
      ASSERT_EQ(edata.from, i);
      ASSERT_EQ(edata.to, o);
    }

    foreach(edge_id_type eid, in_eids) {
      vertex_id_type o = g.source(eid);
      edge_data edata = g.edge_data(eid);
      ASSERT_EQ(edata.from, o);
      ASSERT_EQ(edata.to, i);
    }

    for(size_t j = 0; j < out_eids.size(); ++j) {
      edge_id_type eid = out_eids[j];
      vertex_id_type o = g.target(eid);
      edge_data edata = g.edge_data(eid);
      ASSERT_EQ(edata.from, i);
      ASSERT_EQ(edata.to, o);
    }

    for(size_t j = 0; j < in_eids.size(); ++j) {
      edge_id_type eid = in_eids[j];
      vertex_id_type o = g.source(eid);
      edge_data edata = g.edge_data(eid);
      ASSERT_EQ(edata.from, o);
      ASSERT_EQ(edata.to, i);
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
  graphlab::core<graph_type1, update_functor1> glcore1;
  graphlab::core<graph_type2, update_functor2> glcore2;

  // Initialize the core with the command line arguments
  glcore1.set_options(opts);
  glcore2.set_options(opts);

  size_test(glcore1.graph(), glcore2.graph(), 1000);
  //sparseGraphtest(glcore2.graph());
 // grid_graph_test(glcore2.graph());
}
#include <graphlab/macros_undef.hpp>
