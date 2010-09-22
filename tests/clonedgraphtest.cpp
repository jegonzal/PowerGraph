
// Test the graph class

#include <vector>
#include <string>
#include <cmath>
#include <iostream>


#include <logger/logger.hpp>
#include <graphlab/distributed/graph/cloned_graph.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <serialization/serialization_includes.hpp>

#include <graphlab/macros_def.hpp>

using namespace graphlab;

struct vertex_data {
  size_t bias;
  size_t sum;
  void load(iarchive &arc) {
    arc >> bias >> sum;
  } 

  void save(oarchive &arc) const {
    arc << bias << sum;
  }

};

struct edge_data {
  size_t weight;
  size_t sum;
  void load(iarchive &arc) {
    arc >> weight >> sum;
  } 

  void save(oarchive &arc) const {
    arc << weight << sum;
  }
};

int main(int argc, char** argv) {
  const size_t N = 100000;
  vertex_data *verts = new vertex_data[N];
  edge_data *edges = new edge_data[N];

  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(1);
  dc.barrier();
  
  graphlab::cloned_graph<vertex_data, edge_data> g;
  // processor 0 creates graph
  if (dc.procid() == 0) {
    logger(LOG_INFO, "Checking Add Vertex");
    for(size_t i = 0; i < N; ++i) {     
      verts[i].bias = i; 
      verts[i].sum = 0;
      g.add_vertex(verts[i]);
      vertex_data& vdata = g.vertex_data(i);
      ASSERT_EQ(vdata.bias, verts[i].bias);
      ASSERT_EQ(vdata.sum, verts[i].sum);
      vdata.sum = 3;
      ASSERT_EQ(vdata.sum, g.vertex_data(i).sum);
      ASSERT_NE(vdata.sum, verts[i].sum);      
    }
    // Make a ring
    logger(LOG_INFO,"Checking Add Edge");
    for(size_t i = 0; i < N; ++i) {
      edges[i].weight = i * i;
      edges[i].sum = 0;
      size_t j = (i+1) % N;
      g.add_edge(i, j, edges[i]);
      g.finalize();
      edge_data& edata = g.edge_data(i,j);
      ASSERT_EQ(edata.weight, i * i);
      ASSERT_EQ(edata.sum, (size_t)0);
      edata.sum = 3;
      ASSERT_EQ(edata.sum, g.edge_data(i,j).sum);
      ASSERT_NE(edata.sum, edges[i].sum);    
    }   

    logger(LOG_INFO,"Checking Num vertices");
    ASSERT_EQ(g.num_vertices(),  N);
  }
  
  dc.barrier();
  g.distributed_partition(dc, partition_method::PARTITION_METIS,1);
  g.distribute(dc);
  dc.barrier();
  
  logstream(LOG_INFO) << "Proc " << dc.procid() << " has " 
                      << g.my_vertices().size() << " vertices" << std::endl;
  // everyone checks graph content
  logger(LOG_INFO, "Checking vertices");
  for(size_t i = 0; i < N; ++i) {     
    vertex_data& vdata = g.vertex_data(i);
    ASSERT_EQ(vdata.bias, i);
    ASSERT_EQ(vdata.sum, 3);
  }
  
  dc.barrier();
  // lets change something
  if (dc.procid() == 0) {
    g.vertex_data(0).bias = 100;
    g.vertex_data(0).sum = 100;
    g.edge_data(0).weight = -1;
    g.update_vertex(0);
    g.update_edge(0);
  }
  else {
    sleep(1);
  }
  dc.barrier();
  // everyone check the change
  ASSERT_EQ(g.vertex_data(0).bias, 100);
  ASSERT_EQ(g.vertex_data(0).sum, 100);
  ASSERT_EQ(g.edge_data(0).weight, -1);
}