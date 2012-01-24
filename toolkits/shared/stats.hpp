#ifndef SHARED_STATS
#define SHARED_STATS

#include "types.hpp"

uint * histogram(uint * edge_count, int len, int howmnay){
  uint * ret = new uint[howmnay];
  memset(ret, 0, sizeof(int)*howmnay);
  for (int i=0; i<len; i++){
      assert(edge_count[i] < (uint)howmnay);
      ret[edge_count[i]]++;
   }

  uint sum = 0;
  for (int i=0; i< howmnay; i++)
    sum += ret[i];
  assert(sum == (uint)len);
  return ret;
}


template<typename graph_type>
void calc_stats_and_exit(const graph_type * g, bipartite_graph_descriptor & desc){
  int active = 0;
  int max_in_degree = 0;
  int max_out_degree = 0;
  int inedges = 0;
  int outedges = 0;
  for (int i=0; i< desc.total(); i++){
     int degree = g->out_edges(i).size() + g->in_edges(i).size();
     if (degree > 0)
       active++;
     inedges += g->in_edges(i).size();
     outedges += g->out_edges(i).size();
     max_in_degree= std::max(max_in_degree, (int)g->in_edges(i).size());
     max_out_degree = std::max(max_out_degree, (int)g->out_edges(i).size());
  }
  printf("Number of nodes with edges is %d\n", active);
  printf("Number of nodes without edges is %d\n", desc.total() - active);
  printf("Number of total edges %d\n", (int)g->num_edges());
  printf("Number of in edges %d\n", inedges);
  printf("Number of out edges %d\n", outedges);
  printf("Max out degree  %d\n", max_out_degree);
  printf("Max in degree %d\n", max_in_degree);

  exit(1);
}
template<typename multigraph>
void calc_multigraph_stats_and_exit(multigraph * g, bipartite_graph_descriptor & desc){
typedef unsigned long long ulong64;
typedef typename multigraph::vertex_data_type vertex_data;
  unsigned long active = 0;
  ulong64 max_in_degree = 0;
  ulong64 max_out_degree = 0;
  ulong64 inedges = 0;
  ulong64 outedges = 0;

  for (int j=0; j<g->num_graphs(); j++){
    g->doload(j);
    graphlab::timer mt; mt.start();
//omp_set_num_threads(8);
//#pragma omp parallel for
    for (int i=0; i< desc.total(); i++){
     inedges += g->graph(0)->num_in_edges(i);
     outedges += g->graph(0)->num_out_edges(i);
     max_in_degree= std::max(max_in_degree, (ulong64)g->graph(0)->num_in_edges(i));
     max_out_degree = std::max(max_out_degree, (ulong64)g->graph(0)->num_out_edges(i));
     }
     logstream(LOG_INFO) <<"Time taken to calc stats on graph: " << mt.current_time() << std::endl;
    g->unload(0);
  }
  for (uint j=0; j<g->num_vertices(); j++){
     if (g->vertex_data(j).active)
       active++;
  }
  printf("Number of nodes with edges is %lu\n", active);
  printf("Number of nodes without edges is %lu\n", g->num_vertices() - active);
  printf("Number of total edges %llu\n", (ulong64)g->num_edges());
  printf("Number of in edges %llu\n", inedges);
  printf("Number of out edges %llu\n", outedges);
  printf("Max out degree  %llu\n", max_out_degree);
  printf("Max in degree %llu\n", max_in_degree);

  exit(1);
}


#endif
