#ifndef SHARED_STATS
#define SHARED_STATS

#include "types.hpp"

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
  int active = 0;
  int max_in_degree = 0;
  int max_out_degree = 0;
  int inedges = 0;
  int outedges = 0;
  for (int j=0; j<g->num_graphs(); j++){
    for (int i=0; i< desc.total(); i++){
     
     int degree = g->graph(j)->out_edges(i).size() + g->graph(j)->in_edges(i).size();
     if (degree > 0)
       active++;
     inedges += g->graph(j)->in_edges(i).size();
     outedges += g->graph(j)->out_edges(i).size();
     max_in_degree= std::max(max_in_degree, (int)g->graph(j)->in_edges(i).size());
     max_out_degree = std::max(max_out_degree, (int)g->graph(j)->out_edges(i).size());
   }
  }
  printf("Number of nodes with edges is %d\n", active);
  printf("Number of nodes without edges is %d\n", g->num_vertices() - active);
  printf("Number of total edges %d\n", (int)g->num_edges());
  printf("Number of in edges %d\n", inedges);
  printf("Number of out edges %d\n", outedges);
  printf("Max out degree  %d\n", max_out_degree);
  printf("Max in degree %d\n", max_in_degree);

  exit(1);
}


#endif
