#include <iostream>



#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>


#include <distributed_graphlab.hpp>
#include <graphlab/engine/distributed_engine2.hpp>

#include <graphlab/macros_def.hpp>
using namespace graphlab;

typedef size_t vertex_data;
typedef size_t edge_data;




typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


class color_update : 
        public graphlab::iupdate_functor<graph_type, color_update> {
 public:
  std::set<size_t> colors;
  
  void save(oarchive &oarc) const {
    oarc << colors;
  }
  void load(iarchive &iarc) {
    iarc >> colors;
  }
  
  double priority() const { return 1; }
  void operator+=(const color_update& other) { }
  void operator()(icontext_type& context) {} 
  
  edge_set gather_edges() const { return graphlab::IN_EDGES; }
  edge_set scatter_edges() const {
    return graphlab::NO_EDGES;
  }

  void gather(icontext_type& context, const edge_type& edge) {
    if (edge.source() == context.vertex_id()) {
      colors.insert(context.const_vertex_data(edge.target()));
    }
    else {
      colors.insert(context.const_vertex_data(edge.source()));
    }
  } // end of gather

  // Merge two pagerank_update accumulators after running gather
  void merge(const color_update& other) { 
    foreach(size_t c, other.colors) colors.insert(c);
  }

  // Update the center vertex
  void apply(icontext_type& context) {
    // find the smallest value not in colors
    size_t c = 0;
    while(colors.find(c) != colors.end()) ++c;
    context.vertex_data() = c;
  } // end of apply

  // Reschedule neighbors 
  void scatter(icontext_type& context, const edge_type& edge) {
  } // end of scatter

}; // end of pagerank update functor



typedef graphlab::distributed_engine<graph_type, color_update> dengine_type;

size_t index(const size_t nrows, const size_t i, const size_t j) {
  return i*nrows + j;
}



int main(int argc, char** argv) {
  std::cout << "Testing graph ingress." << std::endl;
  global_logger().set_log_level(LOG_DEBUG);
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
 
  size_t i = 1;
  dc.all_reduce(i);
  std::cout << "Reduction: " << i << std::endl;
  
  std::cout << dc.procid() << ": Starting." << std::endl;
  graph_type graph(dc);
  dengine_type engine(dc, graph, 1);
  const size_t nrows = 3;
  size_t eid = 0;
  for(size_t i = 0; i < nrows; ++i) {
    for(size_t j = 0; j < nrows; ++j) {
      if(eid++ % dc.numprocs() == dc.procid()) {
        const size_t source = index(nrows, i, j);
        if(i+1 < nrows) {
          const size_t target = index(nrows, i+1, j);
          graph.add_edge(source, target);
          graph.add_edge(target, source);
        }
        if(j+1 < nrows) {
          const size_t target = index(nrows, i, j+1);
          graph.add_edge(source, target);
          graph.add_edge(target, source);
        }      
      } // end of if(eid ...
    }
  } // end of for(i ...

  std::cout << dc.procid() << ": Enter Finalize" << std::endl;
  graph.finalize();
  
  std::cout << dc.procid() << ": Finished" << std::endl;
  std::cout << "Finished!" << std::endl;
  
  engine.initialize();
  engine.schedule_all(color_update());
  engine.start();
  
  
  for (vertex_id_type i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
    if (graph.l_get_vertex_record(i).owner == dc.procid()) {
      std::cout << graph.global_vid(i) << ": " << graph.get_local_graph().vertex_data(i) << "\n";
    }
  }
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main
