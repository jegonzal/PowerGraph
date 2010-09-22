#include <graphlab/distributed/lock_manager.hpp>
#include <graphlab/distributed/graph_lock_manager.hpp>
#include <graphlab/distributed/graph/cloned_graph.hpp>
#include <graphlab/distributed/distributed_engine.hpp>
#include <graphlab/schedulers/fifo_scheduler.hpp>
#include <logger/logger.hpp>
#include <logger/assertions.hpp>
#include <graphlab/util/timer.hpp>
#include <vector>
#include <algorithm>

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



typedef graphlab::cloned_graph<vertex_data, edge_data> graph_type;

void create_graph(distributed_control &dc,
                  graph_type& g,
                  size_t N) {
  // processor 0 creates graph
  if (dc.procid() == 0) {
    for(size_t i = 0; i < N; ++i) {
      vertex_data vert;
      vert.bias = i;
      vert.sum = 0;
      g.add_vertex(vert);
    }
    // Make a ring
    for(size_t i = 0; i < N; ++i) {
      edge_data edge;
      edge.weight = i * i;
      edge.sum = 0;
      size_t j = (i+1) % N;
      g.add_edge(i, j, edge);
    }
    g.finalize();
  }
  dc.barrier();
  g.distributed_partition(dc, partition_method::PARTITION_METIS,1);
  g.distribute(dc);
  dc.barrier();
}




int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(2);

  const size_t N = 100000;
  
  graph_type g;
  
  create_graph(dc, g, N);
  typedef distributed_scheduler_wrapper<graph_type, fifo_scheduler<graph_type> > scheduler_type;
  distributed_engine<graph_type, scheduler_type> engine(dc, g, 2);
  
  dc.barrier();
}
