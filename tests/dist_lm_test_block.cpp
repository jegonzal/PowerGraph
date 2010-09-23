#include <graphlab/distributed/lock_manager.hpp>
#include <graphlab/distributed/graph_lock_manager.hpp>
#include <graphlab/distributed/graph/cloned_graph.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/timer.hpp>
#include <vector>
#include <algorithm>
using namespace graphlab;

std::vector<size_t> locallockrequests;
std::vector<size_t> done;



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

void lock_response(distributed_control& dc, size_t source,
                  void* unused, size_t len, handlerarg_t requestid,
                  std::map<vertex_id_t, graph_type::vertex_data_type> &vdata,
                  std::map<vertex_id_t, graph_type::edge_data_type> &edata) {
/*  logstream(LOG_INFO) << dc.procid() << ": Acquired lock " << requestid
                      << " from " << source << std::endl;*/
  // ask to release the lock
  done[requestid] = 1;

}

void generate_random_locks(distributed_control &dc,
                           graph_lock_manager<graph_type> &glm,
                           size_t numvertices,
                           size_t numreqs) {
  std::vector<dist_scope_request> reqs;
  for (size_t reqid = 0;reqid < numreqs; ++reqid) {
    // pick a vertex
    size_t v = rand() % numvertices;
    reqs.push_back(dist_scope_request(v, scope_range::VERTEX_CONSISTENCY));
  }
  size_t blockid = glm.block_deferred_lock(reqs);
  size_t numlocksdone = 0;
  while(1) {
    std::vector<dist_scope_request> newlocks;
    size_t locksremaining = glm.block_status(blockid, newlocks);
    if (newlocks.size() > 0) {
      numlocksdone += newlocks.size();
      glm.block_release_partial(blockid, newlocks);
    }
    if (locksremaining == 0) break;
    else sched_yield();
  }
  glm.block_release(blockid);
  ASSERT_EQ(numlocksdone, numreqs);
}

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

  const size_t N = 1000;
  graphlab::cloned_graph<vertex_data, edge_data> g;
  create_graph(dc, g, N);

  distributed_lock_manager<graph_type> lockmanager(dc, g);
  graph_lock_manager<graph_type> graphlock(dc, lockmanager, g);
  dc.barrier();

  timer ti;
  ti.start();
  size_t numreqs = 1000;
  generate_random_locks(dc, graphlock, N, numreqs);
  double t = ti.current_time();
  //sleep(10);
  dc.barrier();
  sleep(2);
  if (dc.procid() == 0) {
    logstream(LOG_INFO) << dc.numprocs() * numreqs << " locks acquired in " << t << " seconds" << std::endl;
  }

}
