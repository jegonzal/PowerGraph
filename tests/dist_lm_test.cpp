#include <graphlab/distributed/lock_manager.hpp>
#include <graphlab/distributed/graph/cloned_graph.hpp>
#include <logger/logger.hpp>
#include <logger/assertions.hpp>
#include <graphlab/util/timer.hpp>
#include <vector>
#include <algorithm>

using namespace graphlab;

std::vector<procid_t> targetproc;
std::vector<std::vector<lockrequest> > locallockrequests;
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
                  std::vector<graph_type::edge_data_type> &edata) {
/*  logstream(LOG_INFO) << dc.procid() << ": Acquired lock " << requestid
                      << " from " << source << std::endl;*/
  // ask to release the lock
  done[requestid] = 1;
  
}

void generate_random_locks(distributed_control &dc,
                           size_t numvertices,
                           size_t numreqs) {
  double d = (double)(numvertices) / dc.numprocs();
  for (size_t reqid = 0;reqid < numreqs; ++reqid) {
    // pick a machine
    size_t proc = rand() % dc.numprocs();
    size_t low = (size_t)(std::ceil(proc * d));
    size_t high = std::floor((proc+1) * d) - 1;
  
    size_t numlocks = rand() % ((high - low) / 2);
    std::vector<lockrequest> lockreqs;
    for (size_t i = 0; i < numlocks; ++i) {
      lockrequest req;
      req.vertex = (rand() % (high - low)) + low;
      if (rand() % 2 == 0) req.locktype = RDLOCK;
      else req.locktype = WRLOCK;

      
      lockreqs.push_back(req);
    }
    std::sort(lockreqs.begin(), lockreqs.end());
    
    std::vector<lockrequest> uniquelockreqs;
    uniquelockreqs.resize(lockreqs.size());
    
    std::vector<lockrequest>::iterator it =
        unique_copy(lockreqs.begin(), lockreqs.end(), uniquelockreqs.begin());
    uniquelockreqs.resize( it - uniquelockreqs.begin() );
    locallockrequests.push_back(uniquelockreqs);
    done.push_back(0);
    targetproc.push_back(proc);
  }
  
  for (size_t reqid = 0;reqid < numreqs; ++reqid) {
    dc.remote_callxs(targetproc[reqid],
                     distributed_lock_manager<graph_type>::request_lock_handler,
                     NULL, 0, locallockrequests[reqid], reqid, (size_t)lock_response,
                     std::vector<edge_id_t>(), true);
  }

  while(1) {
    bool alldone = true;
    for (size_t reqid = 0; reqid < numreqs; ++reqid) {
    
      if (done[reqid] == 1) {
        dc.remote_callxs(targetproc[reqid],
                         distributed_lock_manager<graph_type>::request_unlock_handler,
                         NULL, 0, locallockrequests[reqid]);
        done[reqid] = 2;
      }
      else if (done[reqid] == 0){
        alldone = false;
      }
    }
    if (alldone) break;
    sched_yield();
  }
}

void create_graph(distributed_control &dc,
                  graph_type& g,
                  size_t N) {
  // processor 0 creates graph
  if (dc.procid() == 0) {
    logger(LOG_INFO, "Checking Add Vertex");
    for(size_t i = 0; i < N; ++i) {
      vertex_data vert;
      vert.bias = i;
      vert.sum = 0;
      g.add_vertex(vert);
    }
    // Make a ring
    logger(LOG_INFO,"Checking Add Edge");
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
  if (dc.procid() == 0) {
    std::vector<procid_t> vtx2owner(N);
    // equipartition
    double d = (double)(N) / dc.numprocs();
    for (size_t v = 0;v < N; ++v) {
      procid_t proc = v / d;
      if (proc >= dc.numprocs()) proc = dc.numprocs() - 1;
      vtx2owner[v] = proc;
    }
    g.set_partition(dc, vtx2owner);
  }
  else {
    std::vector<procid_t> unused;
    g.set_partition(dc, unused);
  }
  g.distribute(dc);
  dc.barrier();
}

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(1);

  const size_t N = 10000;
  graphlab::cloned_graph<vertex_data, edge_data> g;
  create_graph(dc, g, N);

  distributed_lock_manager<graph_type> lockmanager(dc, g);
  dc.barrier();

  timer ti;
  ti.start();
  
  generate_random_locks(dc, N, 1000);
  double t = ti.current_time();
  //sleep(10);
  dc.barrier();
  sleep(2);
  size_t count = 0;
  for (size_t i = 0;i < done.size(); ++i) {
    count += locallockrequests[i].size();
    ASSERT_EQ(done[i], 2);
  }
  logstream(LOG_INFO) << count << " locks acquired in " << t << " seconds" << std::endl;
  
}
