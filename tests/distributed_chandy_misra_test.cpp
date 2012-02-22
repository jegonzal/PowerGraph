/**
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */



#include <vector>
#include <string>
#include <boost/unordered_set.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/distributed_chandy_misra.hpp>
#include <graphlab/graph/distributed_graph.hpp>
// This is messy
#include <graphlab/../../toolkits/shared/io.hpp>


#include <graphlab/macros_def.hpp>

graphlab::mutex mt;
boost::unordered_set<graphlab::vertex_id_type> lockedset;

struct vertex_data {
  uint32_t nupdates;
  double value, old_value;
  vertex_data(double value = 1) :
    nupdates(0), value(value), old_value(0) { }
}; // End of vertex data
SERIALIZABLE_POD(vertex_data);
std::ostream& operator<<(std::ostream& out, const vertex_data& vdata) {
  return out << "Rank=" << vdata.value;
}

struct edge_data { }; // End of edge data
SERIALIZABLE_POD(edge_data);

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


graphlab::distributed_chandy_misra<graph_type> *locks;
graph_type *ggraph;

void callback(graphlab::vertex_id_type v) {
  logstream(LOG_INFO) << "Locked " << ggraph->global_vid(v) << std::endl;
  mt.lock();
  ASSERT_TRUE(lockedset.find(v) == lockedset.end());
  lockedset.insert(v);
  mt.unlock();
//  graphlab::my_sleep(1);
  locks->philosopher_stops_eating(v);
}



int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);


  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_dir = "/mnt/bigbrofs/usr10/haijieg/domain_graph/edata_splits/";
  std::string format = "adj";
  clopts.attach_option("graph", &graph_dir, graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", &format, format,
                       "The graph file format: {metis, snap, tsv, adj, bin}");
  size_t ring = 0;
  clopts.attach_option("ring", &ring, ring,
                       "The size of the ring. "
                       "If ring=0 then the graph file is used.");
  size_t randomconnect = 0;
  clopts.attach_option("randomconnect", &randomconnect, randomconnect,
                       "The size of a randomly connected network. "
                       "If randomconnect=0 then the graph file is used.");

  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << dc.procid() << ": Starting." << std::endl;
  graphlab::timer timer; timer.start();
  graph_type graph(dc);
  ggraph = &graph;
  if(ring > 0) {
    if(dc.procid() == 0) {
      for(size_t i = 0; i < ring; ++i) graph.add_edge(i, i + 1);
      graph.add_edge(ring, 0);
    }
  } else if(randomconnect > 0) {
    if(dc.procid() == 0) {
      for(size_t i = 0; i < randomconnect; ++i) {
        std::vector<bool> v(randomconnect, false);
        v[i] = true;
        for (size_t r = 0; r < randomconnect /2 ; ++r) {
          size_t t = graphlab::random::rand() % randomconnect;
          if (v[t] == false && t > i) {
            graph.add_edge(i, t);
            std::cout << i << "->" << t << "\n";
            v[t] = true;
          }
        }
      }
    }
  } else {
    std::vector<std::string> graph_files;
    graphlab::fs_util::list_files_with_prefix(graph_dir, "", graph_files);
    const bool gzip = false;
    for(size_t i = 0; i < graph_files.size(); ++i) {
      if (i % dc.numprocs() == dc.procid()) {
        const std::string graph_fname = graph_dir + graph_files[i];
        std::cout << "Loading graph from structure file: " << graph_fname
                  << std::endl;
        graphlab::graph_ops<graph_type>::load_structure(graph_fname,
                                                        format, graph, gzip);
      }
    }
  }
  std::cout << dc.procid() << ": Enter Finalize" << std::endl;
  graph.finalize();
  std::cout << " ==============================================================="
            << std::endl;
  std::cout << dc.procid() << ": Finished in " << timer.current_time() << std::endl;

  std::cout
    << "========== Graph statistics on proc " << dc.procid()
    << " ==============="
    << "\n Num vertices: " << graph.num_vertices()
    << "\n Num edges: " << graph.num_edges()
    << "\n Num replica: " << graph.num_replica()
    << "\n Replica to vertex ratio: "
    << (float)graph.num_replica()/graph.num_vertices()
    << "\n --------------------------------------------"
    << "\n Num local own vertices: " << graph.num_local_own_vertices()
    << "\n Num local vertices: " << graph.num_local_vertices()
    << "\n Replica to own ratio: "
    << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
    << "\n Num local edges: " << graph.num_local_edges()
    << "\n Begin edge id: " << graph.global_eid(0)
    << "\n Edge balance ratio: " << (float)graph.num_local_edges()/graph.num_edges()
    << std::endl;

  for (graphlab::vertex_id_type v = 0; v < graph.num_local_vertices(); ++v) {
    std::cout << graph.l_get_vertex_record(v).gvid << ": " << graph.l_get_vertex_record(v).owner << ":";
    foreach(graphlab::procid_t pid,  graph.l_get_vertex_record(v).get_replicas()) {
      std::cout << pid << " ";
    }
    std::cout << "\n";
  }
  dc.barrier();
  locks = new graphlab::distributed_chandy_misra<graph_type>(dc, graph, callback);

  dc.full_barrier();
  for (graphlab::vertex_id_type v = 0; v < graph.num_local_vertices(); ++v) {
    if (graph.l_get_vertex_record(v).owner == dc.procid()) {
      std::cout << dc.procid() << ": Lock Req for " << graph.l_get_vertex_record(v).gvid << std::endl;
      locks->make_philosopher_hungry(v);
    }
  }
  if (dc.procid() == 0) {
    getchar();
  }
  dc.barrier();
  for (graphlab::vertex_id_type v = 0; v < graph.num_local_vertices(); ++v) {
    if (graph.l_get_vertex_record(v).owner == dc.procid()) {
      if(lockedset.find(v) == lockedset.end()) {
          std::cout << graph.l_get_vertex_record(v).gvid << " not acquired\n";
      }
    }
  }
  dc.barrier();
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


