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
#include <fstream>


#include <distributed_graphlab.hpp>
#include <graphlab/engine/distributed_synchronous_engine.hpp>


#include <graphlab/util/stl_util.hpp>
#include <graphlab/macros_def.hpp>



struct vertex_data : public graphlab::IS_POD_TYPE {
  uint32_t dist;
  vertex_data() : dist(std::numeric_limits<uint32_t>::max()) { }
}; // end of vertex_data
//struct edge_data : public graphlab::IS_POD_TYPE { }; // End of edge data
typedef char edge_data;


typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



class delta_sssp :   
  public graphlab::iupdate_functor<graph_type, delta_sssp>,
  public graphlab::IS_POD_TYPE {
private:
  uint32_t dist;
public:
  delta_sssp(const float& dist = 0) : dist(dist) { }
  double priority() const { return -double(dist); }
  void operator+=(const delta_sssp& other){ dist = std::min(dist, other.dist); }
  bool is_factorizable() const { return true; }
  consistency_model consistency() const { return graphlab::DEFAULT_CONSISTENCY; }
  consistency_model gather_consistency() { return graphlab::EDGE_CONSISTENCY; }
  consistency_model scatter_consistency() { return graphlab::NULL_CONSISTENCY; }
  edge_set gather_edges() const { return graphlab::NO_EDGES; }
  edge_set scatter_edges() const { 
    return dist == uint32_t(-1)? graphlab::NO_EDGES : graphlab::OUT_EDGES; 
  }
  void apply(icontext_type& context) {  
    vertex_data& vdata = context.vertex_data();
    if (dist < vdata.dist){  
      vdata.dist = dist; ++dist;
    } else { dist = uint32_t(-1); }    
  }
  // Reschedule neighbors 
  void scatter(icontext_type& context, const edge_type& edge) {
    const vertex_id_type neighbor_id = edge.source() == context.vertex_id()?
      edge.target() : edge.source();
    if(context.const_vertex_data(neighbor_id).dist > dist)
      context.schedule(neighbor_id, delta_sssp(dist));
  } // end of scatter
}; // end of shortest path update functor






/**
 * This aggregator finds the number of touched vertices
 */       
class finite_distance_aggregator :
  public graphlab::iaggregator<graph_type, delta_sssp, finite_distance_aggregator>,
  public graphlab::IS_POD_TYPE {
private:
  size_t count;
public:
  finite_distance_aggregator() : count(0) { }
  void operator()(icontext_type& context) {
    count += context.const_vertex_data().dist < std::numeric_limits<uint32_t>::max();
  } // end of operator()
  void operator+=(const finite_distance_aggregator& other) {
    count += other.count;
  }
  void finalize(iglobal_context_type& context) {
    std::cout << "Touched:\t\t" << count << std::endl;
  }
}; //



#if defined(SYNCHRONOUS_ENGINE)
typedef graphlab::distributed_synchronous_engine<graph_type, delta_sssp> engine_type;
#else
typedef graphlab::distributed_fscope_engine<graph_type, delta_sssp> engine_type;
#endif




int main(int argc, char** argv) {
  //global_logger().set_log_level(LOG_DEBUG);
  //global_logger().set_log_to_console(true);

 

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("SSSP algorithm.");
  clopts.use_distributed_options();
  std::string graph_dir; 
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
  bool output = false;
  clopts.attach_option("output", &output, output,
                       "Output results");

  bool savebin = false;
  clopts.attach_option("savebin", &savebin, savebin,
                       "Option to save the graph as binary\n");

  std::string binpath = "./";
  std::string binprefix = "x";
  clopts.attach_option("binpath", &binpath, binpath,
                       "The path for save binary file\n");
  clopts.attach_option("binprefix", &binprefix, binprefix,
                       "The prefix for load/save binary file\n");

  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);


  std::cout << dc.procid() << ": Starting." << std::endl;
  graphlab::timer timer; timer.start();
  graph_type graph(dc, clopts);
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
          if (v[t] == false) {
            graph.add_edge(i, t);
            v[t] = true;
          }
        }
      }
    }
  } if (format=="bin") {
    logstream(LOG_INFO) << "Load graph from binary." << std::endl;
    graph.load(graph_dir, binprefix);
    dc.barrier();
  } else {
    std::vector<std::string> graph_files;
    if(boost::starts_with(graph_dir, "hdfs://")) {
      graphlab::hdfs hdfs;
      graph_files = hdfs.list_files(graph_dir);
    } else {
      graphlab::fs_util::list_files_with_prefix(graph_dir, "", graph_files);
      for(size_t i = 0; i < graph_files.size(); ++i)
        graph_files[i] = graph_dir + graph_files[i];
    }
    std::sort(graph_files.begin(), graph_files.end());
    for(size_t i = 0; i < graph_files.size(); ++i) {
      if (i % dc.numprocs() == dc.procid()) {
        std::cout << "Loading graph from structure file: " 
                  << graph_files[i] << std::endl;
        const bool success = 
          graphlab::graph_ops::load_structure(graph_files[i], format, graph);
        ASSERT_TRUE(success);
      }
    }
  }
  std::cout << dc.procid() << ": Enter Finalize" << std::endl;
  graph.finalize();
  std::cout << " ===============================================================" 
            << std::endl;
  std::cout << dc.procid() << ": Finished in " << timer.current_time() << std::endl;

  if(dc.procid() == 0){
    std::cout
      << "========== Complete Graph Statistics " << dc.procid() 
      << " ==============="
      << "\n Num vertices: " << graph.num_vertices()
      << "\n Num edges: " << graph.num_edges()
      << "\n Num replica: " << graph.num_replicas()
      << "\n Replica to vertex ratio: " 
      << float(graph.num_replicas())/graph.num_vertices()
      << "\n --------------------------------------------" 
      << "\n Num local own vertices: " << graph.num_local_own_vertices()
      << "\n Num local vertices: " << graph.num_local_vertices()
      << "\n Replica to own ratio: " 
      << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
      << "\n Num local edges: " << graph.num_local_edges()
      //<< "\n Begin edge id: " << graph.global_eid(0)
      << "\n Edge balance ratio: " << (float)graph.num_local_edges()/graph.num_edges()
      << std::endl;
  }


  if (savebin) {
    graph.save(binpath, binprefix);
  }


  std::cout << dc.procid() << ": Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts.get_ncpus());
  std::cout << dc.procid() << ": Intializing engine" << std::endl;
  engine.set_options(clopts);
  engine.initialize();
  engine.add_aggregator("count", finite_distance_aggregator(), 1000);
  std::cout << "Determing the highest degree vertex" << std::endl;
  const graphlab::vertex_id_type max_vid = graph.max_degree_vertex();
  if(graph.is_master(max_vid)) {
    std::cout << "Max degree vertex " << max_vid 
              << " found on " << dc.procid() 
              << " with in-degree " << graph.num_in_edges(max_vid)
              << " and out-degree " << graph.num_out_edges(max_vid) 
              << std::endl;
    engine.schedule(max_vid, delta_sssp(0));
  }
  dc.full_barrier();
  
  // Run the Sssp ---------------------------------------------------------
  std::cout << "Running sssp!" << std::endl;
  timer.start();
  engine.start();  // Run the engine

  const double runtime = timer.current_time();
  std::cout << "Graphlab finished, runtime: " << runtime << " seconds." 
            << std::endl
            << "Updates executed: " << engine.last_update_count() 
            << std::endl
            << "Update Rate (updates/second): " 
            << engine.last_update_count() / runtime
            << std::endl;

  engine.aggregate_now("count");
  
  if (output) {
    std::string fname = "results_";
    fname = fname + graphlab::tostr((size_t)dc.procid());
    std::ofstream fout(fname.c_str());
    for (size_t i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
      if (graph.l_get_vertex_record(i).owner == dc.procid()) {
        fout << graph.l_get_vertex_record(i).gvid << "\t" 
             << graph.l_get_vertex_record(i).num_in_edges + 
          graph.l_get_vertex_record(i).num_out_edges << "\t" 
             << graph.get_local_graph().vertex_data(i).dist << "\n";
      }
    }
  }
  /*
  if (output) {
    std::string fname = "results_local_";
    fname = fname + graphlab::tostr((size_t)dc.procid());
    std::ofstream fout(fname.c_str());
    for (size_t i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
      fout << graph.l_get_vertex_record(i).gvid << "\t" 
           << graph.l_get_vertex_record(i).num_in_edges + 
        graph.l_get_vertex_record(i).num_out_edges << "\t" 
           << graph.get_local_graph().vertex_data(i).value << "\t"
           << graph.get_local_graph().vertex_data(i).nupdates << "\n";
    }
  } */

/*  if (output) {
    std::string fname = "adj_";
    fname = fname + graphlab::tostr((size_t)dc.procid());
    std::ofstream fout(fname.c_str());
    typedef graph_type::local_graph_type::edge_type etype;

    for (size_t i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
      foreach(graph_type::edge_type e, graph.l_in_edges(i)) {
        fout << e.source() << "\t" << e.target() << "\n";
      }
    }
  } */
  

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main





