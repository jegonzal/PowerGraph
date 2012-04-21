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



struct vertex_data {
  float value;
  float old_value;
  vertex_data(float value = 1) : value(value),old_value(0) { }
  void save(graphlab::oarchive &oarc) const { oarc << value << old_value; }
  void load(graphlab::iarchive &iarc) { iarc >> value >> old_value; }
}; // End of vertex data

std::ostream& operator<<(std::ostream& out, const vertex_data& vdata) {
  return out << "Rank=" << vdata.value;
}

struct edge_data : public graphlab::IS_POD_TYPE { }; // End of edge data


typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


//! Global random reset probability
double RESET_PROB = 0.15;

//! Global accuracy tolerance
double ACCURACY = 1e-5;

/**
 * The factorized page rank update function
 */
class factorized_pagerank : 
  public graphlab::iupdate_functor<graph_type, factorized_pagerank> {
private:
  float accum;
public:
  factorized_pagerank(const float& accum = 0) : accum(accum) { }
  double priority() const { return accum; }
  void operator+=(const factorized_pagerank& other) { accum += other.accum; }
  bool is_factorizable() const { return true; }
  consistency_model consistency() const { return graphlab::DEFAULT_CONSISTENCY; }
  consistency_model gather_consistency() { return graphlab::EDGE_CONSISTENCY; }
  consistency_model scatter_consistency() { return graphlab::NULL_CONSISTENCY; }
  edge_set gather_edges() const { return graphlab::IN_EDGES; }
  edge_set scatter_edges() const {
    return (accum != 0) ? graphlab::OUT_EDGES : graphlab::NO_EDGES;
  }

  // Reset the accumulator before running the gather
  void init_gather(icontext_type& context) { accum = 0; }


  void save(graphlab::oarchive &oarc) const {
    oarc << accum;
  };
  void load(graphlab::iarchive &iarc) {
    iarc >> accum;
  };
  
  // Run the gather operation over all in edges
  void gather(icontext_type& context, const edge_type& edge) {
    const size_t num_out_edges = context.num_out_edges(edge.source());
    const float weight =  1.0 / float(num_out_edges);
    accum += context.const_vertex_data(edge.source()).value * weight;
  } // end of gather

  // Merge two factorized_pagerank accumulators after running gather
  void merge(const factorized_pagerank& other) { accum += other.accum; }

  // Update the center vertex
  void apply(icontext_type& context) {
    vertex_data& vdata = context.vertex_data();
    vdata.value = RESET_PROB + (1 - RESET_PROB) * accum;
    if (std::fabs(vdata.value - vdata.old_value) >= ACCURACY) {
      accum = std::fabs(vdata.value - vdata.old_value);
      vdata.old_value = vdata.value;
    }
    else {
      accum = 0;
    }
  } // end of apply

  // Reschedule neighbors 
  void scatter(icontext_type& context, const edge_type& edge) {
    if (accum != 0) context.schedule(edge.target(), factorized_pagerank(accum));
  } // end of scatter
}; // end of factorized_pagerank update functor



/**
 * The factorized page rank update function
 */
class synchronous_pagerank :
  public graphlab::iupdate_functor<graph_type, synchronous_pagerank>  {
private:
  float accum;
public:
  synchronous_pagerank(const float& accum = 0) : accum(accum) { }
  double priority() const { return accum; }
  void operator+=(const synchronous_pagerank& other) { accum += other.accum; }
  bool is_factorizable() const { return true; }
  consistency_model consistency() const { return graphlab::DEFAULT_CONSISTENCY; }
  consistency_model gather_consistency() { return graphlab::EDGE_CONSISTENCY; }
  consistency_model scatter_consistency() { return graphlab::NULL_CONSISTENCY; }
  edge_set gather_edges() const { return graphlab::IN_EDGES; }
  edge_set scatter_edges() const {
    return graphlab::NO_EDGES;
  }
  void save(graphlab::oarchive &oarc) const {
    oarc << accum;
  };
  void load(graphlab::iarchive &iarc) {
    iarc >> accum;
  };

  // Reset the accumulator before running the gather
  void init_gather(icontext_type& context) { accum = 0; }

  // Run the gather operation over all in edges
  void gather(icontext_type& context, const edge_type& edge) {
    const size_t num_out_edges = context.num_out_edges(edge.source());
    ASSERT_EQ(edge.target(), context.vertex_id());
    ASSERT_GT(num_out_edges, 0);
    const float weight =  1.0 / float(num_out_edges);
    accum += context.const_vertex_data(edge.source()).value * weight;
  } // end of gather

  // Merge two factorized_pagerank accumulators after running gather
  void merge(const synchronous_pagerank& other) { accum += other.accum; }

  // Update the center vertex
  void apply(icontext_type& context) {
    vertex_data& vdata = context.vertex_data();
    vdata.old_value = vdata.value;
    vdata.value =  RESET_PROB + (1 - RESET_PROB) * accum;
  } // end of apply

  // Reschedule neighbors
  void scatter(icontext_type& context, const edge_type& edge) {
  } // end of scatter
}; // end of factorized_pagerank update functor










#if defined(FSCOPE)
typedef graphlab::distributed_fscope_engine<graph_type, factorized_pagerank> engine_type;
#elif defined(SYNCHRONOUS_ENGINE)
typedef graphlab::distributed_synchronous_engine<graph_type, factorized_pagerank> engine_type;
#else
typedef graphlab::distributed_engine<graph_type, factorized_pagerank> engine_type;
#endif


typedef graphlab::distributed_synchronous_engine<graph_type, synchronous_pagerank> synchronous_engine_type;


/**
 * This aggregator finds the number of touched vertices
 */
class sum_residual_aggregator :
  public graphlab::iaggregator<graph_type, factorized_pagerank, sum_residual_aggregator>,
  public graphlab::IS_POD_TYPE {
private:
  float error_norm;
  float vector_norm;
public:
  sum_residual_aggregator() : error_norm(0),vector_norm(0) { }
  void operator()(icontext_type& context) {
    float e = context.const_vertex_data().value * context.const_vertex_data().old_value;
    error_norm += e;
    vector_norm += context.const_vertex_data().old_value * context.const_vertex_data().old_value;
  } // end of operator()
  void operator+=(const sum_residual_aggregator& other) {
    error_norm += other.error_norm;
    vector_norm += other.vector_norm;
  }
  void finalize(iglobal_context_type& context) {
    std::cout << "|x'Ax|/|x'x| :\t\t" << (error_norm) / (vector_norm) << std::endl;
  }
}; //


int main(int argc, char** argv) {
  //global_logger().set_log_level(LOG_DEBUG);
  //global_logger().set_log_to_console(true);

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
 

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
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
  size_t lognormal = 0;
  clopts.attach_option("lognormal", &lognormal, lognormal,
                       "Generate a synthetic lognormal out-degree graph. ");

  size_t powerlaw = 0;
  clopts.attach_option("powerlaw", &powerlaw, powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");


  size_t randomconnect = 0;
  clopts.attach_option("randomconnect", &randomconnect, randomconnect,
                       "The size of a randomly connected network. "
                       "If randomconnect=0 then the graph file is used.");
  bool output = false;
  clopts.attach_option("output", &output, output,
                       "Output results");
  clopts.attach_option("accuracy", &ACCURACY, ACCURACY,
                       "residual termination threshold");
  clopts.attach_option("resetprob", &RESET_PROB, RESET_PROB,
                       "Random reset probability"); 

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

  std::cout << dc.procid() << ": Starting." << std::endl;
  graphlab::timer timer; timer.start();
  graph_type graph(dc, clopts);
  if(powerlaw > 0) {
    graph.build_powerlaw(powerlaw);
  } else if(lognormal > 0) {
    graph.build_lognormal(lognormal);
  } else if(ring > 0) {
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
  engine.add_aggregator("residual", sum_residual_aggregator(), 3);
  std::cout << dc.procid() << ": Scheduling all" << std::endl;
  engine.schedule_all(factorized_pagerank(1.0));
  dc.full_barrier();
  
  // Run the PageRank ---------------------------------------------------------

  std::cout << "Running pagerank!" << std::endl;
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

  std::cout << "Running one synchronous iteration for the error...";
  // run one more synchronous round to check the error
  synchronous_engine_type sync_engine(dc, graph, clopts.get_ncpus());
  sync_engine.set_options(clopts);
  sync_engine.initialize();
  sync_engine.schedule_all(synchronous_pagerank(1.0));
  dc.full_barrier();
  sync_engine.start();
  
  engine.aggregate_now("residual");
  
  if (output) {
    std::string fname = "results_";
    fname = fname + graphlab::tostr((size_t)dc.procid());
    std::ofstream fout(fname.c_str());
    for (size_t i = 0;i < graph.get_local_graph().num_vertices(); ++i) {
      if (graph.l_get_vertex_record(i).owner == dc.procid()) {
        fout << graph.l_get_vertex_record(i).gvid << "\t" 
             << graph.get_local_graph().vertex_data(i).value << "\n";
        //             << graph.get_local_graph().vertex_data(i).nupdates << "\n";
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


