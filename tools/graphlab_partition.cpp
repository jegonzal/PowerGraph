#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>


#include <graphlab.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/distributed_chromatic_engine.hpp>
#include <graphlab/distributed2/distributed_glshared.hpp>
using namespace graphlab;


#include <graphlab/macros_def.hpp>

/// GLOBAL CONSTANTS
const size_t MAX_CHANGES(2);
const size_t MAX_ITERATIONS(1000);
const size_t SYNC_INTERVAL(100);
const size_t NUM_COLORS(10);

///////////////////////////////////////////////////////////////////////////////
///////////////////////// Types ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
struct vertex_data_type {
  procid_t atomid;
  procid_t num_changes;
  bool       is_set;
  bool       is_seed;
};
SERIALIZABLE_POD(vertex_data_type);
struct edge_data_type { };
SERIALIZABLE_POD(edge_data_type);


typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;
typedef distributed_chromatic_engine< graph_type > engine_type;
typedef engine_type::iscope_type iscope_type;
typedef engine_type::icallback_type icallback_type;
typedef ishared_data<graph_type> ishared_data_type;
typedef engine_type::icallback_type icallback_type;
typedef engine_type::update_task_type update_task_type;


size_t num_atoms(10);



// global data
typedef std::map<procid_t, double> atom2count_type;
typedef distributed_glshared<atom2count_type> shared_atom2count_type;
shared_atom2count_type global_atom2count;
// update the counts in the appropriate table
void atom2count_sum_fun(iscope_type& iscope,
                  any& acc) {
  if(iscope.const_vertex_data().is_set) 
    acc.as<atom2count_type>()[iscope.const_vertex_data().atomid]++;
}
// Identity apply
void atom2count_apply_fun(any& current_data, 
                           const any& acc) { 
  current_data.as<atom2count_type>() = 
    acc.as<atom2count_type>();
  atom2count_type& count_map(current_data.as<atom2count_type>());
  typedef atom2count_type::value_type pair_type;
  double sum = 0;
  foreach(pair_type& pair, count_map)  sum += pair.second;
  ASSERT_GT(sum, 0);
  foreach(pair_type& pair, count_map) pair.second /= sum;
} 
// Sum the two maps
void atom2count_merge_fun(any& any_dest,
                           const any& any_src) {
  typedef atom2count_type::value_type pair_type;
  const atom2count_type& src(any_src.as<atom2count_type>());
  atom2count_type& dest(any_dest.as<atom2count_type>());
  foreach(const pair_type& pair, src)
    dest[pair.first] += pair.second;
}










procid_t find_best_atom(atom2count_type& local_atom2count,
                        const atom2count_type& global_atom2count) {
  double best_score = 0;
  procid_t best_atomid = local_atom2count.begin()->first;
  ASSERT_LT(best_atomid, num_atoms);
  typedef atom2count_type::value_type pair_type;
  // normalize the local_atom2count map
  double sum(0);
  foreach(pair_type& pair, local_atom2count) sum += pair.second;
  ASSERT_GT(sum, 0);
  foreach(pair_type& pair, local_atom2count) pair.second /= sum;

  foreach(const pair_type& pair, local_atom2count) {
    const procid_t atomid(pair.first);
    ASSERT_LT(atomid, num_atoms);
    const double local_count(pair.second);
    ASSERT_GT(local_count, 0);
    const double global_count(safe_get(global_atom2count, atomid, double(0)) );
    // auto join
    if(global_count == 0) { return atomid; }
    ASSERT_GT(global_count, 0);
    // otherwise compute the 
    const double score = local_count / global_count;
    if(score > best_score)  {
      best_atomid = atomid;
      best_score = score;
    }
  }
  return best_atomid;
} // end of best join atoms


void partition_update_function(iscope_type& scope,
                               icallback_type& callback,
                               ishared_data_type* unused) {
  atom2count_type local_atom2count;
  // Get the number of neighbor assignments
  foreach(const edge_id_t eid, scope.in_edge_ids()) {
    const vertex_id_t vid(scope.source(eid));
    const vertex_data_type& vdata = 
      scope.const_neighbor_vertex_data(vid);
    if(vdata.is_set) ++local_atom2count[vdata.atomid];
  }
  foreach(const edge_id_t eid, scope.out_edge_ids()) {
    const vertex_id_t vid(scope.target(eid));
    const vertex_data_type& vdata = 
      scope.const_neighbor_vertex_data(vid);
    if(vdata.is_set) ++local_atom2count[vdata.atomid];
  }

  // Get the vertex data
  const vertex_data_type& vdata(scope.const_vertex_data());


  // If the neighbor change has not reached this machine yet then
  // reschedule self
  if(!vdata.is_seed && local_atom2count.empty()) {
    callback.add_task(scope.vertex(), partition_update_function);
    return;
  }

  bool changed(false);
  if(!vdata.is_seed) {
    ASSERT_GT(local_atom2count.size(), 0);
    // Get the new atomid assignment for this vertex
    typedef shared_atom2count_type::const_ptr_type shared_ptr_type;
    shared_ptr_type global_a2c_ptr(global_atom2count.get_ptr());
    const procid_t new_atomid = 
      find_best_atom(local_atom2count, *global_a2c_ptr);    

    if(!vdata.is_set ||
       (vdata.num_changes < MAX_CHANGES && 
        vdata.atomid != new_atomid) ) {
      vertex_data_type& vdata(scope.vertex_data());
      vdata.atomid = new_atomid;
      vdata.is_set = true;
      vdata.num_changes++;
      changed = true;
    }
  } // end update assig
  // Reschedule the neighbors
  if(changed || vdata.is_seed) {
    // Schedule all in neighbors
    foreach(const edge_id_t eid, scope.in_edge_ids()) {
      const vertex_id_t vid(scope.source(eid));
      const vertex_data_type& vdata = 
        scope.const_neighbor_vertex_data(vid);
      if(vdata.num_changes < MAX_CHANGES) 
        callback.add_task(vid, partition_update_function);
    }
    // Schedule all out neighbors
    foreach(const edge_id_t eid, scope.out_edge_ids()) {
      const vertex_id_t vid(scope.target(eid));
      const vertex_data_type& vdata = 
        scope.const_neighbor_vertex_data(vid);
      if(vdata.num_changes < MAX_CHANGES) 
        callback.add_task(vid, partition_update_function);
    }
  }
} // end of partition_update_function






int main(int argc, char** argv) {
  // set the global logger
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  // Initialize the mpi tools
  graphlab::mpi_tools::init(argc, argv);

  // Parse the command lines
  std::string aindex("atom_index.txt");
  std::string partfile("partitioning.txt");
  graphlab::command_line_options 
    clopts("Partition the graph using the GraphLab engine.");
  clopts.attach_option("aindex", &aindex, aindex,
                       "The atom index file.");
  clopts.attach_option("nparts", &num_atoms, num_atoms,
                       "The number of parts to create.");
  clopts.attach_option("partfile", &partfile, partfile,
                       "[output] file containing the partitioning.");
  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }
  logstream(LOG_INFO) << "Partitioning into " << num_atoms
                      << " parts." << std::endl;

  // Initialize the distributed control plane
  dc_init_param param;
  if( ! init_param_from_mpi(param) ) 
    logstream(LOG_FATAL) 
      << "Failed MPI laucher!" << std::endl; 
  param.initstring = "buffered_queued_send=yes, ";
  param.numhandlerthreads = 8;
  global_logger().set_log_level(LOG_DEBUG);
  distributed_control dc(param);


  logstream(LOG_INFO) 
    << "Loading graph from atom index file: " << aindex << std::endl;
  const bool NO_LOAD_DATA(true);
  graph_type  graph(dc, aindex, NO_LOAD_DATA);
 

  logstream(LOG_INFO)
    << "Artificially color the graph" << std::endl;
  foreach(const vertex_id_t vid, graph.owned_vertices()) {
    graph.color(vid) =  rand() % NUM_COLORS;
  }


  logstream(LOG_INFO)  
    << "Initializing engine with " << clopts.get_ncpus() 
    << " local threads." <<std::endl;
  engine_type engine(dc, graph, clopts.get_ncpus());


  logstream(LOG_INFO)  
    << "Set the scheduler options." << std::endl;  
  scheduler_options schedopts;
  schedopts.add_option("update_function", partition_update_function);
  //  schedopts.add_option("max_iterations", MAX_ITERATIONS);
  engine.set_scheduler_options(schedopts);

  logstream(LOG_INFO) << "Register a sync." << std::endl;
  engine.set_sync(global_atom2count,
                  atom2count_sum_fun,
                  atom2count_apply_fun,
                  any(atom2count_type()), 
                  SYNC_INTERVAL,
                  atom2count_merge_fun);

  logstream(LOG_INFO) << "Scheduling tasks." << std::endl;  

  // Scheduling tasks
  if(dc.procid() == 0) {
    for(size_t i = 0; i < num_atoms; ++i) {
      const vertex_id_t vid(rand() % graph.num_vertices());
      vertex_data_type vdata;
      vdata.atomid = i;
      vdata.is_set = true;  
      vdata.is_seed = true;
      graph.set_vertex_data(vid, vdata);
      logstream(LOG_INFO) << "Adding seed: " << vid;
      engine.add_vtask(vid, partition_update_function);
    }
  }

  logstream(LOG_INFO) << "Running partitioner." << std::endl;
  engine.start();
  logstream(LOG_INFO) << "Finished partitioning." << std::endl;

  logstream(LOG_INFO) << "Gathering partitioning." << std::endl;

  typedef std::vector< std::pair<vertex_id_t, procid_t> > vector_of_pairs;
  std::vector<vector_of_pairs> proc2pairs(dc.numprocs());
  foreach(const vertex_id_t vid, graph.owned_vertices()) {
    const vertex_data_type& vdata(graph.vertex_data(vid));
    // Require all vertices to be assinged a class
    ASSERT_TRUE(vdata.is_set);
    proc2pairs[dc.procid()].
      push_back(std::make_pair(vid, vdata.atomid));
  }
  const size_t ROOT_NODE(0);
  dc.gather(proc2pairs, ROOT_NODE);
  if (dc.procid() == ROOT_NODE) {
    // construct final map
    std::vector<procid_t> result(graph.num_vertices());
    std::vector<size_t> counts(num_atoms);
    std::vector<size_t> vertex2proc(graph.num_vertices());
    for (size_t i = 0; i < dc.numprocs(); ++i) {
      for(size_t j = 0; j < proc2pairs[i].size(); ++j) {
        result.at(proc2pairs[i][j].first) = proc2pairs[i][j].second;
        counts.at(proc2pairs[i][j].second)++;
        vertex2proc.at(proc2pairs[i][j].first) = i;
      }
    }
    {
      std::ofstream fout(partfile.c_str());
      ASSERT_TRUE(fout.good());
      for(size_t i = 0; i < result.size(); ++i) 
        fout << result[i] << "\n";    
      fout.close();
    }
    {
      std::string fname = "machine_" + partfile;
      std::ofstream fout(fname.c_str());
      ASSERT_TRUE(fout.good());
      for(size_t i = 0; i < result.size(); ++i) 
        fout << vertex2proc[i] << "\n";    
      fout.close();
    }

    
    std::cout << "\n\n\n\n" << std::endl 
              <<  "======================================"
              << "\n\n" << std::endl;

    std::cout << "Counts:  ";
    size_t max_counts(0);
    for(size_t i = 0; i < counts.size(); ++i) {
      std::cout << counts[i]  << '\t';
      max_counts = std::max(max_counts, counts[i]);
    }
    std::cout << std::endl;

    std::cout << "ECounts: ";
    atom2count_type atom2counts(global_atom2count.get_val());
    typedef atom2count_type::value_type pair_type;
    foreach(pair_type pair, atom2counts) 
      std::cout << pair.second  << '\t';
    std::cout << std::endl;

    const double imbalance = 
      double(max_counts) * double(counts.size()) / 
      double(graph.num_vertices());
    std::cout << "Imbalance max/average: " << imbalance << std::endl;

    std::cout << "\n\n" << std::endl 
              <<  "======================================"
              << "\n\n\n\n" << std::endl;




  }
 
  // Gather metrics
  dc.fill_metrics();
  graph.fill_metrics();

  
  if(dc.procid() == ROOT_NODE) {
    basic_reporter reporter;
    metrics::report_all(reporter);
  } 
  
  logstream(LOG_INFO) << "Finished " << dc.procid() << std::endl;
  dc.full_barrier();
  graphlab::mpi_tools::finalize();  
  return EXIT_SUCCESS;
}

