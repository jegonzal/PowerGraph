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
const size_t MAX_CHANGES(100);
const size_t MAX_ITERATIONS(20);
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
  vertex_data_type() : 
    atomid(0), num_changes(0), 
    is_set(false), is_seed(false) { }
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


size_t NVERTS(-1);
size_t NEDGES(-1);
size_t NATOMS(-1);
// size_t MAX_SAMPLES(-1);

struct statistics {
  typedef std::vector<size_t> atom2count_type;

  
  // Counts of vertices assigned to different atoms
  atom2count_type atom2vcount;
  // counts of edges assigned to different atoms
  atom2count_type atom2ecount;

  size_t vset;
  size_t eset;
  size_t visited;
  size_t edges_cut;
  size_t nchanges;

  // Random sample of unvisited vertices
  std::vector<vertex_id_t> samples;  

  statistics(size_t natoms = 0) : 
    atom2vcount(natoms, 0), atom2ecount(natoms, 0),
    vset(0), eset(0), visited(0),
    edges_cut(0), nchanges(0) { }

  
  void operator+=(const iscope_type& iscope) {
    visited++;
    const vertex_data_type& vdata(iscope.const_vertex_data());
    nchanges += vdata.num_changes;
    if(vdata.is_set)  {
      ASSERT_LT(vdata.atomid, NATOMS);
      // update counts
      atom2vcount[vdata.atomid]++;
      atom2ecount[vdata.atomid] += iscope.in_edge_ids().size();
      // compute edges cut
      foreach(const edge_id_t eid, iscope.in_edge_ids()) {
        const vertex_id_t nvid(iscope.source(eid));
        const vertex_data_type& nvdata = 
          iscope.const_neighbor_vertex_data(nvid);
        if(nvdata.is_set && nvdata.atomid != vdata.atomid) edges_cut++;
      } // end of loop over edges
    } else {
      // Update uniform samples
      if(samples.size() < NATOMS) {
        samples.push_back(iscope.vertex());
      } else if( random::rand01() < double(samples.size())/visited) {
        samples[rand() % samples.size()] = iscope.vertex();
      }
    }
  }
  void operator+=(const statistics& other) {
    // Merge maps
    for(size_t i = 0; i < NATOMS; ++i) atom2vcount[i] += other.atom2vcount[i];
    for(size_t i = 0; i < NATOMS; ++i) atom2ecount[i] += other.atom2ecount[i];
    // Merge Samples into new samples
    std::vector<vertex_id_t> new_samples;
    // fill out the rest of this vector
    size_t i(0), j(0);
    while(new_samples.size() < NATOMS && 
          i < samples.size() && 
          j < other.samples.size()) {
      ASSERT_LT(j, other.visited);
      ASSERT_LT(i, visited);
      const double accept_prob = 
        double(other.visited - j) / 
        double(other.visited + visited - double(i+j));
      ASSERT_GE(accept_prob, 0);
      if(random::rand01() < accept_prob) 
        new_samples.push_back(other.samples[j++]);
      else new_samples.push_back(samples[i++]);      
    }
    while(new_samples.size() < NATOMS && i < samples.size()) 
      new_samples.push_back(samples[i++]);      
    while(new_samples.size() < NATOMS && j < other.samples.size()) 
      new_samples.push_back(other.samples[j++]);      
    samples.swap(new_samples);

    // Merge basic counters
    visited += other.visited;
    edges_cut += other.edges_cut;
    nchanges += other.nchanges;

    // check local data structure
    ASSERT_EQ(atom2vcount.size(), NATOMS);
    ASSERT_EQ(atom2ecount.size(), NATOMS);
  }
  void print() {
    std::cout 
      << "------------------------------------------------------------\n"
      << "Vset:         " << vset << '\n'
      << "Eset:         " << eset << '\n'
      << "Visited:      " << visited << '\n'
      << "Edges Cut:    " << edges_cut << '\n'
      << "Nchanges:     " << nchanges << '\n'
      << "VBalance:     " << vertex_imbalance() << '\n'
      << "EBalance:     " << edge_imbalance() << '\n'
      << "Counts ------ " << '\n';
    for(size_t i = 0; i < atom2vcount.size(); ++i) {
      const double VBAL(vertex_prop(i));
      const double EBAL(edge_prop(i));
      const size_t VSTARS(std::min(size_t(8), 
                                   size_t(std::ceil(VBAL*NATOMS))));
      const size_t ESTARS(std::min(size_t(8), 
                                   size_t(std::ceil(EBAL*NATOMS))));
      std::cout << std::right << std::setw(5) << i << "|"
                << std::right << std::setw(10) << atom2vcount[i] << "|"
                << std::right << std::setw(10) <<  VBAL << "|"
                << std::left  << std::setw(10) << std::string(VSTARS, '*')
                << "|"
                << std::right << std::setw(10) << atom2ecount[i] << "|"
                << std::right << std::setw(10) << EBAL << "|"
                << std::left  << std::setw(10) << std::string(ESTARS, '*')
                << '\n';
    }
    std::cout 
      << "------------------------------------------------------------"
      << std::endl;      
  }
  double vertex_prop(size_t i) const {
    ASSERT_LT(i, NATOMS);
    ASSERT_LE(vset, NVERTS);
    if(vset == 0) return 0;
    return double(atom2vcount[i]) / vset;
  }
  double edge_prop(size_t i) const {
    ASSERT_LT(i, NATOMS);
    ASSERT_LE(eset, NEDGES);
    if(eset == 0) return 0;
    return double(atom2ecount[i]) / eset;
  }

  size_t edge_min() const {
    return
      std::min_element(atom2ecount.begin(), atom2ecount.end()) -
      atom2ecount.begin();
  }
  size_t vertex_min() const {
    return
      std::min_element(atom2vcount.begin(), atom2vcount.end()) -
      atom2vcount.begin();
  }

  double vertex_imbalance() const {
    if(vset == 0) return 0;
    const size_t max_ind = 
      std::max_element(atom2vcount.begin(), atom2vcount.end()) -
      atom2vcount.begin();
    return double(atom2vcount[max_ind] * NATOMS) / double(vset);
  }


  double edge_imbalance() const {
    if(eset == 0) return 0;
    const size_t max_ind = 
      std::max_element(atom2ecount.begin(), atom2ecount.end()) -
      atom2ecount.begin();
    return double(atom2ecount[max_ind] * NATOMS) / double(eset);
  }



  void finalize() {
    vset = 0; 
    foreach(size_t count, atom2vcount) vset+=count;
    eset = 0;
    foreach(size_t count, atom2ecount) eset+=count;    
    ASSERT_LE(vset, NVERTS);
    ASSERT_LE(eset, NEDGES);
    ASSERT_EQ(visited, NVERTS);
    print(); 
  }   
  void load(iarchive& iarc) {
    iarc >> atom2vcount 
         >> atom2ecount 
         >> vset
         >> eset
         >> visited
         >> edges_cut
         >> nchanges
         >> samples;
  }
  void save(oarchive& oarc) const {
    oarc << atom2vcount 
         << atom2ecount 
         << vset
         << eset
         << visited
         << edges_cut
         << nchanges
         << samples;
  }
};

typedef distributed_glshared<statistics> shared_statistics_type;
// global data
shared_statistics_type shared_statistics;

// update the counts in the appropriate table
void statistics_sum_fun(iscope_type& iscope,  any& acc) {
  acc.as<statistics>() += iscope;
}
// Identity apply
void statistics_apply_fun(any& current_data, 
                           const any& acc) { 
  current_data.as<statistics>() = acc.as<statistics>();
  current_data.as<statistics>().finalize();
} 
// Sum the two maps
void statistics_merge_fun(any& any_dest,  const any& any_src) {
  any_dest.as<statistics>() += any_src.as<statistics>();
}


















void partition_update_function(iscope_type& scope,
                               icallback_type& callback,
                               ishared_data_type* unused) {
  statistics::atom2count_type nbr_a2c(NATOMS);
  // Get the number of neighbor assignments ===============================
  size_t nbr_sum(0);
  foreach(const edge_id_t eid, scope.in_edge_ids()) {
    const vertex_data_type& nbr_vdata = 
      scope.const_neighbor_vertex_data(scope.source(eid));
    if(nbr_vdata.is_set) { 
      ++nbr_a2c[nbr_vdata.atomid];
      ++nbr_sum;
    }
  }
  foreach(const edge_id_t eid, scope.out_edge_ids()) {
    const vertex_data_type& nbr_vdata = 
      scope.const_neighbor_vertex_data(scope.target(eid));
    if(nbr_vdata.is_set) { 
      ++nbr_a2c[nbr_vdata.atomid];
      ++nbr_sum;
    }
  }
  // Determine the new atomid to join (if any) ============================
  bool changed(false);
  // Get the vertex data (as a constant)
  const vertex_data_type& vdata(scope.const_vertex_data());
  if(!vdata.is_seed && !vdata.is_set) {
    // Get the shared statistics
    typedef shared_statistics_type::const_ptr_type shared_ptr_type;
    shared_ptr_type shared_statistics_ptr(shared_statistics.get_ptr());
    const statistics& stats(*shared_statistics_ptr);
    // Compute the best atomid to join
    double best_score(0);
    procid_t best_atomid(random::rand_int(NATOMS-1));
    // if(stats.eset == 0) best_atomid = rand() % NATOMS;
    // else best_atomid = stats.edge_min();    
    //    best_atomid = stats.edge_min();    
    for(size_t i = 0; i < NATOMS; ++i) {
      double new_imbalance(0);
      if(!vdata.is_set)
        new_imbalance = 
          double((stats.atom2ecount[i]+scope.in_edge_ids().size()) * NATOMS) 
          / double(stats.eset + 1 + scope.in_edge_ids().size());   
      else if(vdata.is_set && vdata.atomid == i) 
        new_imbalance = 
          double(stats.atom2ecount[i] * NATOMS) / double(stats.eset + 1);     
      else 
        new_imbalance = 
          double((stats.atom2ecount[i]+scope.in_edge_ids().size()) * NATOMS) 
          / double(stats.eset + 1);   


      //    const double score = nbr_a2c[i];
      // if(new_imbalance < 3.0 && score > best_score) {
      //   best_score = score; best_atomid = i;
      // }

      const double score = 
        (double(nbr_a2c[i]) / (nbr_sum + 1)) / 
        (double(stats.atom2ecount[i]) / (stats.eset + 1)) ;
      if(score > best_score) {
        best_score = score; best_atomid = i;
      }

      // const double score = nbr_a2c[i];
      // if(score > best_score) {
      //   best_score = score; best_atomid = i;
      // }

    }
    // if we failed try again 
    if(best_atomid == procid_t(-1)) {
      callback.add_task(scope.vertex(), partition_update_function); 
      return;
    }
    
    if(best_atomid != procid_t(-1) &&
       vdata.num_changes < MAX_CHANGES && 
       (vdata.atomid != best_atomid || !vdata.is_set ))  {
      vertex_data_type& vdata(scope.vertex_data());
      vdata.atomid = best_atomid;
      vdata.is_set = true;
      vdata.num_changes++;
      changed = true;
    }
  } // end update assign
  // Reschedule the neighbors
  if(changed || vdata.is_seed) {
    // Schedule all in neighbors
    foreach(const edge_id_t eid, scope.in_edge_ids()) {
      const vertex_id_t vid(scope.source(eid));
      const vertex_data_type& vdata = 
        scope.const_neighbor_vertex_data(vid);
      if(vdata.num_changes < MAX_CHANGES &&
         !vdata.is_seed)
        callback.add_task(vid, partition_update_function);
    }
    // Schedule all out neighbors
    foreach(const edge_id_t eid, scope.out_edge_ids()) {
      const vertex_id_t vid(scope.target(eid));
      const vertex_data_type& vdata = 
        scope.const_neighbor_vertex_data(vid);
      if(vdata.num_changes < MAX_CHANGES &&
         !vdata.is_seed)
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
  clopts.attach_option("nparts", &NATOMS, NATOMS,
                       "The number of parts to create.");
  clopts.attach_option("partfile", &partfile, partfile,
                       "[output] file containing the partitioning.");
  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }
  logstream(LOG_INFO) << "Partitioning into " << NATOMS
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
  
  // Initialize global quantites
  NVERTS = graph.num_vertices();
  NEDGES = graph.num_edges();

  shared_statistics.set(statistics(NATOMS));

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
  engine.set_sync(shared_statistics,
                  statistics_sum_fun,
                  statistics_apply_fun,
                  any(statistics(NATOMS)), 
                  SYNC_INTERVAL,
                  statistics_merge_fun);

  if(dc.procid() == 0) {
    std::cout << "Schedule initial set of vertices" << std::endl;
    for(size_t i = 0; i < NATOMS; ++i) {
      const vertex_id_t vid(rand() % graph.num_vertices());
      vertex_data_type vdata;
      vdata.atomid = i;
      vdata.is_set = true;
      vdata.is_seed = true;
      graph.set_vertex_data(vid, vdata);
      engine.add_vtask(vid, partition_update_function);
    }
  }


  // Scheduling tasks
  if(dc.procid() == 0)
    std::cout << "Running partitioner." << std::endl;
  bool finished(false);
  for(size_t iteration_counter = 0; 
        iteration_counter < MAX_ITERATIONS; iteration_counter++) {
    if(dc.procid() == 0)
      std::cout << "Starting iteration: " << iteration_counter
                << std::endl;        
    engine.start();
    statistics stats(shared_statistics.get_val());    
    if(dc.procid() == 0)
      std::cout << "Finished iteration: " << iteration_counter
              << std::endl;
    
    if(stats.vset == NVERTS){
      finished = true;
      break;
    }
    // Gather unset vertices
    if(dc.procid() == 0) stats.print();
    ASSERT_LE(stats.samples.size(), NATOMS);     
    // schedule local vertices
    foreach(const vertex_id_t vid, stats.samples) {
      ASSERT_NE(vid, vertex_id_t(-1));
      if(graph.vertex_is_local(vid)) {
        vertex_data_type vdata;
        vdata.atomid = rand() % NATOMS;
        vdata.is_set = true;
        vdata.is_seed = true;
        graph.set_vertex_data(vid, vdata);        
        engine.add_vtask(vid, partition_update_function);
      } // end of if 0
    } // end of loop over unset vertices
    dc.full_barrier();

  } // end of for iteration_counter

  if(!finished) {
    if(dc.procid() == 0) 
      std::cout  << "Running one last iteration. " << std::endl;
    // do one extra pass
    foreach(const vertex_id_t& vid, graph.owned_vertices()) {
      vertex_data_type& vdata = graph.vertex_data(vid);
      if(!vdata.is_set) {
        vdata.atomid = rand() % NATOMS;
        vdata.is_set = true;
        graph.set_vertex_data(vid, vdata);
        engine.add_vtask(vid, partition_update_function);       
      }
      //        vdata.is_seed = true;
      // graph.vertex_is_modified(vid);
      // graph.increment_vertex_version(vid);
    }
    if(dc.procid() == 0)
      std::cout  << "Running final run." << std::endl;
    engine.start();
  } // end of if finished
  

  if(dc.procid() == 0)
    std::cout  << "Gathering partitioning." << std::endl;

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
    std::vector<size_t> counts(NATOMS);
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
    statistics  stats(shared_statistics.get_val());
    stats.print();
    std::cout << "Totals:  ";
    size_t max_counts(0);
    for(size_t i = 0; i < counts.size(); ++i) {
      std::cout << counts[i]  << '\t';
      max_counts = std::max(max_counts, counts[i]);
    }
    std::cout << std::endl;

    std::cout << "VCounts: ";
    foreach(size_t count, stats.atom2vcount) 
      std::cout << count  << '\t';
    std::cout << std::endl;

    std::cout << "ECounts: ";
    typedef statistics::atom2count_type::value_type pair_type;
    foreach(size_t count, stats.atom2ecount) 
      std::cout << count  << '\t';
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

