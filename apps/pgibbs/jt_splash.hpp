#ifndef PGIBBS_JT_SPLASH_HPP
#define PGIBBS_JT_SPLASH_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <string>
#include <cassert>
#include <algorithm>

#include <boost/unordered_set.hpp>


// Including Standard Libraries

#include <graphlab.hpp>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>


#include "factorized_model.hpp"
#include "mrf.hpp"
#include "junction_tree.hpp"






#include <graphlab/macros_def.hpp>

class jt_worker : public graphlab::runnable {
public:
  
  typedef graphlab::general_scope_factory<mrf::graph_type>
  scope_factory_type;
  typedef scope_factory_type::iscope_type iscope_type;
  
public:
  
  scope_factory_type* scope_factory;
  size_t worker_id;
  size_t worker_count;
  size_t max_tree_size;
  size_t max_tree_width;
  size_t max_factor_size;
  size_t max_tree_height;

  size_t tree_count;

  size_t total_samples;
  size_t collisions;
  size_t changes;

  float finish_time_seconds;

  bool use_priorities;

  std::vector<bool> assigned_factors;


  const factorized_model::factor_map_t* factors_ptr;

  // Tree building data structures 
  size_t root_index;
  const std::vector<vertex_id_t>* roots;
  vertex_id_t current_root;

  elim_map_t elim_time_map;
  clique_vector cliques;
  std::deque<vertex_id_t> bfs_queue;
  graphlab::mutable_queue<size_t, double> priority_queue;
  boost::unordered_set<vertex_id_t> visited;

  // Local junction tree graphlab core
  junction_tree::gl::core jt_core;

 
  jt_worker() : 
    scope_factory(NULL), 
    worker_id(0),
    worker_count(0),
    max_tree_size(0),
    max_tree_width(0),
    max_factor_size(0),
    tree_count(0),
    total_samples(0),
    collisions(0),
    changes(0),
    finish_time_seconds(0),
    use_priorities(false) { }


  void init(size_t wid,
            scope_factory_type& sf, 
            const factorized_model::factor_map_t& factors,
            const std::vector<vertex_id_t>& root_perm, 
            size_t ncpus,         
            size_t treesize,
            size_t treewidth,
            size_t factorsize,
            size_t max_height,
            size_t internal_threads,
            bool priorities);

  
  //! set the runtime for this worker
  void set_runtime(float runtime_seconds) {
    assert(runtime_seconds >= 0);
    finish_time_seconds = 
      graphlab::lowres_time_seconds() + runtime_seconds;          
  }


  //! Get the next root
  void move_to_next_root() {  
    root_index += worker_count;
    if(root_index >= roots->size()) root_index = worker_id;
    current_root = roots->at(root_index);
  }


  //! The main loop
  void run();



  /**
   * Grab this vertex into the tree owned by worker id
   */
  bool quick_try_vertex(vertex_id_t vid);



  /**
   * Grab this vertex into the tree owned by worker id
   */
  bool try_grab_vertex(iscope_type& scope);


  /**
   * Release the vertex
   */
  void release_vertex(iscope_type& scope);




  /**
   * This function computes the value of adding the vertex to the tree
   *
   */
  factor_t clique_factor;
  factor_t product_of_marginals_factor;
  factor_t conditional_factor;
  factor_t marginal_factor;

  double score_vertex(vertex_id_t vid);

  double score_vertex_l1_diff(vertex_id_t vid);

  double score_vertex_log_odds(vertex_id_t vid);

  double score_vertex_lik(vertex_id_t vid);
  
  void grow_bfs_jtree();

  void grow_prioritized_jtree();


  size_t sample_once() {
    if(use_priorities) {
      // grow the prioritized junction tree data structure
      grow_prioritized_jtree();
    } else {
      // grow the bfs junction tree data structure
      grow_bfs_jtree();
    }


    assert(scope_factory != NULL);
    // Get the scope factory
    mrf::graph_type& mrf = scope_factory->get_graph();
    
    // If we failed to build a tree return failure
    if(cliques.empty()) return 0;

    //  std::cout << "Varcount: " << cliques.size() << std::endl;  

        // ///////////////////////////////////
        // // plot the graph
        // if(worker_id == 0) {
        //   std::cout << "Saving treeImage:" << std::endl;
        //   size_t rows = std::sqrt(mrf.num_vertices());
        //   image img(rows, rows);
        //   for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {
        //     vertex_id_t tree_id = mrf.vertex_data(vid).tree_id;
        //     img.pixel(vid) = 
        //         tree_id == vertex_id_t(-1)? 0 : tree_id + worker_count;
        //   }
        //   img.save(make_filename("tree", ".pgm", tree_count).c_str());
        // }



    // Build the junction tree and sample
    jt_core.graph().clear();
    jtree_from_cliques(mrf,
                       cliques,
		       assigned_factors,
                       elim_time_map,
                       jt_core.graph());

    // jtree_from_cliques(mrf,  
    //                    cliques.begin(), cliques.end(), 
    //                    jt_core.graph());

    // Rebuild the engine (clear the old scheduler)
    jt_core.rebuild_engine();
    // add tasks to all vertices
    jt_core.add_task_to_all(junction_tree::calibrate_update, 1.0);

    // Run the core
    jt_core.start();



    // Check that the junction tree is sampled

    size_t actual_tree_width = 0;
    size_t local_changes = 0;
    for(vertex_id_t vid = 0; 
        vid < jt_core.graph().num_vertices(); ++vid) {
      const junction_tree::vertex_data& vdata = 
        jt_core.graph().vertex_data(vid);
      assert(vdata.sampled);
      assert(vdata.calibrated);
      actual_tree_width = 
        std::max(vdata.variables.num_vars(), actual_tree_width); 
      local_changes += vdata.changes;
    } 
    changes += local_changes;
    
    // std::cout << "Treewidth: " << actual_tree_width << std::endl;
    // std::cout << "Local Changes: " << local_changes << std::endl;
      
    // Return the number of variables in the tree
    return elim_time_map.size();
  } // end of sample once



};  // End of JT work



class parallel_sampler {

  std::vector<jt_worker> workers;
  graphlab::general_scope_factory<mrf::graph_type> scope_factory;
  std::vector< vertex_id_t > roots;
  bool use_cpu_affinity;

public:

  parallel_sampler(const factorized_model& fmodel,
                   mrf::graph_type& mrf,
                   const graphlab::engine_options& eopts,
                   size_t max_tree_size = 1000,
                   size_t max_tree_width = MAX_DIM,
                   size_t max_factor_size = (1 << MAX_DIM),
                   size_t max_tree_height = 0,
                   size_t internal_threads = 1,
                   bool use_priorities = false) :
    workers(eopts.ncpus),
    scope_factory(mrf, eopts.ncpus, 
                  graphlab::scope_range::EDGE_CONSISTENCY),
    roots(mrf.num_vertices()),
    use_cpu_affinity(eopts.enable_cpu_affinities) { 

    // Shuffle ther oot ordering 
    for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid)
      roots[vid] = vid;
    std::random_shuffle(roots.begin(), roots.end());
       
    for(size_t i = 0; i < eopts.ncpus; ++i) {
      // Initialize the worker
      workers[i].init(i, 
                      scope_factory, 
                      fmodel.factors(),
                      roots,
                      eopts.ncpus,    
                      max_tree_size,
                      max_tree_width,
                      max_factor_size,
                      max_tree_height,
                      internal_threads,
                      use_priorities);    
    }


    
  } // end of constructor


  size_t total_changes() const {
    size_t total_changes = 0;
    // Record the total number of samples
    foreach(const jt_worker& worker, workers)      
      total_changes += worker.num_changes(); 
    return total_changes;
  }


  size_t total_samples() const {
    size_t total_samples = 0;
    foreach(const jt_worker& worker, workers) 
      total_samples += worker.num_samples();    
    return total_samples;
  }


  size_t total_collisions() const {
    size_t total_collisions = 0;
    foreach(const jt_worker& worker, workers) 
      total_collisions += worker.num_collisions();
    return total_collisions;
  }




  size_t total_trees() const {
    size_t total_trees = 0;
    foreach(const jt_worker& worker, workers) 
      total_trees += worker.num_trees();
    return total_trees;
  }





  
  void sample_once(float runtime_secs) {
    // create workers
    graphlab::thread_group threads;
    
    for(size_t i = 0; i < workers.size(); ++i) {
      workers[i].set_runtime(runtime_secs);
      // Launch the threads
      if(use_cpu_affinity) threads.launch(&(workers[i]), i);
      else threads.launch(&(workers[i]));            
    }
 
    // Wait for all threads to finish
    threads.join();

  }                   

 

};
















#include <graphlab/macros_undef.hpp>
#endif
