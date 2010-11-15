#ifndef PGIBBS_JT_SPLASH_SAMPLER_HPP
#define PGIBBS_JT_SPLASH_SAMPLER_HPP


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



//! used to calibrate and sample a junctiont ree
void jtree_sample_update(jtree_gl::iscope& scope,
                      jtree_gl::icallback& callback,
                      jtree_gl::ishared_data* shared_data);





//! The jt worker executes splashes sequential within each thread.
class jt_worker : public graphlab::runnable {
public:
  //! The scope factory type
  typedef graphlab::general_scope_factory<mrf_graph_type>
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
  bool use_priorities;

  size_t ntrees;
  size_t ncollisions;

  float finish_time_seconds;


  // Tree building data structures 
  size_t root_index;
  const std::vector<vertex_id_t>* roots;
  vertex_id_t current_root;


  //! Local junction tree graphlab core
  jtree_gl::core jt_core;

  /**
   *  Local jt list used to build on the structure of the
   * jt_core.grap()
   */
  jtree_list jt_list;


  jt_worker() : 
    scope_factory(NULL), 
    worker_id(0),
    worker_count(0),
    max_tree_size(0),
    max_tree_width(0),
    max_factor_size(0),
    max_tree_height(0),
    use_priorities(false),
    ntrees(0),
    ncollisions(0),
    finish_time_seconds(0), 
    root_index(0),
    roots(NULL),
    current_root(0){ }


  void init(size_t wid,
            scope_factory_type& sf, 
            const factorized_model::factor_map_t& factors,
            const std::vector<vertex_id_t>& root_perm, 
            size_t ncpus,         
            size_t treesize,
            size_t treewidth,
            size_t factorsize,
            size_t max_height,
            bool priorities,
            size_t internal_threads);

  
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
  

  /**
   * Local data structures to reduce thread contention
   */
  std::deque<vertex_id_t> bfs_queue;
  graphlab::mutable_queue<size_t, double> priority_queue;
  boost::unordered_set<vertex_id_t> visited;

  
  void grow_bfs_jtree();
  void grow_prioritized_jtree();
  size_t sample_once();


};  // End of JT worker















class jt_splash_sampler {

  std::vector<jt_worker> workers;
  jt_worker::scope_factory_type scope_factory;
  std::vector< vertex_id_t > roots;
  bool use_cpu_affinity;

public:

  jt_splash_sampler(mrf_gl::core& mrf_core,
                    const graphlab::engine_options& eopts,
                    size_t max_tree_size,
                    size_t max_tree_width,
                    size_t max_factor_size,
                    size_t max_tree_height, 
                    size_t internal_threads,
                    bool use_priorities);

  size_t total_collisions() const;
  size_t total_trees() const;

  
  void sample_once(float runtime_secs);
 
};







#endif
