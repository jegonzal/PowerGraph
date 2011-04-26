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


/**
 * The settings for the jt_splash_sampler.  Originally these formed a
 * long list of arguments but since the order can easily introduce
 * bugs we switched to a struct.
 */
struct splash_settings {
  size_t ntrees;
  size_t max_tree_size;
  size_t max_tree_height;
  size_t max_tree_width;
  size_t max_factor_size;
  bool   priorities;
  size_t vanish_updates; 
  size_t subthreads;

  splash_settings() : 
    ntrees(2), 
    max_tree_size(std::numeric_limits<size_t>::max()),
    max_tree_height(std::numeric_limits<size_t>::max()),
    max_tree_width(2),
    max_factor_size(std::numeric_limits<size_t>::max()),
    priorities(false),
    vanish_updates(10),
    subthreads(1) { }
};



/**
 * Run the jtsplash sampler
 *
 */
void run_jtsplash_sampler(mrf_graph_type& mrf_graph,
                          const std::string& jtsplash_results_fn,
                          const std::vector<double>& runtimes,
                          const bool draw_images,
                          const splash_settings& settings);







/**
 * This fairly complex update function assembles the clique factors by
 * conditioning on variables not in the tree.  Then it computes
 * messages at each clique to calibrate the junction tree.  Finally,
 * using the messages and the conditioned parents, it samples each
 * clique constructing new assignments to each variable.
 */
void jtree_sample_update(jtree_gl::iscope& scope,
                         jtree_gl::icallback& callback);



// Termination management
struct termination_condition {
  bool   error;
  float  finish_time_seconds;
  size_t target_nsamples;
  size_t target_ntrees;
  graphlab::atomic<size_t> atomic_nsamples;
  graphlab::atomic<size_t> atomic_ntrees;
  termination_condition();
  bool finished() const;
  void reset();
};







//! Predecleration 
class jt_worker;

/**
 * The jt_splash_sampler implements the junction tree based Gibbs
 * sampler defined in:
 *
 *  Parallel Gibbs Sampling: From Colored Fields to Think Junction Trees
 *   by Joseph Gonzalez, Yucheng Low, Arthur Gretton, and Carlos Guestrin
 *  
 */
class jt_splash_sampler {
public:
  typedef graphlab::general_scope_factory<mrf_graph_type>
  scope_factory_type;

private:
  std::vector<jt_worker*>     workers;
  scope_factory_type          scope_factory;
  std::vector< vertex_id_t >  root_perm;
  termination_condition       terminator;
public:
  jt_splash_sampler(mrf_graph_type& mrf_core,
                    const splash_settings& settings);
  ~jt_splash_sampler();


  /**
   * Get the number of times the splash sampler collided on a root.
   * This minor race event can lead to wasted cpu cycles but does not
   * affect the quality of the samples.
   */
  size_t total_collisions() const;

  /**
   * Get the total number of trees constructed on the last run
   */
  size_t total_trees() const;

  /**
   * Get the total number of single variable updates on the last run.
   */
  size_t total_samples() const;

  /** Run the splash sampler for a fixed number of seconds */
  void sample_seconds(float runtime_secs);
  /** Run the splash sampler for a fixed number of trees */
  void sample_trees(size_t total_trees);
  /**  
   * Run the splash sampler for a fixed number of single variable *
   * updates
   */
  void sample_updates(size_t total_updates);
private:
  void run();
};








//! The jt worker executes splashes sequential within each thread.
class jt_worker {
public:
  //! The scope factory type
  typedef jt_splash_sampler::scope_factory_type 
  scope_factory_type;
  typedef scope_factory_type::iscope_type iscope_type;
  
public:
  size_t worker_id;
  splash_settings settings;
  scope_factory_type* scope_factory_ptr;
  // Tree building data structures 
  size_t root_index;
  const std::vector<vertex_id_t>* root_perm_ptr;
  vertex_id_t current_root;
  //! track termination
  termination_condition* terminator_ptr;
  //! Track the collisions with the roots
  size_t ncollisions;
  //! Local junction tree graphlab core
  jtree_gl::core jt_core;

  /**
   * Local jt list used to build on the structure of the
   * jt_core.graph()
   */
  jtree_list jt_list;


  jt_worker(size_t worker_id, 
            const splash_settings& spsettings,
            scope_factory_type& sf, 
            const std::vector<vertex_id_t>& root_perm,
            termination_condition& terminator);
private:
  //! jt_worker is not copyable
  jt_worker(const jt_worker& other);
  //! jt_worker is not copyable
  jt_worker& operator=(const jt_worker& other);
public:

  //! The main loop
  void run();
  
private:
  //! Construct a single splash
  size_t splash_once();
  //! advance the root
  void advance_root();
  /**
   * Grab this vertex into the tree owned by worker id
   */
  bool is_vertex_available(vertex_id_t vid);
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
};  // End of JT worker




#endif
