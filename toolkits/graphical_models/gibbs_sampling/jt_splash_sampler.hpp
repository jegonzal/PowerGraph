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


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>

#include <graphlab.hpp>

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
class jtree_update :
  public graphlab::iupdate_functor<jtree_graph_type, jtree_update> {
public:  
  typedef graphlab::iupdate_functor<jtree_graph_type, jtree_update> base;
  jtree_update(mrf_graph_type* mrf_ptr = NULL) : mrf_ptr(mrf_ptr) { }
  mrf_graph_type* mrf_ptr;
  void operator()(base::icontext_type& context);
}; // end of class jtree_update



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
//! The jt worker executes splashes sequential within each thread.
class jt_builder : 
  public graphlab::iupdate_functor<mrf_graph_type, jt_builder>{
public:
  typedef graphlab::iupdate_functor<mrf_graph_type, jt_builder> base;

  struct splash_state {
    size_t worker_id;
    splash_settings settings;
    // Tree building data structures 
    size_t root_index;
    const std::vector<vertex_id_t>* root_perm_ptr;
    vertex_id_t current_root;
    //! track termination
    termination_condition* terminator_ptr;
    mrf_graph_type* graph_ptr;
    //! Track the collisions with the roots
    size_t ncollisions;
    //! Local junction tree graphlab core
    jtree_gl::core jt_core;
    /**
     * Local jt list used to build on the structure of the
     * jt_core.graph()
     */
    jtree_list jt_list;
    /**
     * Local data structures to reduce thread contention
     */
    std::deque<vertex_id_t> bfs_queue;
    graphlab::mutable_queue<size_t, double> priority_queue;
    boost::unordered_set<vertex_id_t> visited;
    factor_t clique_factor;
    factor_t product_of_marginals_factor;
    factor_t conditional_factor;
    factor_t marginal_factor;
  };

  std::set<splash_state*> state_set;
  jt_worker(splash_state* state_ptr = NULL);
  void operator+=(const jt_builder& other);






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


  double score_vertex(vertex_id_t vid);
  double score_vertex_l1_diff(vertex_id_t vid);
  double score_vertex_log_odds(vertex_id_t vid);
  double score_vertex_lik(vertex_id_t vid);

  void grow_bfs_jtree();
  void grow_prioritized_jtree();
};  // End of JT worker





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












#endif
