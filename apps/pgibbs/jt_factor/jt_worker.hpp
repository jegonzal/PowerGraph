#ifndef JT_WORKER_HPP
#define JT_WORKER_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <string>
#include <cassert>
#include <algorithm>



// Including Standard Libraries

#include <graphlab.hpp>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>

#include "data_structures.hpp"


#include "sequential_jt_gibbs.hpp"



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

  float finish_time_seconds;

  bool use_priorities;



  const factorized_model::factor_map_t* factors_ptr;

  // Tree building data structures 
  size_t root_index;
  const std::vector<vertex_id_t>* roots;
  vertex_id_t current_root;

  std::map<vertex_id_t, vertex_id_t> elim_time_map;
  clique_vector cliques;
  std::deque<vertex_id_t> bfs_queue;
  graphlab::mutable_queue<size_t, double> priority_queue;
  std::set<vertex_id_t> visited;

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
            bool priorities) {
    // Initialize parameters
    scope_factory = &sf;
    worker_id = wid;
    worker_count = ncpus;
    max_tree_size = treesize;
    max_tree_width = treewidth;
    max_tree_height = max_height;


    if(factorsize <= 0) {
      max_factor_size = std::numeric_limits<size_t>::max();
    } else {
      max_factor_size = factorsize;
    }
    factors_ptr = &factors;
    use_priorities = priorities;

    roots = &root_perm;    
    root_index = root_perm.size();
    current_root = worker_id;

    // Initialize local jtcore
    if(internal_threads > 1) {
      jt_core.set_scheduler_type("fifo");
      jt_core.set_scope_type("edge");
      jt_core.set_ncpus(internal_threads);
      jt_core.set_engine_type("async");
    } else {
      jt_core.set_scheduler_type("fifo");
      jt_core.set_scope_type("none");
      jt_core.set_ncpus(1);
      jt_core.set_engine_type("async_sim");
    }


    // Initialize the shared data in the factorized model
    jt_core.shared_data().set_constant(junction_tree::FACTOR_KEY, 
                                       factors_ptr);
    mrf::graph_type* mrf_graph_ptr = &scope_factory->get_graph();
    jt_core.shared_data().set_constant(junction_tree::MRF_KEY, 
                                       mrf_graph_ptr); 
  }


  size_t num_samples() const { return total_samples; }
  size_t num_collisions() const { return collisions; }
  size_t num_trees() const { return tree_count; }

  void set_runtime(float runtime_seconds) {
    assert(runtime_seconds >= 0);
    finish_time_seconds = 
      graphlab::lowres_time_seconds() + runtime_seconds;          
  }


  
  void move_to_next_root() {
    // current_root = 
    //   graphlab::random::rand_int(scope_factory->num_vertices() - 1);
    root_index += worker_count;
    if(root_index >= roots->size()) root_index = worker_id;
    current_root = roots->at(root_index);
  }


  // get a root
  void run() {   
    // looup until runtime is reached
    while(graphlab::lowres_time_seconds() < finish_time_seconds) {
      /////////////////////////////////////////////////////////
      // Construct one tree (we must succeed in order to count a tree
      size_t sampled_variables = 0;
      move_to_next_root();
      while(sampled_variables == 0 && 
            graphlab::lowres_time_seconds() < finish_time_seconds) {
        //  move_to_next_root();
        sampled_variables = sample_once();
        if(sampled_variables == 0) {
          collisions++;
          // sched_yield();
        }
      }


//       // Get a local copy of the graph
//       mrf::graph_type& mrf(scope_factory->get_graph());
//       if(worker_id == 0) {
//         std::cout << "Saving sample: " << std::endl;
//         size_t rows = std::sqrt(mrf.num_vertices());
//         image img(rows, rows);
//         for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {   
//           img.pixel(vid) = mrf.vertex_data(vid).asg.asg_at(0);
//         }
//         img.pixel(0) = 0;
//         img.pixel(1) = mrf.vertex_data(0).variable.arity -1;
//         img.save(make_filename("sample", ".pgm", tree_count).c_str());
//       }

      tree_count++;
      total_samples += sampled_variables;

    } 
  }

  /**
   * Grab this vertex into the tree owned by worker id
   */
  bool try_grab_vertex(iscope_type& scope) {
    // Check that this vertex is not already in a tree
    bool in_tree = scope.vertex_data().tree_id != NULL_VID;
    if(in_tree) return false;

    // check that the neighbors are not in any other trees than this
    // one
    foreach(edge_id_t in_eid, scope.in_edge_ids()) {
      vertex_id_t neighbor_vid = scope.source(in_eid);
      const mrf::vertex_data& vdata = 
        scope.const_neighbor_vertex_data(neighbor_vid);
      bool in_tree = vdata.tree_id != NULL_VID;
      // if the neighbor is in a tree other than this one quit
      if(in_tree && worker_id != vdata.tree_id) return false;
    }
    // Assert that this vertex is not in a tree and that none of the
    // neighbors are in other trees
    // This vertex does not neighbor any other trees than this one
    scope.vertex_data().tree_id = worker_id;
    return true;
  } // end of try grab vertex


  /**
   * Release the vertex
   */
  void release_vertex(iscope_type& scope) {
    // This vertex does not neighbor any other trees than this one
    scope.vertex_data().tree_id = NULL_VID;
  } // release the vertex



  factor_t clique_factor;
  factor_t conditional_factor;
  factor_t marginal_factor;
  /**
   * This function computes the value of adding the vertex to the tree
   *
   */
  double score_vertex(vertex_id_t vid) {
    // Get the scope factory
    const mrf::graph_type& mrf = scope_factory->get_graph();
    const mrf::vertex_data& vdata = mrf.vertex_data(vid);

    // Construct the domain of neighbors that are already in the tree
    domain_t in_tree_vars = vdata.variable;
    foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
      const vertex_id_t neighbor_vid = mrf.source(ineid);
      const mrf::vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
      // test to see if the neighbor is in the tree by checking the
      // elimination time map
      if(elim_time_map.find(neighbor_vid) != elim_time_map.end()) {
        // otherwise add the tree variable
        in_tree_vars += neighbor.variable;
        // If this vertex has too many tree neighbor than the priority
        // is set to 0;
        if(in_tree_vars.num_vars() > max_tree_width) return 0;
        if(in_tree_vars.size() > max_factor_size) return 0;
      } 
    }

    
    // Compute the clique factor
    clique_factor.set_args(in_tree_vars);
    clique_factor.uniform();
    // get all the factors
    const factorized_model::factor_map_t& factors(*factors_ptr);
    // Iterate over the factors and multiply each into this factor
    foreach(size_t factor_id, vdata.factor_ids) {
      const factor_t& factor = factors[factor_id];      
      // Build up an assignment for the conditional
      domain_t conditional_args = factor.args() - in_tree_vars;
      if(conditional_args.num_vars() > 0) {
        assignment_t conditional_asg;
        for(size_t i = 0; i < conditional_args.num_vars(); ++i)
          conditional_asg &= mrf.vertex_data(conditional_args.var(i).id).asg;      
        // set the factor arguments
        conditional_factor.set_args(factor.args() - conditional_args);
        conditional_factor.condition(factor, conditional_asg);        
        // Multiply the conditional factor in
        clique_factor *= conditional_factor;
        //        clique_factor.normalize();
      } else {
        clique_factor *= factor;
      }
    } // end of loop over factors
    // Compute the conditional factor and marginal factors
    conditional_factor.set_args(in_tree_vars - vdata.variable);
    conditional_factor.condition(clique_factor, vdata.asg);        
    marginal_factor.set_args(in_tree_vars - vdata.variable);
    marginal_factor.marginalize(clique_factor);
    
    // Compute metric
    conditional_factor.normalize();
    marginal_factor.normalize();
    double residual = conditional_factor.log_residual(marginal_factor);
    
    residual = (std::tanh(residual) + 1) / (vdata.updates + 1);

    // rescale by updates
    //    residual = residual / (vdata.updates + 1);

    assert( residual >= 0);
    assert( !std::isnan(residual) );
    assert( std::isfinite(residual) );

    return residual;
  } // end of score vertex

  

  void grow_bfs_jtree() {
    assert(scope_factory != NULL);
    // Get the scope factory
    mrf::graph_type& mrf = scope_factory->get_graph();
    // Clear local data structures
    elim_time_map.clear();
    cliques.clear();
    bfs_queue.clear();
    visited.clear();
     
    // add the root
    bfs_queue.push_back(current_root);
    visited.insert(current_root);

    while(!bfs_queue.empty()) {
      // Take the top element
      const vertex_id_t next_vertex = bfs_queue.front();
      bfs_queue.pop_front();

      // Get the scope
      iscope_type* scope_ptr = 
        scope_factory->get_edge_scope(worker_id, next_vertex);
      assert(scope_ptr != NULL);
      iscope_type& scope(*scope_ptr);

      // See if we can get the vertex for this tree
      bool grabbed = try_grab_vertex(scope);

      // If we failed to grab the scope then skip this vertex
      if(grabbed) {

        // if this is the root vertex
        vertex_id_t min_height = 0;
        if(max_tree_height != 0 && !cliques.empty()) {
          min_height = max_tree_height;
          // find the closest vertex to the root
          foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
            vertex_id_t neighbor_vid = mrf.target(eid);
            // if the neighbor is already in the tree
            if(elim_time_map.find(neighbor_vid) != elim_time_map.end()) {
              min_height = 
                std::min(min_height, mrf.vertex_data(neighbor_vid).height + 1);
            } 
          }
        } 
        bool safe_extension = 
          (max_tree_height == 0) ||
          (min_height < max_tree_height);

        // test that its safe  
        safe_extension = safe_extension && 
          extend_clique_list(mrf,
                             next_vertex,
                             elim_time_map,
                             cliques,
                             max_tree_width,
                             max_factor_size);

        if(safe_extension) {   
          // set the height
          mrf.vertex_data(next_vertex).height = min_height;

          // add the neighbors to the search queue
          foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
            vertex_id_t neighbor_vid = mrf.target(eid);
            if(visited.count(neighbor_vid) == 0) {
              bfs_queue.push_back(neighbor_vid);
              visited.insert(neighbor_vid);
            }
          }
        } else {
          // release the vertex since it could not be used in the tree
          release_vertex(scope);
        }
      } // end of grabbed
      // release the scope
      scope_factory->release_scope(&scope);        
      // Limit the number of variables
      if(cliques.size() > max_tree_size) break;
    } // end of while loop
  } // end grow_bfs_jtree



  void grow_prioritized_jtree() {
    assert(scope_factory != NULL);
    // Get the scope factory
    mrf::graph_type& mrf = scope_factory->get_graph();
    // Clear local data structures
    elim_time_map.clear();
    cliques.clear();
    priority_queue.clear();
    visited.clear();
     
    // add the root
    priority_queue.push(current_root, 1.0);
    visited.insert(current_root);

    while(!priority_queue.empty()) {
      // Take the top element
      const vertex_id_t next_vertex = priority_queue.pop().first;
      // Get the scope
      iscope_type* scope_ptr = 
        scope_factory->get_edge_scope(worker_id, next_vertex);
      assert(scope_ptr != NULL);
      iscope_type& scope(*scope_ptr);

      // See if we can get the vertex for this tree
      bool grabbed = try_grab_vertex(scope);

      // If we failed to grab the scope then skip this vertex
      if(grabbed) {
        // compute the tree height of the new vertex
        vertex_id_t min_height = 0;        
        // if this is not the root and we care about tree height
        if(max_tree_height != 0 && !cliques.empty()) {
          min_height = max_tree_height;
          // find the closest vertex to the root
          foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
            vertex_id_t neighbor_vid = mrf.target(eid);
            // if the neighbor is already in the tree
            if(elim_time_map.find(neighbor_vid) != elim_time_map.end()) {
              min_height = 
                std::min(min_height, 
                         mrf.vertex_data(neighbor_vid).height + 1);
            } 
          }
        } // end of tree height check for non root vertex 
          
        // test the 
        bool safe_extension = 
          (max_tree_height == 0) ||
          (min_height < max_tree_height);

        safe_extension = safe_extension &&
          extend_clique_list(mrf,
                             next_vertex,
                             elim_time_map,
                             cliques,
                             max_tree_width,
                             max_factor_size);

        // If the extension was safe than the elim_time_map and
        // cliques data structure are automatically extended
        if(safe_extension) {
          // set the height
          mrf.vertex_data(next_vertex).height = min_height;

          // add the neighbors to the search queue or update their priority
          foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
            vertex_id_t neighbor_vid = mrf.target(eid);          
            if(visited.count(neighbor_vid) == 0) {
              // Vertex has not yet been visited
              double score = score_vertex(neighbor_vid);
              // if the score is greater than zero then add the
              // neighbor to the priority queue.  The score is zero if
              // there is no advantage or the treewidth is already too
              // large
              if(score > 0) priority_queue.push(neighbor_vid, score);
              visited.insert(neighbor_vid);

            } 
            // else if(priority_queue.contains(neighbor_vid)) {
            //   // vertex is still in queue we may need to recompute
            //   // score
            //   double score = score_vertex(neighbor_vid);
            //   if(score > 0) {
            //     // update the priority queue with the new score
            //     priority_queue.update(neighbor_vid, score);
            //   } else {
            //     // The score computation revealed that the clique
            //     // would be too large so simply remove the vertex from
            //     // the priority queue
            //     priority_queue.remove(neighbor_vid);
            //   }
            // } // otherwise the vertex has been visited and processed
          }
        } else {
          // release the vertex since it could not be used in the tree
          release_vertex(scope);
        }
      } // end of grabbed      
      // release the scope
      scope_factory->release_scope(&scope);        
      // Limit the number of variables
      if(cliques.size() > max_tree_size) break;
    } // end of while loop

  } // end grow_prioritized_jtree



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

    //        std::cout << "Varcount: " << cliques.size() << std::endl;  

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
                       elim_time_map,
                       cliques.begin(), cliques.end(), 
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
    for(vertex_id_t vid = 0; 
        vid < jt_core.graph().num_vertices(); ++vid) {
      const junction_tree::vertex_data& vdata = 
        jt_core.graph().vertex_data(vid);
      assert(vdata.sampled);
      assert(vdata.calibrated);
      actual_tree_width = 
        std::max(vdata.variables.num_vars(), actual_tree_width);     
    } 
    
    //    std::cout << "Treewidth: " << actual_tree_width << std::endl;

      
    // Sampled root successfully
    return cliques.size();
  } // end of sample once



};  // End of JT work



class parallel_sampler {

  std::vector<jt_worker> workers;
  graphlab::general_scope_factory<mrf::graph_type> scope_factory;
  std::vector< vertex_id_t > roots;
  size_t total_samples;
  size_t total_collisions;


public:

  parallel_sampler(const factorized_model& fmodel,
                   mrf::graph_type& mrf,
                   size_t ncpus,
                   size_t max_tree_size = 1000,
                   size_t max_tree_width = MAX_DIM,
                   size_t max_factor_size = (1 << MAX_DIM),
                   size_t max_tree_height = 1000,
                   size_t internal_threads = 1,
                   bool use_priorities = false) :
    workers(ncpus),
    scope_factory(mrf, ncpus, 
                  graphlab::scope_range::EDGE_CONSISTENCY),
    roots(mrf.num_vertices()),
    total_samples(0),
    total_collisions(0) { 

    // Shuffle ther oot ordering 
    for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid)
      roots[vid] = vid;
    std::random_shuffle(roots.begin(), roots.end());
       
    for(size_t i = 0; i < ncpus; ++i) {
      // Initialize the worker
      workers[i].init(i, 
                      scope_factory, 
                      fmodel.factors(),
                      roots,
                      ncpus,    
                      max_tree_size,
                      max_tree_width,
                      max_factor_size,
                      max_tree_height,
                      internal_threads,
                      use_priorities);    
    }


    
  } // end of constructor



  
  void sample_once(float runtime_secs) {
    // create workers
    graphlab::thread_group threads;
    
    for(size_t i = 0; i < workers.size(); ++i) {
      workers[i].set_runtime(runtime_secs);
      // Launch the threads
      bool use_cpu_affinity = false;
      if(use_cpu_affinity) threads.launch(&(workers[i]), i);
      else threads.launch(&(workers[i]));            
    }
 
    // Wait for all threads to finish
    threads.join();

    // Record the total number of samples
    foreach(const jt_worker& worker, workers) {
      total_samples += worker.num_samples();
      total_collisions += worker.num_collisions();
    }
    std::cout << "Total samples: " << total_samples << "\n";
    std::cout << "Total collisions: " << total_collisions << "\n";

  }                   

 

};
















#include <graphlab/macros_undef.hpp>
#endif
