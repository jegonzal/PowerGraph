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


 

private:
  
  scope_factory_type* scope_factory;
  size_t worker_id;
  size_t worker_count;
  size_t max_tree_size;
  size_t total_samples;
  size_t collisions;

  float finish_time_seconds;

  bool use_priorities;


  const factorized_model::factor_map_t* factors_ptr;

  // Tree building data structures 
  vertex_id_t current_root;
  std::map<vertex_id_t, vertex_id_t> elim_time_map;
  clique_vector cliques;
  std::deque<vertex_id_t> bfs_queue;
  graphlab::mutable_queue<vertex_id_t, double> priority_queue;
  std::set<vertex_id_t> visited;


  // Local junction tree graphlab core
  junction_tree::gl::core jt_core;
 
  


public:

  jt_worker() : 
    scope_factory(NULL), 
    worker_id(0),
    worker_count(0),
    max_tree_size(0),
    total_samples(0),
    collisions(0),
    finish_time_seconds(0),
    use_priorities(false) { }

  void init(size_t wid,
            scope_factory_type& sf, 
            const factorized_model::factor_map_t& factors,
            size_t ncpus,
            float finish_time_secs,
            size_t treesize,
            bool priorities) {
    // Initialize parameters
    scope_factory = &sf;
    worker_id = wid;
    worker_count = ncpus;
    max_tree_size = treesize;
    factors_ptr = &factors;
    use_priorities = priorities;
    finish_time_seconds = finish_time_secs;


    current_root = worker_id;

    // Initialize local jtcore
    // jt_core.set_scheduler_type("fifo");
    // jt_core.set_scope_type("edge");
    // jt_core.set_ncpus(2);
    // jt_core.set_engine_type("async");


    jt_core.set_scheduler_type("fifo");
    jt_core.set_scope_type("none");
    jt_core.set_ncpus(1);
    jt_core.set_engine_type("async_sim");


    // Initialize the shared data in the factorized model
    jt_core.shared_data().set_constant(junction_tree::FACTOR_KEY, 
                                       factors_ptr);
    mrf::graph_type* mrf_graph_ptr = &scope_factory->get_graph();
    jt_core.shared_data().set_constant(junction_tree::MRF_KEY, 
                                       mrf_graph_ptr); 
  }

  size_t num_samples() const { return total_samples; }
  size_t num_collisions() const { return collisions; }

  // get a root
  void run() {
    
    // Track the number of samples
    total_samples = 0;
    collisions = 0;
    // End of for loop
    //size_t round = 0;
    while(graphlab::lowres_time_seconds() < finish_time_seconds) {
      /////////////////////////////////////////////////////////
      // Construct one tree (we must succeed in order to count a tree
      size_t sampled_variables = 0;
      while(sampled_variables == 0) {
        current_root = 
          graphlab::random::rand_int(scope_factory->num_vertices() - 1);
        sampled_variables = sample_once();
        if(sampled_variables == 0) collisions++;
      }


      // // Get a local copy of the graph
      // mrf::graph_type& mrf(scope_factory->get_graph());
      // if(worker_id == 0) {
      //   std::cout << "Saving Image: " << std::endl;
      //   size_t rows = std::sqrt(mrf.num_vertices());
      //   image img(rows, rows);
      //   for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {
      //     vertex_id_t tree_id = mrf.vertex_data(vid).tree_id;
      //     img.pixel(vid) = 
      //       tree_id == vertex_id_t(-1)? 0 : tree_id + worker_count;
      //   }
      //   img.save(make_filename("tree", ".pgm", round++).c_str());
      // }

      // std::cout << "Worker " << worker_id 
      //           << " sampled " << current_root
      //           << " a tree of size " << sampled_variables
      //           << std::endl;        
      total_samples += sampled_variables;
    } 
  }

  /**
   * Grab this vertex into the tree owned by worker id
   */
  bool try_grab_vertex(iscope_type& scope) {
    // Check that this vertex is not already in a tree
    bool in_tree = scope.vertex_data().tree_id != vertex_id_t(-1);
    if(in_tree) return false;

    // check that the neighbors are not in any other trees than this
    // one
    foreach(edge_id_t in_eid, scope.in_edge_ids()) {
      vertex_id_t neighbor_vid = scope.source(in_eid);
      const mrf::vertex_data& vdata = 
        scope.const_neighbor_vertex_data(neighbor_vid);
      bool in_tree = vdata.tree_id != vertex_id_t(-1);
      // if the neighbor is in a tree other than this one quit
      if(in_tree && worker_id != vdata.tree_id) return false;
    }
    // Assert that this vertex is not in a tree and that none of the
    // neighbors are in other trees

    // This vertex does not neighbor any other trees than this one
    scope.vertex_data().tree_id = worker_id;
    return true;
  
  }


  /**
   * Release the vertex
   */
  void release_vertex(iscope_type& scope) {
    // This vertex does not neighbor any other trees than this one
    scope.vertex_data().tree_id = vertex_id_t(-1);  
  }



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

    // start with the variable
    domain_t in_tree_vars = vdata.variable;
    assignment_t out_of_tree_asg;
    foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
      const vertex_id_t neighbor_vid = mrf.source(ineid);
      const mrf::vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
      // test to see if the neighbor is in the tre
      if(elim_time_map.find(neighbor_vid) != elim_time_map.end()) {
        // If this vertex has too many neighbor than the priority is
        // set to 0;
        if(in_tree_vars.num_vars() > MAX_DIM) {
          std::cout << "clique too big: " 
                    << in_tree_vars.num_vars() << std::endl;
          return 0;
        }
        // otherwise add the tree variable
        in_tree_vars += neighbor.variable;
      } else { //outside the tree
        // record the assignment (we will use this condition)
        out_of_tree_asg &= neighbor.asg; 
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

    // std::cout << "Residual: " << residual << std::endl;
    // if(residual > 1000) {
    //   std::cout << clique_factor << std::endl
    //             << conditional_factor  << std::endl
    //             << marginal_factor << std::endl;
    // }

    assert( residual >= 0);

    return residual;
  }

  

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
        // test the 
        bool safe_extension = 
          extend_clique_list(mrf,
                             next_vertex,
                             elim_time_map,
                             cliques);

        if(safe_extension) {   
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
        // test the 
        bool safe_extension = 
          extend_clique_list(mrf,
                             next_vertex,
                             elim_time_map,
                             cliques);
        // If the extension was safe than the elim_time_map and
        // cliques data structure are automatically extended

        if(safe_extension) {                 
          // add the neighbors to the search queue or update their priority
          foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
            vertex_id_t neighbor_vid = mrf.target(eid);
            if(visited.count(neighbor_vid) == 0) {
              double score = score_vertex(neighbor_vid);
              if(score > 0) priority_queue.push(neighbor_vid, score);
              visited.insert(neighbor_vid);
            } else if(priority_queue.contains(neighbor_vid)) {
              // vertex is still in queue we may need to recompute
              // score
              double score = score_vertex(neighbor_vid);
              if(score > 0) {
                // update the priority queue with the new score
                priority_queue.update(neighbor_vid, score);
              } else {
                // The score computation revealed that the clique
                // would be too large so simply remove the vertex from
                // the priority queue
                priority_queue.remove(neighbor_vid);
              }
            } // otherwise the vertex has been visited and processed
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


    std::cout << "Varcount: " << cliques.size() << std::endl;  
    
    // If we failed to build a tree return failure
    if(cliques.empty()) return 0;
    

    // Build the junction tree and sample
    jt_core.graph().clear();
    jtree_from_cliques(mrf, 
                       elim_time_map,
                       cliques.begin(), cliques.end(), 
                       jt_core.graph());
    // Rebuild the engine (clear the old scheduler)
    jt_core.rebuild_engine();
    // add tasks to all vertices
    jt_core.add_task_to_all(junction_tree::calibrate_update, 1.0);
    // Run the core
    std::cout << "Starting engine: " << worker_id << std::endl;
    jt_core.start();
    std::cout << "Finished engine: " << worker_id << std::endl;
    //    std::cout << "Last Update count: " << jt_core.last_update_count() << std::endl;

    // Check that the junction tree is sampled
    for(vertex_id_t vid = 0; 
        vid < jt_core.graph().num_vertices(); ++vid) {
      assert(jt_core.graph().vertex_data(vid).sampled);
    }
      
    // Sampled root successfully
    return cliques.size();
  }



};  // End of JT work









void parallel_sample(const factorized_model& fmodel,
                     mrf::graph_type& mrf,
                     size_t ncpus,
                     float runtime_secs,
                     size_t max_tree_size = 1000,
                     bool use_priorities = false) {
  // create workers
  graphlab::thread_group threads;
  std::vector<jt_worker> workers(ncpus);

  // Create a scope factor
  graphlab::general_scope_factory<mrf::graph_type>
    scope_factory(mrf, ncpus,
                  graphlab::scope_range::EDGE_CONSISTENCY);
  
  float finish_time_secs = 
    graphlab::lowres_time_seconds() + runtime_secs;

  for(size_t i = 0; i < ncpus; ++i) {
    // Initialize the worker
    workers[i].init(i, 
                    scope_factory, 
                    fmodel.factors(),
                    ncpus,
                    finish_time_secs,
                    max_tree_size,
                    use_priorities);    
    // Launch the threads
    bool use_cpu_affinity = false;
    if(use_cpu_affinity) threads.launch(&(workers[i]), i);
    else threads.launch(&(workers[i]));            
  }
  
  // Wait for all threads to finish
  threads.join();

  // Record the total number of samples
  size_t total_samples = 0;
  size_t total_collisions = 0;
  foreach(const jt_worker& worker, workers) {
    total_samples += worker.num_samples();
    total_collisions += worker.num_collisions();
  }
  std::cout << "Total samples: " << total_samples << std::endl;
  std::cout << "Total collisions: " << total_collisions << std::endl;

}










#include <graphlab/macros_undef.hpp>
#endif
