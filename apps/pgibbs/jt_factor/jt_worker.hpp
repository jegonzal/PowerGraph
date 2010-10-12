#ifndef JT_WORKER_HPP
#define JT_WORKER_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cassert>



// Including Standard Libraries

#include <graphlab.hpp>

#include <graphlab/parallel/pthread_tools.hpp>


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
  bool active;

  // Tree building data structures 
  vertex_id_t current_root;
  std::map<vertex_id_t, vertex_id_t> elim_time_map;
  clique_vector cliques;
  std::queue<vertex_id_t> bfs_queue;
  std::set<vertex_id_t> visited;

  // Local junction tree graphlab core
  junction_tree::gl::core jt_core;
 
  


public:

  jt_worker() : scope_factory(NULL) { }

  void init(size_t wid,
            size_t wcount,
            scope_factory_type& sf, 
            const factorized_model::factor_map_t& factors) {
    // Initialize parameters
    scope_factory = &sf;
    worker_id = wid;
    worker_count = wcount;
    active = true;

    current_root = worker_id;

    // Initialize local jtcore
    jt_core.set_scheduler_type("fifo");
    jt_core.set_scope_type("none");
    jt_core.set_ncpus(1);
    jt_core.set_engine_type("async_sim");

    // Initialize the shared data in the factorized model
    const  factorized_model::factor_map_t* factors_ptr = &factors;
    jt_core.shared_data().set_constant(junction_tree::FACTOR_KEY, 
                                       factors_ptr);
    mrf::graph_type* mrf_graph_ptr = &scope_factory->get_graph();
    jt_core.shared_data().set_constant(junction_tree::MRF_KEY, 
                                       mrf_graph_ptr); 
  }

  // get a root
  void run() {
    size_t total_samples = 0;
    for(size_t i = 0; i < 10; ++i) {
      size_t sampled_variables = sample_once();

      // Update root
      // current_root += worker_count;
      // if(current_root >= scope_factory->num_vertices()) {
      //   current_root = worker_id;
      // }
      current_root = 
        graphlab::random::rand_int(scope_factory->num_vertices() - 1);

      std::cout << "Worker " << worker_id 
                << " sampled " << current_root
                << " a tree of size " << sampled_variables
                << std::endl;        
      total_samples += sampled_variables;
    }      
  }

  /**
   * Grab this vertex into the tree owned by worker id
   */
  bool try_grab_vertex(iscope_type& scope) {
    // check that the neighbors are not in any other trees than this one
    foreach(edge_id_t in_eid, scope.in_edge_ids()) {
      vertex_id_t neighbor_vid = scope.source(in_eid);
      const mrf::vertex_data& vdata = 
        scope.const_neighbor_vertex_data(neighbor_vid);
      bool in_tree = vdata.tree_id != vertex_id_t(-1);
      if(in_tree && worker_id != vdata.tree_id) return false;
    }
    // This vertex does not neighbor any other trees than this one
    scope.vertex_data().tree_id = worker_id;
    return true;
  
  }


  /**
   * Release the vertex
   */
  void release_vertex(iscope_type& scope) {
    // This vertex does not neighbor any other trees than this one
    scope.vertex_data().tree_id = -1;  
  }



  

  size_t sample_once() {
    std::cout << "Sample once: " << worker_id << std::endl;
    assert(scope_factory != NULL);

    // Get the scope factory
    mrf::graph_type& mrf = scope_factory->get_graph();

    // Clear local data structures
    elim_time_map.clear();
    cliques.clear();
    assert(bfs_queue.empty());
    visited.clear();
     
    // add the root
    bfs_queue.push(current_root);
    visited.insert(current_root);

    while(!bfs_queue.empty()) {
      // Take the top element
      const vertex_id_t next_vertex = bfs_queue.front();
      bfs_queue.pop(); 

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
              bfs_queue.push(neighbor_vid);
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
    } // end of while loop
  
    std::cout << "Varcount: " << cliques.size() << std::endl;  
    // If we failed to build a tree return failure
    if(cliques.empty()) return 0;

    // Build the junction tree and sample
    jt_core.graph().clear();
    junction_tree_from_cliques(mrf, 
                               cliques.begin(), cliques.end(), 
                               jt_core.graph());
    // add tasks to all vertices
    jt_core.add_task_to_all(junction_tree::calibrate_update, 1.0);
    // Run the core
    jt_core.start();
    std::cout << "Last Update count: " << jt_core.last_update_count() << std::endl;

    // Check that the entire tree is cleared
    foreach(elim_clique clique, cliques) {
      const mrf::vertex_data& vdata = 
        mrf.vertex_data(clique.elim_vertex);
      assert(vdata.tree_id == vertex_id_t(-1));
    } // end of foreach
      
    // Sampled root successfully
    return cliques.size();
  }



};  // End of JT work









void parallel_sample(const factorized_model& fmodel,
                     mrf::graph_type& mrf,
                     size_t ncpus) {
  // create workers
  graphlab::thread_group threads;
  std::vector<jt_worker> workers(ncpus);

  // Create a scope factor
  graphlab::general_scope_factory<mrf::graph_type>
    scope_factory(mrf, ncpus,
                  graphlab::scope_range::EDGE_CONSISTENCY);
  
  for(size_t i = 0; i < ncpus; ++i) {
    // Initialize the worker
    workers[i].init(i, ncpus, scope_factory, fmodel.factors());    
    // Launch the threads
    bool use_cpu_affinity = false;
    if(use_cpu_affinity) threads.launch(&(workers[i]), i);
    else threads.launch(&(workers[i]));            
  }
  
  // Wait for all threads to finish
  threads.join();
}










#include <graphlab/macros_undef.hpp>
#endif
