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


/*
 *  PageRank tutorial application.
 *  pagerank.cpp
 *
 *  For description of the PageRank algorithm, see Wikipedia article
 *  http://en.wikipedia.org/wiki/Pagerank
 */

#include "pagerank.hpp"



#include <graphlab/macros_def.hpp>

double termination_bound = 1e-5;
double random_reset_prob = 0.15;   // PageRank random reset probability

#ifdef FACTORIZED

/**
 * The factorized page rank update function
 */
class pagerank_update : public gl::iupdate_functor::factorized {
private:
  double prio;
  double accum;
public:
  pagerank_update(const double& prio = 0) : prio(prio), accum(0) { }
  double priority() const { return prio; }
  void operator+=(const pagerank_update& other) { 
    prio += other.prio;
    accum += other.accum;
  }

  void gather(gl::iscope& scope, gl::icallback& callback, 
              gl::edge_id in_eid) {
    // Get the neighobr vertex value
    const vertex_data& neighbor_vdata =
      scope.const_neighbor_vertex_data(scope.source(in_eid));
    const double neighbor_value = neighbor_vdata.value;    
    // Get the edge data for the neighbor
    edge_data& edata = scope.edge_data(in_eid);
    // Compute the contribution of the neighbor
    double contribution = edata.weight * neighbor_value;    
    // Add the contribution to the sum
    accum += contribution;
    // Remember this value as last read from the neighbor
    edata.old_source_value = neighbor_value;
  } // end of gather

  void apply(gl::iscope& scope,
             gl::icallback& callback) {                       
    // Get the data associated with the vertex
    vertex_data& vdata = scope.vertex_data();
    // add the contribution from a self-link.
    accum += vdata.value * vdata.self_weight;
    // add the random reset probability
    accum = random_reset_prob/scope.num_vertices() + 
      (1-random_reset_prob)*accum;
    vdata.value = accum;
  } // end of apply

  void scatter(gl::iscope& scope, gl::icallback& callback, 
               gl::edge_id out_eid) {
    // Get the data associated with the vertex
    const vertex_data& vdata = scope.const_vertex_data();
    // get the data associated with the out edge
    const edge_data& outedgedata = scope.const_edge_data(out_eid);    
    // Compute edge-specific residual by comparing the new value of this
    // vertex to the previous value seen by the neighbor vertex.
    double residual =
      outedgedata.weight *
      std::fabs(outedgedata.old_source_value - vdata.value);
    // If the neighbor changed sufficiently add to scheduler.
    if(residual > termination_bound) {
      callback.schedule(scope.target(out_eid), pagerank_update(residual));
    }
  } // end of scatter
  
}; // end of pagerank update functor

#else

/**
 * The Page rank update function
 */
class pagerank_update : public gl::iupdate_functor {
private:
  double prio;  
public:
  pagerank_update(const double& prio = 0) : prio(prio) { }
  double priority() const { return prio; }
  void operator+=(const pagerank_update& other) { prio += other.prio; }

  void operator()(gl::iscope& scope,
                  gl::icallback& callback) {                       
    // Get the data associated with the vertex
    vertex_data& vdata = scope.vertex_data();
  
    // Sum the incoming weights; start by adding the 
    // contribution from a self-link.
    float sum = vdata.value * vdata.self_weight;
    // Loop over all in edges to this vertex
    const gl::edge_list in_edges = scope.in_edge_ids();
    foreach(gl::edge_id eid, in_edges) {
      // Get the neighobr vertex value
      const vertex_data& neighbor_vdata =
        scope.const_neighbor_vertex_data(scope.source(eid));
      const double neighbor_value = neighbor_vdata.value;    
      // Get the edge data for the neighbor
      edge_data& edata = scope.edge_data(eid);
      // Compute the contribution of the neighbor
      double contribution = edata.weight * neighbor_value;    
      // Add the contribution to the sum
      sum += contribution;
      // Remember this value as last read from the neighbor
      edata.old_source_value = neighbor_value;
    }

    // compute the jumpweight
    sum = random_reset_prob/scope.num_vertices() + 
      (1-random_reset_prob)*sum;
    vdata.value = sum;
   
    // Schedule the neighbors as needed
    foreach(gl::edge_id eid, scope.out_edge_ids()) {
      edge_data& outedgedata = scope.edge_data(eid);    
      // Compute edge-specific residual by comparing the new value of this
      // vertex to the previous value seen by the neighbor vertex.
      double residual =
        outedgedata.weight *
        std::fabs(outedgedata.old_source_value - vdata.value);
      // If the neighbor changed sufficiently add to scheduler.
      if(residual > termination_bound) {
        callback.schedule(scope.target(eid), pagerank_update(residual));
      }
    }
  } // end of operator()
}; // end of pagerank update functor

#endif






int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");

  // Metrics  
  graphlab::metrics app_metrics("app::pagerank");

  // Setup the parser
  graphlab::command_line_options
    clopts("Run the PageRank algorithm.");

  // Add some command line options
  std::string graph_file;
  clopts.attach_option("graph",
                       &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.attach_option("bound",
                       &termination_bound, termination_bound,
                       "residual termination threshold");
  clopts.attach_option("resetprob",
                       &random_reset_prob, random_reset_prob,
                       "Random reset probability");

  // Parse the command line input
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing input." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Using termination bound:  " << termination_bound 
            << std::endl
            << "Using reset probability:  " << random_reset_prob
            << std::endl;
  
  // Create a graphlab core
  gl::core core;

  // Set the engine options
  core.set_options(clopts);
  
  // Create or load graph depending on if the file was set
  app_metrics.start_time("load");
  if(graph_file.empty()) {
    // Create a synthetic graph
    make_toy_graph(core.graph());
  } else {
    // load the graph from the file
    bool success = load_graph_from_metis_file(graph_file, core.graph());
    if(!success) {
      std::cout << "Error in reading file: " << graph_file
                << std::endl;
      return EXIT_FAILURE;
    }
  }
  app_metrics.stop_time("load");

  // std::cout << "Saving graph as tsv" << std::endl;
  // save_edges_as_tsv("edges.tsv", core.graph());
  // std::cout << "Finished saving graph." << std::endl;

  //  app_metrics.report(core.get_reporter());
  // Schedule all vertices to run pagerank update on the
  // first round.
  core.schedule_all(pagerank_update(100.0));
  
  // Run the engine
  double runtime = core.start();
  
  // We are done, now output results.
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  
  // First we need to compute a normalizer. This could be done with
  // the sync facility, but for simplicity, we do it by hand.
  double norm = 0.0;
  for(gl::vertex_id vid = 0; 
      vid < core.graph().num_vertices(); vid++) {
    norm += core.graph().vertex_data(vid).value;
  }
  std::cout << "Total Mass: " << norm << std::endl;

  
  // And output 5 first vertices pagerank after dividing their value
  // with the norm.
  for(gl::vertex_id vid = 0; 
      vid < 5 && vid < core.graph().num_vertices(); vid++) {
    std::cout << "Page " << vid << " pagerank = " <<
      core.graph().vertex_data(vid).value << '\n';
  }    

  // Write the pagerank vector
  std::cout << "Saving pagerank vector." << std::endl;
  save_pagerank("pagerank.tsv", core.graph());
  std::cout << "Finished." << std::endl;

  return EXIT_SUCCESS;
} // End of main


