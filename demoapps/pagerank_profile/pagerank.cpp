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

#include <iomanip>

#include "pagerank.hpp"



#include <graphlab/macros_def.hpp>

//! Global random reset probability
double RANDOM_RESET_PROBABILITY = 0.15;

//! Global accuracy tolerance
double ACCURACY = 1e-5;

/**
 * The factorized page rank update function
 */
class pagerank_update : 
  public graphlab::iupdate_functor<graph_type, pagerank_update> {
  double prio;
public:
  pagerank_update(const double& prio = 0) : prio(prio) { }
  double priority() const { return prio; }
  void operator+=(const pagerank_update& other) { prio += other.prio; }
  void operator()(icontext_type& context) {
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    // Compute weighted sum of neighbors
    double sum = 0;
    /* Iterate over edge_id_list and get source is slow in graph2 */
    foreach( edge_type edge, context.in_edges() ) 
      sum += context.const_edge_data(edge).weight * 
        context.const_vertex_data(edge.source()).value;
    // Add random reset probability
    vdata.old_value = vdata.value;
    vdata.value = RANDOM_RESET_PROBABILITY/context.num_vertices() + 
      (1-RANDOM_RESET_PROBABILITY)*sum; 
    foreach(edge_type edge, context.out_edges()) {    
      const double residual = context.const_edge_data(edge).weight * 
        std::fabs(vdata.value - vdata.old_value);
      // If the neighbor changed sufficiently add to scheduler.
      if(residual > ACCURACY) 
        context.schedule(edge.target(), pagerank_update(residual));      
    }
  } // end of operator()  
}; // end of pagerank update functor




int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_file;
  std::string format = "metis";
  std::string update_type = "basic";
  clopts.attach_option("graph",
                       &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format",
                       &format, format,
                       "The graph file format: {metis, snap, tsv}");
  clopts.attach_option("accuracy",
                       &ACCURACY, ACCURACY,
                       "residual termination threshold");
  clopts.attach_option("resetprob",
                       &RANDOM_RESET_PROBABILITY, RANDOM_RESET_PROBABILITY,
                       "Random reset probability");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Termination bound:  " << ACCURACY 
            << std::endl
            << "Reset probability:  " << RANDOM_RESET_PROBABILITY
            << std::endl;  
  // Setup the GraphLab execution core and load graph -------------------------
  graphlab::core<graph_type, pagerank_update> core;
  core.set_options(clopts); // attach the command line options to the core
  std::cout << "Loading graph from file" << std::endl;
  const bool success = graphlab::graph_ops<graph_type>::
    load_structure(graph_file, format, core.graph());
  if(!success) {
    std::cout << "Error in reading file: " << graph_file << std::endl;
  }
  normalize_graph(core.graph());

  // Run the PageRank ---------------------------------------------------------
  core.schedule_all(pagerank_update(1));
  std::cout << "Running pagerank!" << std::endl;
  const double runtime = core.start();  // Run the engine
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  std::cout << "Update Rate (updates/second): " 
            << core.last_update_count() / runtime
            << std::endl;

 

  // Output Results -----------------------------------------------------------
  // Output the top 5 pages
  std::vector<graph_type::vertex_id_type> top_pages;
  get_top_pages(core.graph(), 5, top_pages);
  for(size_t i = 0; i < top_pages.size(); ++i) {
    std::cout << std::setw(10) << top_pages[i] << ":\t" << std::setw(10) 
              << core.graph().vertex_data(top_pages[i]).value << std::setw(10) 
              << core.graph().vertex_data(top_pages[i]).nupdates << std::endl;
  }

  // Write the pagerank vector
  std::cout << "Saving pagerank vector." << std::endl;
  save_pagerank("pagerank.tsv", core.graph());
  std::cout << "Finished." << std::endl;

  return EXIT_SUCCESS;
} // End of main
