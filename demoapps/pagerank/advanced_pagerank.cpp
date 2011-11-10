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


/**
 * The type of update to use for pagerank
 *  
 *  Basic updates correspond to the classic GraphLab v1 style of
 *  scoped computation
 *
 *  Delta updates use GraphLab v2 functors to send commutative
 *  associative "messages" to vertices
 *
 *  Factorized updates use GraphLab v2 factorized functors to
 *  decompose a PageRank update over the edges of a vertex.
 */
enum update_style {BASIC, DELTA, FACTORIZED};
update_style str2update_style(std::string str);
std::string update_style2str(update_style style);

//! Global variable determining update style
update_style UPDATE_STYLE = BASIC;

//! Global random reset probability
double RANDOM_RESET_PROBABILITY = 0.15;

//! Global accuracy tolerance
double ACCURACY = 1e-5;

/**
 * The factorized page rank update function
 */
class pagerank_update : 
  public graphlab::iupdate_functor<graph_type, pagerank_update> {
  typedef graphlab::iupdate_functor<graph_type, pagerank_update> base;
  typedef base::icontext_type   icontext_type;
  typedef base::edge_list_type  edge_list_type;
  typedef base::edge_id_type    edge_id_type;
  typedef base::vertex_id_type  vertex_id_type;
private:
  float accum;
public:
  pagerank_update(const float& accum = 0) : accum(accum) { }
  double priority() const { return std::fabs(accum); }
  void operator+=(const pagerank_update& other) { accum += other.accum; }
  bool is_factorizable() const { return UPDATE_STYLE == FACTORIZED; }
  graphlab::consistency_model::model_enum consistency() const {
    if(UPDATE_STYLE == DELTA) 
      return graphlab::consistency_model::VERTEX_CONSISTENCY;
    else return graphlab::consistency_model::USE_DEFAULT;
  }
  bool writable_gather() { return false; }
  bool writable_scatter() { return false; }

  void delta_functor_update(icontext_type& context) { 
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    vdata.old_value = vdata.value;
    vdata.value +=
      ((1-RANDOM_RESET_PROBABILITY)/
       (1-(1-RANDOM_RESET_PROBABILITY)*vdata.self_weight)) * accum;
    reschedule_neighbors(context);
  } // end of delta_functor_update


  void operator()(base::icontext_type& context) {
    // if it is a delta function then use the delta function update
    if(UPDATE_STYLE == DELTA) { delta_functor_update(context); return; }      
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    // Compute weighted sum of neighbors
    float sum = 0;
    foreach(base::edge_id_type eid, context.in_edge_ids()) 
      sum += context.edge_data(eid).weight * 
        context.neighbor_vertex_data(context.source(eid)).value;
    const float self_term = 1-(1-RANDOM_RESET_PROBABILITY)*vdata.self_weight; 
    // Add random reset probability
    vdata.old_value = vdata.value;
    vdata.value = 
      RANDOM_RESET_PROBABILITY/(context.num_vertices()*self_term) +
      sum * (1-RANDOM_RESET_PROBABILITY)/self_term;
    reschedule_neighbors(context);
  } // end of operator()  

  // Reset the accumulator before running the gather
  void init_gather() { accum = 0; }

  // Run the gather operation over all in edges
  void gather(icontext_type& context, edge_id_type in_eid) {
    const vertex_data& neighbor_vdata =
      context.const_neighbor_vertex_data(context.source(in_eid));
    const double neighbor_value = neighbor_vdata.value;    
    edge_data& edata = context.edge_data(in_eid);
    accum += edata.weight * neighbor_value;    
  } // end of gather

  // Merge two pagerank_update accumulators after running gather
  void merge(const pagerank_update& other) { accum += other.accum; }

  // Update the center vertex
  void apply(icontext_type& context) {
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    accum += vdata.value * vdata.self_weight;
    accum = RANDOM_RESET_PROBABILITY/context.num_vertices() + 
      (1-RANDOM_RESET_PROBABILITY)*accum;
    vdata.old_value = vdata.value;
    vdata.value = accum;
  } // end of apply

  // Reschedule neighbors
  void scatter(icontext_type& context, edge_id_type out_eid) {
    const vertex_data& vdata = context.const_vertex_data();
    const edge_data& edata   = context.const_edge_data(out_eid);    
    const double residual = 
      edata.weight*(vdata.value - vdata.old_value);
    if(std::fabs(residual) > ACCURACY || vdata.nupdates == 1) {
      context.schedule(context.target(out_eid), pagerank_update(residual));
    }
  } // end of scatter

private:
  void reschedule_neighbors(icontext_type& context) {
    foreach(edge_id_type eid, context.out_edge_ids()) scatter(context, eid);
  } // end of reschedule neighbors
  
}; // end of pagerank update functor




int main(int argc, char** argv) {
  logger(LOG_INFO, "PageRank starting\n");
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_file;
  std::string format = "metis";
  std::string update_type = "basic";
  size_t topk = 5;
  clopts.attach_option("graph",
                       &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.attach_option("format",
                       &format, format,
                       "The graph file format: {metis, jure, tsv}");
  clopts.attach_option("accuracy",
                       &ACCURACY, ACCURACY,
                       "residual termination threshold");
  clopts.attach_option("resetprob",
                       &RANDOM_RESET_PROBABILITY, RANDOM_RESET_PROBABILITY,
                       "Random reset probability");
  clopts.attach_option("type",
                       &update_type, update_type,
                       "The graphlab update type {basic, delta, factorized}");
  clopts.attach_option("topk",
                       &topk, topk,
                       "The number of top pages to display at the end.");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  UPDATE_STYLE = str2update_style(update_type); 
  std::cout << "Termination bound:  " << ACCURACY 
            << std::endl
            << "Reset probability:  " << RANDOM_RESET_PROBABILITY
            << std::endl
            << "Update style:       " << update_style2str(UPDATE_STYLE)
            << std::endl;
  
  // Setup the GraphLab execution core and load graph -------------------------
  graphlab::core<graph_type, pagerank_update> core;
  core.set_options(clopts); // attach the command line options to the core
  const bool success = load_graph(graph_file, format, core.graph());
  if(!success) {
    std::cout << "Error in reading file: " << graph_file << std::endl;
  }

  // Run the PageRank ---------------------------------------------------------
  double initial_delta = 0;
  if(UPDATE_STYLE == DELTA) {
    std::cout << "changing initial delta" << std::endl;
    initial_delta = RANDOM_RESET_PROBABILITY / 
      (core.graph().num_vertices() * (1-RANDOM_RESET_PROBABILITY));
    for(graph_type::vertex_id_type vid = 0; vid < core.graph().num_vertices();
        ++vid) {
      core.graph().vertex_data(vid).value = 0;
    }
  }
  core.schedule_all(pagerank_update(initial_delta));
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
  get_top_pages(core.graph(), topk, top_pages);
  for(size_t i = 0; i < top_pages.size(); ++i) {
    const vertex_data& vdata = core.graph().vertex_data(top_pages[i]);
    std::cout << std::setw(10) << top_pages[i] << ":\t" 
              << std::setw(10) << vdata.value << '\t'
              << std::setw(10) << vdata.nupdates << '\t'
              << std::setw(10) 
              << core.graph().num_in_neighbors(top_pages[i]) << '\t'
              << std::setw(10) 
              << core.graph().num_out_neighbors(top_pages[i])
              << std::endl;
  }


  // Write the pagerank vector
  std::cout << "Saving pagerank vector." << std::endl;
  save_pagerank("pagerank.tsv", core.graph());
  std::cout << "Finished." << std::endl;

  return EXIT_SUCCESS;
} // End of main







update_style str2update_style(std::string str) {
  if(str == "basic") return BASIC;  
  else if(str == "delta") return DELTA;
  else if(str == "factorized") return FACTORIZED;
  else {
    logstream(LOG_WARNING) 
      << "Invalid update style \"" << str 
      <<"\" reverting to basic!" << std::endl;
    return BASIC;
  }  
} // end of str2update_style;



std::string update_style2str(update_style style) {
  switch(style) {
  case BASIC: { return "basic"; }
  case DELTA: { return "delta"; }
  case FACTORIZED: { return "factorized"; }
  }
  return "";
} // end of str2update_style;
