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

#include "dist_pagerank.hpp"



#include <graphlab/macros_def.hpp>



//! Global variable determining update style
update_style UPDATE_STYLE = BASIC;

//! Global random reset probability
double RESET_PROB = 0.15;

//! Global accuracy tolerance
double ACCURACY = 1e-5;

/**
 * The factorized page rank update function
 */
class pagerank_update : 
  public graphlab::iupdate_functor<dist_graph_type, pagerank_update> {
private:
  double accum;
public:
  pagerank_update(const double& accum = 0) : accum(accum) { }
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
  edge_set gather_edges() const { return IN_EDGES; }
  edge_set scatter_edges() const {
    return accum > ACCURACY ? OUT_EDGES : NO_EDGES;
  }

  void delta_functor_update(icontext_type& context) { 
    vertex_data& vdata = context.vertex_data(); 
    if(vdata.nupdates == 0) vdata.value = 0;
    vdata.value += accum; ++vdata.nupdates; 
    if(std::fabs(accum) > ACCURACY || vdata.nupdates == 1) {
      vdata.old_value = vdata.value;
      reschedule_neighbors(context);
    }
  } // end of delta_functor_update


  void operator()(icontext_type& context) {
    // if it is a delta function then use the delta function update
    if(UPDATE_STYLE == DELTA) { delta_functor_update(context); return; }      
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    // Compute weighted sum of neighbors
    double sum = 0;
    foreach(const edge_type& edge, context.in_edges()) {
      sum += context.const_edge_data(edge).weight * 
       context.const_vertex_data(edge.source()).value;
    }
    // Add random reset probability
    vdata.value = RESET_PROB + (1 - RESET_PROB) * sum;
    accum = (vdata.value - vdata.old_value);
    if(std::fabs(accum) > ACCURACY || vdata.nupdates == 1) {
      vdata.old_value = vdata.value;
      reschedule_neighbors(context);
    }
  } // end of operator()  

  // Reset the accumulator before running the gather
  void init_gather(iglobal_context_type& context) { accum = 0; }

  // Run the gather operation over all in edges
  void gather(icontext_type& context, const edge_type& edge) {
    accum +=
      context.const_vertex_data(edge.source()).value *
      context.const_edge_data(edge).weight;
  } // end of gather

  // Merge two pagerank_update accumulators after running gather
  void merge(const pagerank_update& other) { accum += other.accum; }

  // Update the center vertex
  void apply(icontext_type& context) {
    vertex_data& vdata = context.vertex_data(); ++vdata.nupdates;
    vdata.value =  RESET_PROB + (1 - RESET_PROB) * accum;
    accum = vdata.value - vdata.old_value;
    if(std::fabs(accum) > ACCURACY || vdata.nupdates == 1) {
      vdata.old_value = vdata.value;    
      reschedule_neighbors(context);
    }
  } // end of apply

  // Reschedule neighbors 
  void scatter(icontext_type& context, const edge_type& edge) {
    const edge_data& edata = context.const_edge_data(edge);
    const double delta = accum * edata.weight * (1 - RESET_PROB);
    context.schedule(edge.target(), pagerank_update(delta));
  } // end of scatter

private:
  void reschedule_neighbors(icontext_type& context) {
    foreach(const edge_type& edge, context.out_edges()) scatter(context, edge);
  } // end of reschedule neighbors  
}; // end of pagerank update functor
SERIALIZABLE_POD(pagerank_update);





/**
 * This accumulator keeps track of the tpo ranked pages
 */       
class accumulator :
  public graphlab::iaccumulator<dist_graph_type, pagerank_update, 
                                accumulator> {
private:
  typedef std::pair<double, vertex_id_type> pair_type;
  std::priority_queue<pair_type> topk;
  size_t kvalue;
  double total_rank;
public:
  accumulator(size_t kvalue = 5) : kvalue(kvalue), total_rank(0) { }
  void operator()(icontext_type& context) {
    topk.push(pair_type(-context.const_vertex_data().value, 
                       context.vertex_id()));
    if(topk.size() > kvalue) topk.pop();
    total_rank += context.const_vertex_data().value;
  }
  void operator+=(const accumulator& other) { 
    for(std::priority_queue<pair_type> other_topk = other.topk; 
        !other_topk.empty(); other_topk.pop()) {
      topk.push(other_topk.top());
      if(topk.size() > kvalue) topk.pop();
    }
    total_rank += other.total_rank;
  }
  void finalize(iglobal_context_type& context) {
    std::cout << "Pageranks ----------------------------------" << std::endl;
    for( ; !topk.empty(); topk.pop()) 
      std::cout << std::setw(10) << topk.top().second << ":\t"
                << std::setw(10) << -topk.top().first << std::endl;
    std::cout << "Total rank: " << total_rank << std::endl;
  }
  void save(graphlab::oarchive &oarc) const {
    const size_t actual_size = topk.size();
    oarc << actual_size;
    for(std::priority_queue<pair_type> t = topk; !t.empty(); t.pop()) 
      oarc << t.top();
    oarc << kvalue << total_rank; 
  }
  void load(graphlab::iarchive &iarc) { 
    size_t actual_size = 0;
    iarc >> actual_size;
    for(size_t i = 0; i < actual_size; ++i) {
      pair_type pair; iarc >> pair; topk.push(pair);
    }
    iarc >> kvalue >> total_rank; 
  }
}; // end of  accumulator



int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  clopts.use_distributed_options(); 
  std::string graph_file;
  std::string update_type = "basic";
  size_t topk = 5;
  clopts.attach_option("graph", &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("accuracy", &ACCURACY, ACCURACY,
                       "residual termination threshold");
  clopts.attach_option("resetprob", &RESET_PROB, RESET_PROB,
                       "Random reset probability");
  clopts.attach_option("type", &update_type, update_type,
                       "The graphlab update type {basic, delta, factorized}");
  clopts.attach_option("topk", &topk, topk,
                       "The number of top pages to display at the end.");
  if(!clopts.parse(argc, argv) || !clopts.is_set("graph")) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  clopts.print();

  // Compute the update style -------------------------------------------------
  UPDATE_STYLE = str2update_style(update_type); 
  std::cout << "Termination bound:  " << ACCURACY 
            << std::endl
            << "Reset probability:  " << RESET_PROB
            << std::endl
            << "Update style:       " << update_style2str(UPDATE_STYLE)
            << std::endl;

  // Initialize MPI -----------------------------------------------------------
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param param;
  ASSERT_TRUE(graphlab::init_param_from_mpi(param));
  graphlab::distributed_control dc(param);

  // Load the distributed graph -----------------------------------------------
  graphlab::distributed_core<dist_graph_type, pagerank_update> 
    core(dc, graph_file);
  core.set_options(clopts);
  core.build_engine();

  // Initialize the topk accumulator ------------------------------------------
  const size_t sync_interval = core.graph().num_vertices();  
  core.add_sync("topk", accumulator(topk) , sync_interval);
  
  // Run the PageRank ---------------------------------------------------------
  const double initial_delta = RESET_PROB;
  core.schedule_all(pagerank_update(initial_delta));

  std::cout << "Running pagerank!" << std::endl;
  const double runtime = core.start();  // Run the engine
  std::cout << "Graphlab finished, runtime: " << runtime << " seconds." 
            << std::endl
            << "Updates executed: " << core.last_update_count() 
            << std::endl
            << "Update Rate (updates/second): " 
            << core.last_update_count() / runtime
            << std::endl;

  // Display final output -----------------------------------------------------
  std::cout << "Recomputing topk" << std::endl;
  core.sync_now("topk");

  std::cout << "FINISHED!" << std::endl;
  dc.barrier();
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main










