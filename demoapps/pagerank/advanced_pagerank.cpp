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
double RESET_PROB = 0.15;

//! Global accuracy tolerance
double ACCURACY = 1e-5;

/**
 * The factorized page rank update function
 */
class pagerank_update : 
  public graphlab::iupdate_functor<graph_type, pagerank_update> {
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
    ++vdata.nupdates; vdata.value += accum;
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




int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_file;
  std::string format = "metis";
  std::string binfname;
  std::string update_type = "basic";
  size_t topk = 5;
  clopts.attach_option("graph", &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", &format, format,
                       "The graph file format: {metis, snap, tsv}");
  clopts.attach_option("accuracy", &ACCURACY, ACCURACY,
                       "residual termination threshold");
  clopts.attach_option("resetprob", &RESET_PROB, RESET_PROB,
                       "Random reset probability");
  clopts.attach_option("type", &update_type, update_type,
                       "The graphlab update type {basic, delta, factorized}");
  clopts.attach_option("topk", &topk, topk,
                       "The number of top pages to display at the end.");
  clopts.attach_option("binfname", &binfname, binfname,
                       "Optionally save a binary version of the graph");
 
  std::string vinfo_fname;
  clopts.attach_option("vinfo", &vinfo_fname, vinfo_fname,
                       "File containing the vertex info.");
 

  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  UPDATE_STYLE = str2update_style(update_type); 
  std::cout << "Termination bound:  " << ACCURACY 
            << std::endl
            << "Reset probability:  " << RESET_PROB
            << std::endl
            << "Update style:       " << update_style2str(UPDATE_STYLE)
            << std::endl;
  
  // Setup the GraphLab execution core and load graph -------------------------
  graphlab::core<graph_type, pagerank_update> core;
  core.set_options(clopts); // attach the command line options to the core
  if(format == "bin") {
    std::cout << "Loading binary graph." << std::endl;
    core.graph().load(graph_file);
  } else {
    std::cout << "Loading graph from structure file." << std::endl;
    const bool success = graphlab::graph_ops<graph_type>::    
      load_structure(graph_file, format, core.graph());
    if(!success) {
      std::cout << "Error in reading file: " << graph_file << std::endl;
    }
    normalize_graph(core.graph());
  }

  if(!vinfo_fname.empty()) {
    std::cout << "Coloring the graph." << std::endl;
    const size_t num_colors = graphlab::graph_ops<graph_type>::    
      color(core.graph());
    std::ofstream fout(vinfo_fname.c_str());
    for(size_t i = 0; i < core.graph().num_vertices(); ++i) 
      fout 
        << i << '\t' << core.graph().color(i) << '\t'
        << core.graph().num_in_edges(i) << '\t'
        << core.graph().num_out_edges(i) << '\t'
        << graphlab::graph_ops<graph_type>::
        num_neighbors(core.graph(), i)
        << '\n';
    fout.close();
  }

  if(!binfname.empty()) { 
    std::cout << "Saving initial binary version of the graph." << std::endl;
    core.graph().save(binfname);  
    return EXIT_SUCCESS;
  }

  // Run the PageRank ---------------------------------------------------------

  const double initial_delta = RESET_PROB;
  if(UPDATE_STYLE == DELTA) {
    std::cout << "changing initial data" << std::endl;
    for(size_t vid = 0; vid < core.graph().num_vertices(); ++vid) 
      core.graph().vertex_data(vid).value = 0;   
  }

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
  double total_rank = 0;
  for(size_t i = 0; i < core.graph().num_vertices(); ++i) 
    total_rank += core.graph().vertex_data(i).value;
  std::cout << "Total Rank: " << total_rank << std::endl;

  // Output Results -----------------------------------------------------------
  // Output the top 5 pages
  std::vector<graph_type::vertex_id_type> top_pages;
  get_top_pages(core.graph(), topk, top_pages);
  std::cout << std::endl
            << std::setw(10) << "Page Id\t" 
            << std::setw(10) << "Value\t"
            << std::setw(10) << "nupdates\t"
            << std::setw(10) << "indegree\t"
            << std::setw(10) << "outdegree"
            << std::endl;

  for(size_t i = 0; i < top_pages.size(); ++i) {
    const vertex_data& vdata = core.graph().vertex_data(top_pages[i]);
    std::cout << std::setw(10) << top_pages[i] << ":\t" 
              << std::setw(10) << vdata.value << '\t'
              << std::setw(10) << vdata.nupdates << '\t'
              << std::setw(10) 
              << core.graph().in_edges(top_pages[i]).size() << '\t'
              << std::setw(10) 
              << core.graph().out_edges(top_pages[i]).size()
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
