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
 
/**
 * @file dispatcher.cpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */

#include <boost/algorithm/string.hpp>
 
#include <signal.h>
#include "dispatcher.hpp"
#include "process.hpp"
#include "json_message.hpp"

#include <graphlab/macros_def.hpp>

using namespace graphlab;
typedef json_invocation ji;
typedef json_return     jr;

///////////////////////////////// CGI_GATHER ///////////////////////////////////
cgi_gather::cgi_gather(const std::string& state) : mstate(state){}

const char *cgi_gather::c_str() const {
  return mstate.c_str();
}

void cgi_gather::save(oarchive& oarc) const {
  oarc << mstate;
}

void cgi_gather::load(iarchive& iarc) {
  iarc >> mstate;
}

void cgi_gather::operator+=(const cgi_gather& other){
  
  // invoke
  process& p = process::get_process();
  ji invocation("merge", mstate);
  invocation.add_other(other.mstate);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  // update states (if none received, no change)
  const char *cstring = result.result();
  if (NULL != cstring)
    mstate = std::string(cstring);
  
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DISPATCHER //////////////////////////////////
dispatcher::dispatcher(const std::string& state) : mstate(state){}

void dispatcher::save(oarchive& oarc) const {
  oarc << mstate;
}

void dispatcher::load(iarchive& iarc){
  iarc >> mstate;
}

edge_dir_type dispatcher::gather_edges
                             (icontext_type& context,
                              const vertex_type& vertex) const {
  
  // invoke
  process& p = process::get_process();
  ji invocation("gather_edges", mstate);
  invocation.add_vertex(vertex);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  return result.edge_dir();
  
}

edge_dir_type dispatcher::scatter_edges
                            (icontext_type& context,
                             const vertex_type& vertex) const {
  
  // invoke
  process& p = process::get_process();
  ji invocation("scatter_edges", mstate);
  invocation.add_vertex(vertex);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  return result.edge_dir();
  
}

void dispatcher::signal(icontext_type& context,
                        edge_type& edge,
                        const jr& result) const {
  
  // null means no signal                      
  const char *cstring = result.signal();
  if (NULL == cstring) return;

  // either signal source or target
  if (!strcmp("SOURCE", cstring)){
    context.signal(edge.source());
  }else if (!strcmp("TARGET", cstring)){
    context.signal(edge.target());
  }else logstream(LOG_FATAL) << "Unrecognized signal: " << cstring << std::endl;
  
}

dispatcher::gather_type dispatcher::gather
                             (icontext_type& context,
                              const vertex_type& vertex,
                              edge_type& edge) const {

  // invoke
  process& p = process::get_process();
  ji invocation("gather", mstate);
  invocation.add_vertex(vertex);
  invocation.add_edge(edge);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  // signal
  signal(context, edge, result);
  
  // return result
  const char *cstring = result.result();
  return cgi_gather(std::string(cstring));

}

void dispatcher::apply(icontext_type& context,
                       vertex_type& vertex,
                       const gather_type& total) {

  // invoke
  process &p = process::get_process();
  ji invocation("apply", mstate);
  invocation.add_vertex(vertex);
  invocation.add_gather(total);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  // update states (null means no change)
  const char *cstring = result.program();
  if (NULL != cstring)
    mstate = std::string(cstring);                      // copy
  cstring = result.vertex();
  if (NULL != cstring)
    vertex.data() = std::string(cstring);               // copy

}

void dispatcher::scatter(icontext_type& context,
                         const vertex_type& vertex,
                         edge_type& edge) const {
                         
  process& p = process::get_process();
  
  // invoke
  ji invocation("scatter", mstate);
  invocation.add_vertex(vertex);
  invocation.add_edge(edge);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);

  // signal
  signal(context, edge, result);
  
}
 
////////////////////////////////////////////////////////////////////////////////

struct dispatcher_writer {
  std::string save_vertex(dispatcher::vertex_type vertex) {
    std::stringstream strm;
    strm << ":" << vertex.id() << "\t" << vertex.data() << "\n";
    return strm.str();
  }
  std::string save_edge(dispatcher::edge_type edge) {
    std::stringstream strm;
    strm << edge.source().id() << "\t" << edge.target().id() << "\t" << edge.data() << "\n";
    return strm.str();
  }
};

bool dispatcher::read_graph(graph_type& graph,
                           const std::string& filename,
                           const std::string& textline) {
  
  // skip blank lines
  if (textline.empty()) return true;
  
  switch(textline[0]){
    case '#': return true;    // skip comments
    case ':': {               // load vertices
      if (textline.length() <= 1) return false;
      std::istringstream iss(textline.substr(1));
      dispatcher::graph_type::vertex_id_type vid;
      dispatcher::graph_type::vertex_data_type vdata;
      iss >> vid >> vdata;
      graph.add_vertex(vid, vdata);
      break;
    }
    default: {                // load edges
      if (textline.length() <= 2) return false;
      std::istringstream iss(textline);
      dispatcher::graph_type::vertex_id_type source;
      dispatcher::graph_type::vertex_id_type target;
      dispatcher::graph_type::edge_data_type edata;
      iss >> source >> target >> edata;
      if (source == target){
        logstream(LOG_ERROR) << "Self edge not allowed: " << source << std::endl;
      }else {
        graph.add_edge(source, target, edata);
      }
      break;
    }
  }
  
  return true;
  
}

void dispatcher::init_edge(graph_type::edge_type edge) {

  process& p = process::get_process();
  
  // invoke (not passing state, because we are not supposed to be in a vertex program)
  std::string method("init_edge");
  ji invocation(method);
  invocation.add_edge(edge);
  p.send(invocation);
  
  // parse results
  jr result;
  p.receive(result);

  // null means no change
  const char *cstring = result.edge();
  if (NULL == cstring) return;
  
  edge.data() = std::string(cstring);
  
}
 
void dispatcher::init_vertex(graph_type::vertex_type vertex) {
  
  process& p = process::get_process();
  
  // invoke (not passing state, because we are not supposed to be in a vertex program)
  std::string method("init_vertex");
  ji invocation(method);
  invocation.add_vertex(vertex);
  p.send(invocation);
  
  // parse results
  jr result;
  p.receive(result);

  // null means no change
  const char *cstring = result.vertex();
  if (NULL == cstring) return;
  
  vertex.data() = std::string(cstring);
  
}

/////////////////////////// HELPERS FOR MAIN METHOD ////////////////////////////

static void ignore_sigpipe(){
  struct sigaction sa;
  std::memset(&sa, 0, sizeof(sa));
  sa.sa_handler = SIG_IGN;
  CHECK(!::sigaction(SIGPIPE, &sa, NULL));
}

////////////////////////////////// MAIN METHOD /////////////////////////////////
int main(int argc, char** argv) {

  global_logger().set_log_level(LOG_DEBUG);

  // Initialize control plain using mpi
  mpi_tools::init(argc, argv);
  distributed_control dc;
  
  // Parse command line options ------------------------------------------------
  command_line_options clopts("GraphLab Dispatcher");
  
  std::string graph_dir = "toy.tsv";
  clopts.attach_option("graph", &graph_dir, graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  std::string format = "tsv";
  clopts.attach_option("format", &format, format,
                       "The graph file format: {metis, snap, tsv, adj, bin, cgi}");
  std::string program;
  clopts.attach_option("program", &program,
                       "The program to execute (required)");
  std::string saveprefix;
  clopts.attach_option("saveprefix", &saveprefix,
                       "If set, will save the resultant graph to a "
                       "sequence of files with prefix saveprefix");
  std::string signal;
  clopts.attach_option("signal", &signal,
                       "If set, then only vertices with these IDs will be signalled "
                       "when the engine is started.");
                       
  if(!clopts.parse(argc, argv)) {
    logstream(LOG_ERROR) << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  
  if(!clopts.is_set("program")) {
    logstream(LOG_ERROR) << "Program not provided." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  
  // Signal handling -----------------------------------------------------------
  ignore_sigpipe();
  
  // Load graph ----------------------------------------------------------------
  logstream(LOG_INFO) << dc.procid() << ": Starting." << std::endl;
    
  logstream(LOG_EMPH) << "Loading graph from " << graph_dir << std::endl;
  dispatcher::graph_type graph(dc, clopts);
  if (0 == format.compare("cgi")) graph.load(graph_dir, dispatcher::read_graph);
  else graph.load_format(graph_dir, format);
  graph.finalize();
  
  logstream(LOG_EMPH) << "#vertices: " << graph.num_vertices()
                      << " #edges:" << graph.num_edges() << std::endl;
  logstream(LOG_INFO) << "executable: " << program;

  process::set_executable(program);
  
  std::stringstream arg1; arg1 << "--num-vertices" << graph.num_vertices();
  std::stringstream arg2; arg2 << "--num-edges" << graph.num_edges();
  process::add_arg(arg1.str());
  process::add_arg(arg2.str());
  
  // Create engine -------------------------------------------------------------
  logstream(LOG_INFO) << dc.procid() << ": Creating engine" << std::endl;
  omni_engine<dispatcher> engine(dc, graph, "synchronous", clopts);
  
  graph.transform_edges(dispatcher::init_edge);
  graph.transform_vertices(dispatcher::init_vertex);
  
  // Schedule vertices ---------------------------------------------------------
  if (clopts.is_set("signal")){
    std::vector<std::string> ids;
    boost::split(ids, signal, boost::is_any_of(", "), boost::token_compress_on);
    foreach(std::string id, ids){
      dispatcher::graph_type::vertex_id_type vid;
      std::istringstream(id) >> vid;
      engine.signal(vid);
    }
  } else engine.signal_all();
  
  
  // Run the dispatcher --------------------------------------------------------
  engine.start();
  
  const float runtime = engine.elapsed_seconds();
  size_t update_count = engine.num_updates();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;

  // Save output ---------------------------------------------------------------
  if(clopts.is_set("saveprefix"))
    graph.save(saveprefix, dispatcher_writer(), false);
  
  // Teardown communication tools ----------------------------------------------
  mpi_tools::finalize();
  ::kill(-getpid(), SIGTERM);
  return EXIT_SUCCESS;
  
}
////////////////////////////////////////////////////////////////////////////////

#include <graphlab/macros_undef.hpp>
