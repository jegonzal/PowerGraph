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
 
#include "dispatcher.hpp"
#include "process.hpp"
#include "json_message.hpp"
#include "dispatcher_ops.hpp"

#include <graphlab/macros_def.hpp>

using namespace graphlab;
typedef json_schedule   js;
typedef json_invocation ji;
typedef json_return     jr;

/////////////////////////////// CGI_GATHER_TYPE ////////////////////////////////
void cgi_gather_type::save(oarchive& oarc) const {
  oarc << mstate;
}

void cgi_gather_type::load(iarchive& iarc) {
  iarc >> mstate;
}

void operator+=(const cgi_gather_type& other){
  
  // invoke
  process& p = process::get_process();
  ji invocation("merge", mstate);
  invocation.add_other(other.mstate);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  // update states
  const char *cstring = result.updater();
  mstate = std::string(cstring);
  
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DISPATCHER //////////////////////////////////

dp::dispatcher_update(const std::string& state) : mstate(state) {}

dp::dispatcher_update(const dp& other) : mstate(other.mstate) {}

bool dp::is_factorizable() const { return UPDATE_STYLE == FACTORIZED; }

consistency_model dp::consistency() const { return graphlab::DEFAULT_CONSISTENCY; }

consistency_model dp::gather_consistency() {
  
  process& p = process::get_process();
  
  // send invocation to child
  ji invocation("gather_consistency", mstate);
  p.send(invocation);
  
  // parse results
  jr result;
  p.receive(result);
  
  return result.consistency();
  
}

consistency_model dp::scatter_consistency() {
  
  process& p = process::get_process();
  
  // send invocation to child
  ji invocation("scatter_consistency", mstate);
  p.send(invocation);
  
  // parse results
  jr result;
  p.receive(result);
  
  return result.consistency();
  
}

edge_set dp::gather_edges() const {
  
  // invoke
  process& p = process::get_process();
  ji invocation("gather_edges", mstate);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  return result.edge_set();
  
}

edge_set dp::scatter_edges() const {
  
  // invoke
  process& p = process::get_process();
  ji invocation("scatter_edges", mstate);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  return result.edge_set();
  
}

void dp::init_gather(icontext_type& context) {

  // invoke
  process& p = process::get_process();
  ji invocation("init_gather", mstate);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);

  // update states - null means no change
  const char *cstring = result.updater();
  if (NULL != cstring)
    mstate = std::string(cstring);

}

void dp::gather(icontext_type& context, const edge_type& edge) {

  // invoke
  process& p = process::get_process();
  ji invocation("gather", mstate);
  invocation.add_context(context, ji::VERTEX);
  invocation.add_edge(context, edge);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  // update states - null means no change
  const char *cstring = result.updater();
  if (NULL != cstring)
    mstate = std::string(cstring);

}

void dp::apply(icontext_type& context) {
  
  // invoke
  process &p = process::get_process();
  ji invocation("apply", mstate);
  invocation.add_context(context, ji::VERTEX);
  p.send(invocation);
  
  // receive
  jr result;
  p.receive(result);
  
  // update states - null means no change
  const char *cstring = result.updater();
  if (NULL != cstring)
    mstate = std::string(cstring);                      // copy
  cstring = result.vertex();
  if (NULL != cstring)
    context.vertex_data().state = std::string(cstring); // copy
  
}

void dp::scatter(icontext_type& context, const edge_type& edge) {

  process& p = process::get_process();
  
  // invoke
  ji invocation("scatter", mstate);
  invocation.add_context(context, ji::VERTEX);
  invocation.add_edge(context, edge);
  p.send(invocation);
  
  // parse results
  jr result;
  p.receive(result);

  // TODO: edge state

  // update states - null means no schedule
  const char *cstring = result.scatter_schedule();
  if (NULL == cstring) return;

  const graph_type::vertex_id_type neighbor = (context.vertex_id() == edge.source()) ? edge.target() : edge.source();
    
  // schedule
  if (strcmp("self", cstring)) context.schedule(neighbor, dp(std::string(cstring)));
  else context.schedule(neighbor, *this);
  
}

inline void dp::operator+=(const dp& other){
  // TODO
}

void dp::operator()(icontext_type& context){
  
  process& p = process::get_process();
  
  // send invocation to child
  ji invocation("update", mstate);
  invocation.add_context(context, ji::VERTEX | ji::EDGES);
  p.send(invocation);
  
  // parse results
  jr result;
  p.receive(result);

  // update states - null means no change
  const char *cstring = result.updater();
  if (NULL != cstring)
    mstate = std::string(cstring);                      // copy
  cstring = result.vertex();
  if (NULL != cstring)
    context.vertex_data().state = std::string(cstring); // copy
  
  // TODO: edge states
    
  // schedule
  schedule(context, result);

}

void dp::schedule(icontext_type& context, const jr& result){

  const js schedule = result.schedule();
  
  switch (schedule.targets()){
    case js::IN_NEIGHBORS: {
      const dp &updater = (schedule.updater() == "self") ? *this : dp(schedule.updater());
      context.schedule_in_neighbors(context.vertex_id(), updater);
      break;
    }
    case js::OUT_NEIGHBORS: {
      const dp &updater = (schedule.updater() == "self") ? *this : dp(schedule.updater());
      context.schedule_out_neighbors(context.vertex_id(), updater);
      break;
    }
    case js::ALL: {
      const dp &updater = (schedule.updater() == "self") ? *this : dp(schedule.updater());
      context.schedule_neighbors(context.vertex_id(), updater);
      break;
    }
    case js::SOME: {
      foreach(vertex_id_type vertex, schedule.vertices()){
        const dp &updater = (schedule.updater() == "self") ? *this : dp(schedule.updater());
        context.schedule(vertex, updater);
      }
      break;
    }
    case js::NONE: {/* do nothing */}
  }
  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////// MAIN METHOD /////////////////////////////////
int main(int argc, char** argv) {

  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  
  // Parse command line options -----------------------------------------------
  command_line_options clopts("GraphLab Dispatcher");
  std::string graph_file = "toy.tsv";
  std::string updater;
  std::string format = "tsv";
  clopts.attach_option("graph", &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", &format, format,
                       "The graph file format: {metis, snap, tsv}");
  clopts.attach_option("updater", &updater,
                       "The updater to execute (required)");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  
  if(!clopts.is_set("updater")) {
    std::cout << "Updater not provided." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }

  // Setup the GraphLab execution core and load graph -------------------------
  core<dp::graph_type, dp> core;
  core.set_options(clopts);      // attach the command line options to the core
  std::cout << "Loading graph from file" << std::endl;
  const bool success = dispatcher_ops::load_structure(graph_file, core.graph());
  if(!success)
    std::cout << "Error in reading file: " << graph_file << std::endl;
  
  // Signal Handling ----------------------------------------------------------
  struct sigaction sa;
  std::memset(&sa, 0, sizeof(sa));
  sa.sa_handler = SIG_IGN;
  CHECK(!::sigaction(SIGPIPE, &sa, NULL));

  // Run the Dispatcher -------------------------------------------------------
  process::set_executable(updater);
  core.schedule_all(dp());
  std::cout << "Running dispatcher..." << std::endl;
  const double runtime = core.start();  // Run the engine
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  std::cout << "Update Rate (updates/second): " 
            << core.last_update_count() / runtime
            << std::endl;

  // Output Results -----------------------------------------------------------
  dispatcher_ops::save_structure("output.graph", core.graph());
  return EXIT_SUCCESS;
  
}
////////////////////////////////////////////////////////////////////////////////

#include <graphlab/macros_undef.hpp>
