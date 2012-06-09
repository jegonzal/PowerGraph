#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>
#include "../../rapidjson.hpp"

#define INITIAL_LENGTH 256

namespace json = rapidjson;

double stof(const char *str){
  if (0 == strlen(str)) return 1e+99;
  return atof(str);
}

void gather(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  // find distance from source to this vertex thru neighbor
  const json::Value& edge = invocation["params"]["edge"];
  double edge_dist = atof(edge["state"].GetString());
  double nbr_dist = stof(edge["source"]["state"].GetString());
  double vertex_dist = nbr_dist + edge_dist;
  // update updater state
  double best_dist = stof(invocation["state"].GetString());
  if (vertex_dist < best_dist){
    std::ostringstream double_stream;
    double_stream << vertex_dist;
    const char *double_str = double_stream.str().c_str();
    json::Value dist(double_str, allocator);
    return_json.AddMember("updater", dist, allocator);
  }
}

void merge(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  // keep the shortest distance
  double left = stof(invocation["state"].GetString());
  double right = stof(invocation["params"]["other"].GetString());
  double min_dist = std::min(left, right);
  // return
  std::ostringstream double_stream;
  double_stream << min_dist;
  const char *double_str = double_stream.str().c_str();
  json::Value dist(double_str, allocator);
  return_json.AddMember("updater", dist, allocator);
}

void apply(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  // save to vertex if it improves its distance
  const char *vertex_state = invocation["params"]["context"]["vertex"]["state"].GetString();
  double vertex_dist = stof(vertex_state);
  // update updater state
  double best_dist = stof(invocation["state"].GetString());
  if (best_dist < vertex_dist){
    std::ostringstream double_stream;
    double_stream << best_dist;
    const char *double_str = double_stream.str().c_str();
    json::Value dist(double_str, allocator);
    return_json.AddMember("vertex", dist, allocator);
  }
}

void scatter(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  double best_dist = stof(invocation["params"]["context"]["vertex"]["state"].GetString());
  // parse distances
  const json::Value& edge = invocation["params"]["edge"];
  double nbr_dist = stof(edge["target"]["state"].GetString());
  double edge_dist = atof(edge["state"].GetString());
  // relax edge
  if (nbr_dist > (best_dist + edge_dist)){
    return_json.AddMember("schedule", "self", allocator);
  }
}

void init_gather(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("updater", "", allocator);
}

void gather_edges(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("edges", "IN_EDGES", allocator);
}

void scatter_edges(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("edges", "OUT_EDGES", allocator);
}

void gather_consistency(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("consistency", "EDGE", allocator);
}

void scatter_consistency(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("consistency", "EDGE", allocator);
}

const char *handle_invocation(const char *buffer, json::StringBuffer& return_buffer){

  if (NULL == buffer) return NULL;

  json::Document invocation;
  if (invocation.Parse<0>(buffer).HasParseError()){/* TODO: error handling */}
  
  // TODO: error handling for missing elements
  if (!strcmp(invocation["method"].GetString(), "exit")) return NULL;
  
  json::Document return_json;
  return_json.SetObject();
  
  if (!strcmp(invocation["method"].GetString(), "gather"))
    gather(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "init_gather"))
    init_gather(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "gather_edges"))
    gather_edges(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "merge"))
    merge(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "apply"))
    apply(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "scatter"))
    scatter(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "scatter_edges"))
    scatter_edges(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "gather_consistency"))
    gather_consistency(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "scatter_consistency"))
    scatter_consistency(invocation, return_json);
  
  return_buffer.Clear();
  json::Writer<json::StringBuffer> writer(return_buffer);
  return_json.Accept(writer);
  return return_buffer.GetString();
  
}

int main(int argc, char** argv) {

  std::string line;
  std::size_t length = 0;
  std::size_t current_length = INITIAL_LENGTH;
  char *buffer = new char[current_length];
  
  // loop until exit is received
  while (true){
    
    // TODO: assume NULL-terminated string for now
    std::cin >> length;
    std::getline(std::cin, line);
    if (length + 1 > current_length){
      current_length = length + 1;
      delete[] buffer;
      buffer = new char[current_length];
    }
    
    // read message, break if exit
    std::cin.read(buffer, length);
    buffer[length] = NULL;    // terminate string w. null
    json::StringBuffer return_buffer;
    
    const char *return_json = handle_invocation(buffer, return_buffer);
    if (!return_json) break;
    
    // return
    std::cout << strlen(return_json) << "\n";
    std::cout << return_json << std::flush;
    
  }
  
  delete[] buffer;
  return 0;

}
