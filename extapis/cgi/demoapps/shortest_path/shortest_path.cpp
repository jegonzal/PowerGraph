#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <sstream>

#include "../../rapidjson.hpp"

#define INITIAL_LENGTH 256
#define INFINITY       1e+99

namespace json = rapidjson;

double stof(const char *str){
  if (0 == strlen(str)) return INFINITY;
  return atof(str);
}

void init(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("program", "x", allocator);
}

void gather_edges(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("edges", "IN_EDGES", allocator);
}

void scatter_edges(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.AddMember("edges", "OUT_EDGES", allocator);
}

void gather(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  // find distance from source to this vertex thru neighbor
  const json::Value& edge = invocation["params"]["edge"];
  double edge_dist = atof(edge["state"].GetString());
  double nbr_dist = stof(edge["source"]["state"].GetString());
  double vertex_dist = nbr_dist + edge_dist;
  // return gather result
  std::ostringstream double_stream;
  double_stream << vertex_dist;
  const char *double_str = double_stream.str().c_str();
  json::Value dist(double_str, allocator);
  return_json.AddMember("result", dist, allocator);
}

void merge(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  // keep the shortest distance
  double left = stof(invocation["state"].GetString());
  double right = stof(invocation["params"]["other"].GetString());
  double min_dist = std::min(left, right);
  // return merged gather result
  std::ostringstream double_stream;
  double_stream << min_dist;
  const char *double_str = double_stream.str().c_str();
  json::Value dist(double_str, allocator);
  return_json.AddMember("result", dist, allocator);
}

void apply(json::Document& invocation, json::Document& return_json){
  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  // current vertex distance
  const char *vertex_state = invocation["params"]["vertex"]["state"].GetString();
  double vertex_dist = stof(vertex_state);
  // best distance found from gather
  double best_dist = stof(invocation["params"]["gather"].GetString());
  // if improvement, update
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
  // best distance found for this vertex
  double best_dist = stof(invocation["params"]["vertex"]["state"].GetString());
  // parse distances
  const json::Value& edge = invocation["params"]["edge"];
  double nbr_dist = stof(edge["target"]["state"].GetString());
  double edge_dist = atof(edge["state"].GetString());
  // relax edge
  if (best_dist < INFINITY && nbr_dist > (best_dist + edge_dist)){
    return_json.AddMember("signal", "TARGET", allocator);
  }
}

const char *handle_invocation(const char *buffer, json::StringBuffer& return_buffer){

  if (NULL == buffer) return NULL;

  json::Document invocation;
  if (invocation.Parse<0>(buffer).HasParseError()){/* TODO: error handling */}
  
  // TODO: error handling for missing elements
  if (!strcmp(invocation["method"].GetString(), "exit")) return NULL;
  
  json::Document return_json;
  return_json.SetObject();
  
  if (!strcmp(invocation["method"].GetString(), "init"))
    init(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "gather_edges"))
    gather_edges(invocation, return_json);
    
  if (!strcmp(invocation["method"].GetString(), "scatter_edges"))
    scatter_edges(invocation, return_json);
    
  if (!strcmp(invocation["method"].GetString(), "gather"))
    gather(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "merge"))
    merge(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "apply"))
    apply(invocation, return_json);
  
  if (!strcmp(invocation["method"].GetString(), "scatter"))
    scatter(invocation, return_json);
  
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
    
    // TODO: does this process know if the other end closes the pipe?
    
    // TODO: assume NULL-terminated string for now
    std::cin >> length; std::getline(std::cin, line);
    if (length + 1 > current_length){
      current_length = length + 1;
      delete[] buffer;
      buffer = new char[current_length];
    }
    
    // read message, break if exit
    std::cin.read(buffer, length);
    buffer[length] = '\0';    // terminate string w. null
    json::StringBuffer return_buffer;
    
    std::cerr << buffer << std::endl;
    const char *return_json = handle_invocation(buffer, return_buffer);
    if (!return_json) break;
    
    // return
    std::cout << strlen(return_json) << "\n";
    std::cout << return_json << std::flush;
    
  }
  
  delete[] buffer;
  return 0;

}
