#include <graphlab.hpp>
#include <stdio.h>
#include <unistd.h>

#include <graphlab/macros_def.hpp>
#include "../../rapidjson.hpp"

#define INITIAL_LENGTH 256

namespace json = rapidjson;

double stof(const char *str){
  if (0 == strlen(str)) return 1e99;
  return atof(str);
}

void shortest_path_update(json::Document & invocation, json::Document& return_json){

  json::Document::AllocatorType& allocator = return_json.GetAllocator();
  return_json.SetObject();

  // TODO: error handling for missing elements
  const char *vertex_state = invocation["params"]["context"]["vertex"]["state"].GetString();
  double vertex_dist = stof(vertex_state);
  
  // relax all incoming edges
  json::Value& in_edges = invocation["params"]["context"]["in_edges"];
  for (json::SizeType i = 0; i < in_edges.Size(); i++){
    json::Value& edge = in_edges[i];
    double edge_dist = std::max(1.0, atof(edge["state"].GetString()));
    double nbr_dist = stof(edge["source"]["state"].GetString());
    vertex_dist = std::min(vertex_dist, nbr_dist + edge_dist);
  }
  
  std::stringstream doublestream;
  doublestream << vertex_dist;
  const char *doublestr = doublestream.str().c_str();
  
  // add shortest distance to return json
  return_json.AddMember("vertex", doublestr, allocator);

  // construct schedule
  json::Value schedule;
  schedule.SetObject();
  schedule.AddMember("updater", "self", allocator);
  
  json::Value vertices;
  vertices.SetArray();
  
  // schedule affected members
  json::Value& out_edges = invocation["params"]["context"]["out_edges"];
  for (json::SizeType i = 0; i < out_edges.Size(); i++){
  
    json::Value& edge = out_edges[i];
    int nbr_id = edge["target"]["id"].GetInt();
    double nbr_dist = stof(edge["target"]["state"].GetString());
    double edge_dist = std::max(1.0, atof(edge["state"].GetString()));
    
    if (nbr_dist > (vertex_dist + edge_dist))
      vertices.PushBack(nbr_id, allocator);
    
  }
  
  // add to return json
  schedule.AddMember("vertices", vertices, allocator);
  return_json.AddMember("schedule", schedule, allocator);

}

// not thread-safe
const char *handle_invocation(const char *buffer){

  if (NULL == buffer) return NULL;

  json::Document invocation;
  if (invocation.Parse<0>(buffer).HasParseError()){/* TODO: error handling */}
  
  // TODO: error handling for missing elements
  if (!strcmp(invocation["method"].GetString(), "exit")) return NULL;
  if (!strcmp(invocation["method"].GetString(), "update")){
  
    json::Document return_json;
    
    shortest_path_update(invocation, return_json);
    
    // this is not thread-safe
    static json::StringBuffer buffer;
    buffer.Clear();
    json::Writer<json::StringBuffer> writer(buffer);
    return_json.Accept(writer);
    
    return buffer.GetString();
    
  }

  // TODO: error handling - method not found?
  return NULL;
  
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
    const char *return_json = handle_invocation(buffer);
    if (!return_json) break;
    
    // return
    std::cout << strlen(return_json) << "\n";
    std::cout << return_json << std::flush;
    
  }
  
  delete[] buffer;
  return 0;

}







