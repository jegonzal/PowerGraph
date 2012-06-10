#include <unistd.h>

#include "cgi.hpp"

#define INFINITY       1e+99

using namespace graphlab;
namespace json = rapidjson;

/**
 * Converts cstring to double. If cstring is empty, return INFINITY.
 */
static double stof(const char *str){
  if (0 == strlen(str)) return INFINITY;
  return atof(str);
}

/**
 * CGI Handler - contains methods (callbacks) to handle JSON invocations. Refer
 * to main GraphLab documentation for description of these methods.
 */
class cgi_handler : public icgi_handler {
public:
  
  void gather_edges(json::Document& invocation, json::Document& return_json){
    // we want to gather all incoming edges
    json::Document::AllocatorType& allocator = return_json.GetAllocator();
    return_json.AddMember("edges", "IN_EDGES", allocator);
  }
  
  void scatter_edges(json::Document& invocation, json::Document& return_json){
    // scatter to all outgoing edges
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
    if (best_dist < INFINITY && nbr_dist > (best_dist + edge_dist))
      return_json.AddMember("signal", "TARGET", allocator);
  }
  
};

int main(int argc, char** argv) {

  cgi_handler handler;
  graphlab::cgi client(handler); 
  client.listen();
  return 0;

}
