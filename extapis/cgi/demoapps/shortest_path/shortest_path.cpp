#include <graphlab.hpp>
#include <stdio.h>
#include <unistd.h>

#include <graphlab/macros_def.hpp>
#include "../../rapidjson.hpp"

#define INITIAL_LENGTH 256

// class shortest_path_update : public
// graphlab::iupdate_functor<graph_type, shortest_path_update> {
// public:
//   void operator()(icontext_type& context) {
//     vertex_data& vdata = context.vertex_data();    
//     foreach(edge_type edge, context.in_edges()) {
//       const vertex_id_type nbr_id = edge.source();
//       const vertex_data& nbr = context.const_vertex_data(nbr_id);
//       const edge_data& edata = context.const_edge_data(edge);
//       vdata.dist = std::min(vdata.dist, nbr.dist + edata.dist);
//     }
//     // Reschedule any affected neighbors
//     foreach(edge_type edge, context.out_edges()) {
//       const vertex_id_type nbr_id = edge.target();
//       const vertex_data& nbr = context.const_vertex_data(nbr_id);
//       const edge_data& edata = context.const_edge_data(edge);
//       if(nbr.dist > (vdata.dist + edata.dist)) 
//         context.schedule(nbr_id, *this);    
//     }
//   }
// }; // end of shortest path update functor

// not thread-safe
const char *handle_invocation(const char *buffer, std::size_t length){

  if (NULL == buffer || 0 == length) return NULL;

  rapidjson::Document document;
  if (document.Parse<0>(data).HasParseError()){/* TODO: error handling */}
  
  if (!strcmp(document["method"], "exit")) return NULL;
  if (!strcmp(document["method"], "update")){
  
    rapidjson::Document return_json;
    
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    return_json.Accept(writer);
    
    return buffer.GetString();
    
  }

  return NULL;
  
}

int main(int argc, char** argv) {

  std::string line;
  std::size_t length = 0;
  std::size_t current_length = INITIAL_LENGTH;
  char *buffer = new char[current_length];
  
  // loop until exit is received
  while (true){
    
    std::cin >> length;
    std::getline(std::cin, first_line);
    if (length > current_length){
      current_length = length;
      delete[] buffer;
      buffer = new char[current_length];
    }
    
    // read message, break if exit
    std::cin.read(buffer, length);
    const char *return_json = handle_invocation(buffer, length);
    if (!return_json) break;
    
    // return
    std::cout << strlen(return_json) << "\n";
    std::cout << return_json << std::flush;
    
  }
  
  delete[] buffer;
  return 0;

}







