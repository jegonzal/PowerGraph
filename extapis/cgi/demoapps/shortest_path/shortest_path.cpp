#include <graphlab.hpp>
#include <stdio.h>
#include <unistd.h>

#include <graphlab/macros_def.hpp>


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


int main(int argc, char** argv) {
  
  // loop until exit is received
  // while (true){
  //    read message
  //    break if exit
  //    process message, construct return
  // }

}







