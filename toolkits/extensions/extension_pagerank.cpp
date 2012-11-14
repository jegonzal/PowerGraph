#include "extensions.hpp"

namespace graphlab {
namespace extension {

void pagerank(extension_graph& graph, 
              const std::string PR_FIELD,
              double tolerance) {
  const std::string PR_CHANGE = PR_FIELD + "_change";

  graph.transform_field(PR_FIELD, [](var v){ return 0.15; });    
  graph.GAS(
      [](const vars&) { return graphlab::IN_EDGES; },             // gather_edges
      [=](const vars&, vars&, const vars& other, edge_direction) { // gather
          return get<double>(other(PR_FIELD)) / 
                get<double>(other("out_degree")) ;
      }, 
      [](var& a, const var& b) {                                  // combine
          get<double>(a) += get<double>(b); 
      }, 
      [=](vars& v, const var& result) -> bool {                    // apply
          double pr = 0.15 + 0.85 * get<double>(result); 
          v(PR_CHANGE) = 
              std::fabs(pr - get<double>(v(PR_FIELD))) / 
              get<double>(v("out_degree"));
          v(PR_FIELD) = pr;         
          return false; 
      }, 
      [=](const vars& v) {                                        // scatter_edges
          return get<double>(v.field(PR_CHANGE)) > tolerance ? 
                                graphlab::OUT_EDGES : graphlab::NO_EDGES; 
      },
      [](const vars&, const vars&, const vars&, edge_direction) {// scatter 
          return true; 
      }
  ); // scatter
}


} // namespace extension
} // namespace graphlab

