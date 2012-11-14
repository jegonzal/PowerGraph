#include "extensions.hpp"

namespace graphlab {
namespace extension {

void pagerank(extension_graph& graph, 
              const std::string PR_FIELD_NAME,
              double tolerance) {
  const std::string PR_CHANGE_NAME = PR_FIELD_NAME + "_change";

  key_id_type PR_FIELD = get_id_from_name(PR_FIELD_NAME);
  key_id_type PR_CHANGE = get_id_from_name(PR_CHANGE_NAME);
  key_id_type OUT_DEG = get_id_from_name("out_degree");

  graph.transform_field(PR_FIELD, [](var v){ return 0.15; });    
  timer ti;
  graph.GAS(
      [](const vars&) { return graphlab::IN_EDGES; },             // gather_edges
      [=](const vars&, vars&, const vars& other, edge_direction) { // gather
          return get<double>(other(PR_FIELD)) / 
                get<double>(other(OUT_DEG)) ;
      }, 
      [](var& a, const var& b) {                                  // combine
          get<double>(a) += get<double>(b); 
      }, 
      [=](vars& v, const var& result) -> bool {                    // apply
          double pr = 0.15 + 0.85 * get<double>(result); 
          v(PR_CHANGE) = 
              std::fabs(pr - get<double>(v(PR_FIELD))) / 
              get<double>(v(OUT_DEG));
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
  std::cout << "PageRank complete in " << ti.current_time() << "s" << std::endl;
}


} // namespace extension
} // namespace graphlab

