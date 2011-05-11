#include <graphlab/distributed2/graph/dgraph_edge_list.hpp>

namespace graphlab {

namespace dgraph_elist_impl {

edge_id_t eid_identity(edge_id_t eid) {
  return eid;
}


} // namespace dgraph_elist_impl
} // namespace graphlab