#include <vector>
#include <graphlab/distributed2/distributed_glshared.hpp>

namespace graphlab {
namespace distgl_impl {
	

std::vector<distributed_glshared_base*>& get_global_dist_glshared_registry() {
  static std::vector<distributed_glshared_base*> global_registry;
  return global_registry;
}	

void register_dist_glshared(distributed_glshared_base* glsharedobj) {
  get_global_dist_glshared_registry().push_back(glsharedobj);
}



} // namespace distgl_impl
} // namespace graphlab
