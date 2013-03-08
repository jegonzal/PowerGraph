#include <graphlab/util/lock_free_pool.hpp>
#include <vector>

namespace graphlab {
namespace dc_impl {
struct pool_manager {
  lock_free_pool<oarchive> pool;

  pool_manager(){ pool.reset_pool(64); }

  ~pool_manager() {
    std::vector<oarchive>& r = pool.unsafe_get_pool_ref();
    for (size_t i = 0;i < r.size(); ++i) {
      if (r[i].buf) free(r[i].buf);
    }
  }
};

static pool_manager pool;

oarchive* oarchive_from_pool() {
  return pool.pool.alloc();
}

void release_oarchive_to_pool(oarchive* oarc) {
  oarc->off = 0;
  pool.pool.free(oarc);
}




} // dc_impl
} // graphlab
