#include <graphlab/util/lock_free_pool.hpp>
#include <graphlab/rpc/dc_compile_parameters.hpp>
#include <vector>

namespace graphlab {
namespace dc_impl {
struct pool_manager {
  lock_free_pool<oarchive> pool;

  pool_manager(){
    pool.reset_pool(64);
    std::vector<oarchive>& r = pool.unsafe_get_pool_ref();
    for (size_t i = 0;i < r.size(); ++i) {
      r[i].buf = (char*)malloc(BUFFER_RELINQUISH_LIMIT);
      r[i].off = 0;
      r[i].len = BUFFER_RELINQUISH_LIMIT;
    }
  }

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
  if (oarc->buf == NULL) {
    oarc->buf = (char*)malloc(BUFFER_RELINQUISH_LIMIT);
    oarc->len = BUFFER_RELINQUISH_LIMIT;
  }
  pool.pool.free(oarc);
}




} // dc_impl
} // graphlab
