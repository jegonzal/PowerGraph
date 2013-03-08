#ifndef GRAPHLAB_RPC_ARCHIVE_MEMORY_POOL
#define GRAPHLAB_RPC_ARCHIVE_MEMORY_POOL
#include <graphlab/serialization/oarchive.hpp>
namespace graphlab {
namespace dc_impl {

oarchive* oarchive_from_pool();
void release_oarchive_to_pool(oarchive* oarc);

}
}
#endif
