#include <graphlab.hpp>
#include "extension_graph.hpp"

namespace graphlab {
namespace extension {



void extension_graph::synchronous_dispatch_new_engine(size_t desc_id) {
  synchronous_engine<extension_update_functor> sync_engine(rmi.dc(), 
                                                           internal_graph,
                                                           __glopts);
  sync_engine.signal_all(desc_id);
  sync_engine.start();
}


} // extension
} // graphlab
