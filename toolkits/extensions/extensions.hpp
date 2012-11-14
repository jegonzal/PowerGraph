#include "extension_graph.hpp"


#ifndef GRAPHLAB_EXTENSIONS_HPP
#define GRAPHLAB_EXTENSIONS_HPP
#ifdef __APPLE__
// annoyingly there is no main
// we have to use this unreliable hack
#define main __real_main
#endif


namespace graphlab {
namespace extension {

// prototype for all implemented extensions
void pagerank(extension_graph& graph, 
              const std::string PR_FIELD,
              double tolerance);

} // namespace extension
} // namespace graphlab




#endif 
