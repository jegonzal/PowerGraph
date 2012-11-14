#include "extension_graph.hpp"


#ifndef GRAPHLAB_EXTENSIONS_HPP
#define GRAPHLAB_EXTENSIONS_HPP
#if 1
// we have to use this unreliable hack
// don't seem to be able to get -wrap main
// working. TOFIX
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
