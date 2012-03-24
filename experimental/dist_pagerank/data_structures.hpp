#ifndef PAGERANK_DATA_STRUCTURES
#define PAGERANK_DATA_STRUCTURES

#include <cstring>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <sstream>


#include <graphlab.hpp>




double termination_bound = 1e-5;
double random_reset_prob = 0.15;   // PageRank random reset probability


/**
 * Stores the value and the self weight
 */
struct vertex_data_type {
  float value;
  // we save weight of vertex's self-edge in the vertex data
  vertex_data_type(float value = 1) : value(value) { }
}; // End of vertex data
SERIALIZABLE_POD(vertex_data_type);

typedef char edge_data_type;




typedef graphlab::core<vertex_data_type, edge_data_type> core_type;
typedef core_type::types gl_types;
typedef core_type::types::graph local_graph_type;





#endif
