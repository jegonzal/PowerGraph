#ifndef COLORING_H
#define COLORING_H

#include <iostream>
#include <graphlab.hpp>

//---------------- TYPES --------------------

/** Numerical representation for color. */
typedef unsigned long color_type;

/**
 * Vertex representation. Each vertex has a color.
 */
struct vertex_data {
  color_type color;
  vertex_data () : color(-1) {};
};

struct edge_data {
  // no edge data required
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;

//--------------- CONSTANTS ---------------
#define DEFAULT_FORMAT      "tsv"
#define NONE                ""

#define OPT_GRAPH_FILE      "graph"
#define OPT_FORMAT          "format"

//--------------- RETURN CODES ------------
#define ERR_INPUT -1
#define SUCCESS 0

#endif