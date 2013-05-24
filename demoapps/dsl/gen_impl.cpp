#include <cmath>
#include <iostream>
#include "graph_typedefs.gen"

using namespace graphlab;

#define ALPHA 0.87

extern "C" void user_program(user_funs* f) {
    float prev = f->get_vertex_data();
    float neighbors = f->reduce_neighbors(IN_EDGES);
    //neighbors = neighbors/out_edges
    float curr = ALPHA*neighbors + (1-ALPHA);
    f->set_vertex_data(curr);
    float last_change = std::abs(curr - prev);
    std::cout << "last change was: " << last_change << std::endl;
    if (last_change > 0.01) {
	f->signal_neighbors(OUT_EDGES);
    }
    std::cout.flush();
}

extern "C" void vertex_reduce(vertex_data_type& a, const vertex_data_type& b) {
    a += b;
}

