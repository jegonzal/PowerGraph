#include <cmath>
#include <iostream>
#include "graph_typedefs.gen"

using namespace graphlab;

extern "C" void user_program(user_funs* f) {
    std::cout << "executing user program" << std::endl;
    float prev = f->get_vertex_data();
    float neighbors = f->reduce_neighbors(IN_EDGES);
    float curr = 0.85*neighbors + 0.15;
    f->set_vertex_data(curr);
    float last_change = std::abs(curr - prev);
    if (last_change > 0.01) {
	f->signal_neighbors(OUT_EDGES);
	std::cout << "signaling neighbors" << std::endl;
    }
}

extern "C" void vertex_reduce(vertex_data_type& a, const vertex_data_type& b) {
    a += b;
}

