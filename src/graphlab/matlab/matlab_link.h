// unfortunately this file has to be .h for emlc to recognize it
#ifndef EMLC_LINK_FUNCTIONS_H
#define EMLC_LINK_FUNCTIONS_H

#include "gl_emx_graphtypes.hpp"
/**
 * the matlab get_edge_data wrapper calls this function
 */
void emx_get_edge_data(HANDLE_TYPE handle, uint32_T eid, gl_emx_edgetype *edge);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_get_vertex_data(HANDLE_TYPE handle, uint32_T vid, gl_emx_vertextype *vertex);

/**
 * the matlab set_edge_data wrapper calls this function
 */
void emx_set_edge_data(HANDLE_TYPE handle, uint32_T eid, gl_emx_edgetype *edge);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_set_vertex_data(HANDLE_TYPE handle, uint32_T vid, gl_emx_vertextype *vertex);

/**
 * the matlab get_vertex_data wrapper calls this function
 */
void emx_set_vertex_data(HANDLE_TYPE handle, uint32_T vid, gl_emx_vertextype *vertex);

/**
 * The matlab add_task wrapper calls this function
 */
void emx_add_task(HANDLE_TYPE handle, uint32_T vid, const char* fnname, double priority);

/**
 * Draws a random integer value uniformly in a 32-bit range
 */
uint32_t emx_rand_int();

/**
 * Draws a random floating point value in the range [0,1)
 */
double emx_rand_double();

/**
 * Draws a gamma distributed random value
 */
double emx_rand_gamma(double alpha);

/**
 * Draws a random boolean which is true with probability p
 */
bool emx_rand_bernoulli(double p);

/**
 * Draws a random boolean which is true with probability p. Fast version
 */
bool emx_rand_bernoulli_fast(double p);

/**
 * Draws a random Gaussian distributed value
 */
double emx_rand_gaussian(double mean, double var);

/**
 * Draws a random integer in the range [1, high_inclusive]
 * equivalent to randi(high_inclusive)
 */
uint32_t emx_rand_int_uniform(const uint32_t high_inclusive);

/**
 * Draws a random integer in the range [1, high_inclusive]. Fast version.
 * equivalent to randi(high_inclusive).
 */
uint32_t emx_rand_int_uniform_fast(const uint32_t high_inclusive);

/**
 * Draws a random multinomial value. prob will be normalized automatically
 */
uint32_t emx_rand_multinomial(double* prob, uint32_t plength);



#endif