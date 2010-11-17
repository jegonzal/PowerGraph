#ifndef GRAPHLAB_OPTIONS_STRUCT_H
#define GRAPHLAB_OPTIONS_STRUCT_H
#include "updates_types.h"
// unfortunately this file has to be a .h otherwise matlab won't recognize it

// To simplify the process of passing options to mex,
// I make use of my mex<->emx converters.
/**
 * Describes part of a schedule.
 * All vertices in [vertices] should be updated
 * according to "update_function" with priorities
 * in [priorities]
 * update_function: update function to use.
 * vertices: array of vertex ids to update
 * priority: same size as vertices. Priority for each update
 */
typedef struct graphlab_initial_schedule{
  emxArray_char_T update_function;
  emxArray_uint32_T vertices;
  emxArray_real_T priorities;
} graphlab_initial_schedule;


/**
 * Options which control the graphlab start up.
 * scheduler: Scheduler string
 * scope: Scope Type
 * ncpus: number of cpus to use
 * initial_schedule: and array of structs describing the schedule
 */
typedef struct graphlab_options_struct{
  emxArray_char_T *scheduler; 
  emxArray_char_T *scope_type;  
  real_T ncpus;
  graphlab_initial_schedule *initial_schedule;
} graphlab_options_struct;

#endif