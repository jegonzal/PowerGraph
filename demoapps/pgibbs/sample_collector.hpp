#ifndef SAMPLE_COLLECTOR_HPP
#define SAMPLE_COLLECTOR_HPP




/**
 * This saves the current sample to the appropriate consumer
 */
void sample_collector_apply(any& current_data, const any& param);

/**
 * merge function
 */ 
void sample_collector_merge(any& merge_dest, const any& merge_src);







#endif



