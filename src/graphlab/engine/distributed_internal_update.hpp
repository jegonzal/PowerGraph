#ifndef GRAPHLAB_DISTRIBUTED_INTERNAL_UPDATE_HPP
#define GRAPHLAB_DISTRIBUTED_INTERNAL_UPDATE_HPP
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
  
  
/**
 * When a vertex is scheduled, the following sequence is performed.
 * - A WAITING task is added
 * - An begin_task broadcast is sent to all mirrors
 * - Local Forks are requested.
 * - Once all forks have been acquired, a task with the LOCKS_ACQUIRED signal is added 
 * 
 * On an owned vertex, the task follows the following sequence
 * - WAITING: initialize with gather_and_apply_op = user task and count_down = #mirrors. 
 *             - This task must have priority = 0 and cannot be executed
 *               - WAITING  + TRANSMISSION = WAITING with decrement of count_down 
 *               - WAITING  + LOCKS_ACQUIRED = GATHERING with decrement of count_down 
 * - GATHERING: 
 *               - This task must have priority = 1
 *               - Execution of this task will
 *                    - Perform the local graph gather using gather_and_apply_op
 *                    - decrement the count_down variable
 *                    - Reschedule under the APPLYING state
 *               - GATHERING  + TRANSMISSION = GATHERING with decrement of count_down
 *               - GATHERING + any other state = failure
 * - APPLYING: - This task has priority = 1 if count_down = 0, and priority = 0 otherwise
 *             - Execution of this task will 
 *                    - Perform the local vertex apply using gather_and_apply_op
 *                    - Reschedule mirrors under MIRROR_SCATTERING state setting scatter_op
 *                    - Perform local scatter using gather_and_apply_op
 *                    - Release forks
 *             - APPLYING + TRANSMISSION = APPLYING with decrement of count_down
 *             - APPLYING + any other state = failure
 * 
 * 
 * On a mirrored vertex, on receiving a call to acquire forks
 * - add a MIRROR_WAITING task 
 * - Local Forks are requested.
 * - Once all forks have been acquired, a task with the LOCKS_ACQUIRED signal is added 
 * 
 * - MIRROR_WAITING: initialize with gather_and_apply_op = user task 
  *             - This task must have priority = 0 and cannot be executed
 *              - MIRROR_WAITING  + LOCKS_ACQUIRED = MIRROR_GATHERING with decrement of count_down 
 *              - MIRROR_WAITING + MIRROR_SCATTERING = MIRROR_WAITING | MIRROR_SCATTERING, joining tasks accordingly
 * - MIRROR_GATHERING: 
 *               - This task must have priority = 1
 *               - Execution of this task will
 *                    - Perform the local graph gather using gather_and_apply_op
 *                    - Reschedule owner vertex, with a TRANSMISSION state
 *               - MIRROR_GATHERING + any other state = failure 
 * - MIRROR_SCATTERING: scheduled by remote
 *               - This task must have priority = 1
 *               - Execution of this task will
 *                    - Perform the local graph scatter using scatter_op
  *                   - if MIRROR_WAITING, unset MIRROR_SCATTERING, clear scatter_op and reschedule current task
 *                    - Release locks
 *               - MIRROR_WAITING + MIRROR_SCATTERING = MIRROR_WAITING | MIRROR_SCATTERING, joining tasks accordingly
 *               - MIRROR_SCATTERING + any other state = failure
 *
 *   
 * Note that Only the MIRROR_WAITING and MIRROR_SCATTERING states can coexist
 */
template <typename Engine>
class distributed_internal_update {
  private:
    enum vertex_execution_state {
      NONE = 0,
      WAITING = 1,   // state on owner
      GATHERING = 2,   // state on owner
      APPLYING = 4,    // state on owner
      MIRROR_WAITING = 8
      MIRROR_GATHERING = 16, // state on mirror
      MIRROR_SCATTERING = 32, // state on mirror
    };
    
    enum vertex_signals {
      TRANSMISSION = 64
      LOCKS_ACQUIRED = 128
    };
    unsigned char state;
    bool has_gather_and_apply;
    bool has_scatter;
    typename Engine::update_functor_type gather_and_apply_op;
    typename Engine::update_functor_type scatter_op;
    size_t count_down;

};
  
}

#endif