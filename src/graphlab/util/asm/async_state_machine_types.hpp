/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef GRAPHLAB_ASYNC_STATE_MACHINE_TYPES_HPP
#define GRAPHLAB_ASYNC_STATE_MACHINE_TYPES_HPP
#include <graphlab/graph/graph_basic_types.hpp>
namespace graphlab {
  typedef uint32_t asm_state_type;
  typedef uint32_t asm_action_type;
  typedef size_t asm_machine_id_type;
  template <typename Enumerator>
  struct asm_state {
    // internal state
    asm_state_type state;  // an arbitrary integer. uint64_t(-1) is reserved.
    // the meaning of this is entirely user-defined.
    size_t identifier; 
    Enumerator enumerator; // enumerator state
  };
  
  // RESERVED ACTIONS
#define ASM_ACTION_DELETE (asm_action_type)(-1)
  
  struct asm_action_args {
    asm_action_type action; // An action to perform to transit to the next state
    uint64_t action_arg1; // An argument for the action. Meaning depends on the
    // action to be executed
    uint64_t action_arg2; // 2ns argument for the action. Meaning depends on the
    asm_action_args() : action(0), action_arg1(0), action_arg2(0) { }
  };
}

#endif
