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
  };
}

#endif
