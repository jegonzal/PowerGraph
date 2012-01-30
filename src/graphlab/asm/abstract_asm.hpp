#ifndef GRAPHLAB_ASYNC_STATE_MACHINE_HPP
#define GRAPHLAB_ASYNC_STATE_MACHINE_HPP

#include <graphlab/asm/async_state_machine_types.hpp>
#include <graphlab/parallel/atomic.hpp>

namespace graphlab {

/**
 Describes an asynchronous state machine.
 Must be copyable.
*/
template <typename ExecutorImplementation> 
class abstract_asm {
 public:
  typedef typename ExecutorImplementation::enumerator_type enumerator_type;
  typedef typename ExecutorImplementation::asm_callback_type asm_callback_type;
  typedef asm_state<enumerator_type> asm_state_type;
  
  virtual void init(asm_state<enumerator_type>& state, asm_callback_type& callback) = 0;
  virtual asm_action_args advance(asm_state<enumerator_type>& state, asm_callback_type& callback) = 0;
};


class abstract_asm_collection {
 public:
  virtual asm_machine_id_type add_state_machine(size_t identifier) = 0;
  virtual bool delete_machine(asm_machine_id_type machineid) = 0;  
  virtual bool advance(asm_machine_id_type machineid) = 0;  
};


/**
 Describes an asynchronous state machine.
 Must be copyable.
*/
template <typename ExecutorImplementation, typename StateMachine> 
class asm_collection : public abstract_asm_collection {
 public:
  typedef typename ExecutorImplementation::enumerator_type enumerator_type;
  typedef typename ExecutorImplementation::asm_callback_type asm_callback_type;
  typedef asm_state<enumerator_type> asm_state_type;
  
 private:
  ExecutorImplementation& impl;
  StateMachine statemachine;
  
  std::vector<mutex> lockset;
  std::vector<std::vector<std::pair<asm_machine_id_type, asm_state_type> > > stateset;
  
  atomic<asm_machine_id_type> next_machine_id;
 public:
  asm_collection(ExecutorImplementation& impl, 
                 size_t hashpoolsize = 100000):impl(impl), statemachine(statemachine) {
    lockset.resize(100000);
    stateset.resize(100000);
    impl.set_asm_collection(this);
    next_machine_id.value = 0;
  }
  
  asm_machine_id_type add_state_machine(size_t identifier) {
    asm_machine_id_type machineid = next_machine_id.inc_ret_last();
    size_t lockid = machineid % lockset.size();
    asm_state_type state;
    state.identifier = identifier;
    statemachine.init(state, impl.get_callback());
    lockset[lockid].lock();
    stateset[lockid].push_back(std::make_pair(machineid, state));
    lockset[lockid].unlock();
    return machineid;
  }

  bool get_state_of_machine(asm_machine_id_type machineid, asm_state_type &state) {
    bool ret = false;
    size_t lockid = machineid % lockset.size();
    lockset[lockid].lock();
    // look for the machine id
    for (size_t i = 0;i < stateset[lockid].size(); ++i) {
      if (stateset[lockid][i].first == machineid) {
        state = stateset[lockid][i].second;
        ret = true;
        break;
      }
    }
    lockset[lockid].unlock();
    return ret;
  }

  bool set_state_of_machine(asm_machine_id_type machineid, const asm_state_type &state) {
    bool ret = false;
    size_t lockid = machineid % lockset.size();
    lockset[lockid].lock();
    // look for the machine id
    for (size_t i = 0;i < stateset[lockid].size(); ++i) {
      if (stateset[lockid][i].first == machineid) {
        stateset[lockid][i].second = state;
        ret = true;
        break;
      }
    }
    lockset[lockid].unlock();
    return ret;
  }
  
  bool delete_machine(asm_machine_id_type machineid) {
    bool ret = false;
    size_t lockid = machineid % lockset.size();
    lockset[lockid].lock();
    // look for the machine id
    for (size_t i = 0;i < stateset[lockid].size(); ++i) {
      if (stateset[lockid][i].first == machineid) {
        // move the last to the front
        if (i < stateset[lockid].size() - 1) {
          std::swap(stateset[lockid][i], *stateset[lockid].rbegin());
        }
        // erase the last element
        stateset[lockid].erase(stateset[lockid].begin() + stateset[lockid].size() - 1);
        ret = true;
        break;
      }
    }
    lockset[lockid].unlock();
    return ret;
  }
  
  bool advance(asm_machine_id_type machineid) {
    asm_state_type curstate;
    if (get_state_of_machine(machineid,  curstate) == false) return false;
    asm_action_args actions = statemachine.advance(curstate, impl.get_callback());
    if (set_state_of_machine(machineid,  curstate) == false) return false;
    
    if (actions.action == ASM_ACTION_DELETE) {
      delete_machine(machineid);
    }
    else {
      impl.perform_delayed_action(machineid, actions);
    }
    return true;
  }
  
};

} // namespace graphlab

#endif
