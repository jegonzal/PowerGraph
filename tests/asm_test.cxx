#include <graphlab/graph/graph.hpp>
#include <graphlab/asm/abstract_asm.hpp>

using namespace graphlab;

typedef graph<int, int> graph_type;

template <typename GraphType>
class executor { 
 private:
  abstract_asm_collection* collection;
  graph_type* graph;
 public:
  typedef std::pair<typename GraphType::edge_list_type*, 
                    typename GraphType::edge_list_type::const_iterator*> enumerator_type;
  typedef executor<GraphType> asm_callback_type;
  
  executor(graph_type* graph) :graph(graph) { };
  
  void set_asm_collection(abstract_asm_collection* coll) {
    collection = coll;
  }
  
  asm_callback_type& get_callback() {
    return *this;
  }
  
  void perform_delayed_action(asm_machine_id_type machineid,
                              asm_action_args action) {
//    std::cout << "Delaying: " << action.action << "\n";
  }
  
  
  bool ready(asm_machine_id_type machineid) {
    if (collection->advance(machineid)) {
//      std::cout << "Advancing: " << machineid << "\n";
      return true;
    }
    else {
//      std::cout << machineid << " no longer exists\n";
      return false;
    }
  }
  
  enumerator_type init_enumerator(vertex_id_type vid) {
    enumerator_type et;
    et.first = new typename GraphType::edge_list_type(graph->in_edges(vid));
    et.second = new typename GraphType::edge_list_type::const_iterator(et.first->begin());
    return et;
  }
  
  void delete_enumerator(enumerator_type enumerator) {
    delete enumerator.second;
    delete enumerator.first;
  }

};


class state_machine: public abstract_asm<executor<graph_type> > {
 public:
  void init(asm_state_type& state, asm_callback_type& callback) {
    state.state = 0;
  }
  
  asm_action_args advance(asm_state_type& state, asm_callback_type& callback) {
    asm_action_args action;
    if (state.state == 0) {
      state.enumerator = callback.init_enumerator(state.identifier);
      state.state++;
      action.action = 0;
    }
    else {
      if (*state.enumerator.second == state.enumerator.first->end()) {
        std::cout << "Done \n";
        callback.delete_enumerator(state.enumerator);
        action.action = ASM_ACTION_DELETE;
      }
      else {
        std::cout << (*state.enumerator.second)->source() << "\t";
        state.state++;
        (*state.enumerator.second)++;
        action.action = state.state;
      }
    }
    return action;
  }
};



class AsmTestSuite: public CxxTest::TestSuite {
 public:  
  void test_basic_asm() {
    graph_type g;
    g.add_vertex();   g.add_vertex();   g.add_vertex();
    g.add_edge(1, 0);
    g.add_edge(2, 0);
    g.add_edge(1, 2);
    g.finalize();
    executor<graph_type> exec(&g);
    
    asm_collection<executor<graph_type>, state_machine> coll(exec);
    asm_machine_id_type id[3];
    id[0] = coll.add_state_machine(0);
    id[1] = coll.add_state_machine(1);
    id[2] = coll.add_state_machine(2);
    std::cout << "\n";
    TS_ASSERT(exec.ready(id[0])); // init
    TS_ASSERT(exec.ready(id[0])); // print 1
    TS_ASSERT(exec.ready(id[0])); // print 2
    TS_ASSERT(exec.ready(id[0])); // finish and delete
    TS_ASSERT(!exec.ready(id[0])); // fail
    
    TS_ASSERT(exec.ready(id[1])); // init
    TS_ASSERT(exec.ready(id[1])); // finish and delete
    TS_ASSERT(!exec.ready(id[1])); // fail
    
    TS_ASSERT(exec.ready(id[2])); // init
    TS_ASSERT(exec.ready(id[2])); // print 1
    TS_ASSERT(exec.ready(id[2])); // finish and delete
    TS_ASSERT(!exec.ready(id[2])); // fail
  }
};
