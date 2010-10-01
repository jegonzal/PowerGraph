#ifndef JT_WORKER_HPP
#define JT_WORKER_HPP


#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cassert>



// Including Standard Libraries

#include <graphlab.hpp>


#include "data_structure.hpp"



#include <graphlab/macros_def.hpp>

class jt_worker : public graphlab::runnable {
public:
  
  typedef graphlab::general_scope_factory<mrf::graph_type>
  scope_factory_type;
  typedef scope_factor_type::iscope_type iscope_type;

private:
  scope_factory_type* scope_factory;
  size_t worker_id;
  bool active;

public:

  jt_worker() : scope_factory(NULL) { }

  void init(scope_factory_type* sf, size_t id) {
    scope_factory = sf;
    worker_id = id;
    active = true;
  }

  // get a root
  void run() {    
    jt::graph_type junction_tree;
  }


};  // End of JT work













#include <graphlab/macros_undef.hpp>
#endif
