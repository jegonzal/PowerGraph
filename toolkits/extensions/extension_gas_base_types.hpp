#ifndef GRAPHLAB_EXTENSION_GAS_BASE_TYPES_HPP
#define GRAPHLAB_EXTENSION_GAS_BASE_TYPES_HPP
#include "extension_data.hpp"

namespace graphlab {
namespace extension {
// edge direction enum
enum edge_direction {
  IN_EDGE, OUT_EDGE
};


//////////////////////////////////////////////////////////////
// Here we will produce base types of all the functors      //
//////////////////////////////////////////////////////////////

// the combiner base type
struct combiner_functor {
  virtual void operator()(var& , const var&) = 0;
};


// the gather base type
struct gather_functor {
  virtual var operator()(const vars& center, 
                         vars& edge, 
                         const vars& other,
                         edge_direction direction) = 0;
};

//gather select base type
struct gather_select_functor {
  virtual edge_dir_type operator()(const vars& center) {
    return ALL_EDGES;
  }
};

// the apply base type
struct apply_functor {
  // return true to schedule self
  virtual bool operator()(vars& center, 
                          const var& gather_result) = 0;
};

// scatter select base type
struct scatter_select_functor {
  virtual edge_dir_type operator()(const vars& center) {
    return ALL_EDGES;
  }
};


// the scatter base type
struct scatter_functor {
  // return true to schedule other
  virtual bool operator()(const vars& center, 
                          vars& edge, 
                          const vars& other,
                          edge_direction direction) = 0;
};



} // extension
} // graphlab 
#endif 
