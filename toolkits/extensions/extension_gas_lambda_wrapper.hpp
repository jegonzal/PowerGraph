#ifndef GRAPHLAB_EXTENSION_GAS_LAMBDA_WRAPPER
#define GRAPHLAB_EXTENSION_GAS_LAMBDA_WRAPPER


#include "extension_gas.hpp"
/*
  Implements a collection of lightweight generic wrappers around the GAS, 
  transform and map operations.
  These wrappers do nothing but store a functor which is compatible with the 
  base types in extension_gas_base_types.hpp, and call them.
*/

namespace graphlab {
namespace extension {


template <typename FieldType, typename Functor>
struct transform_field_wrapper {
  Functor* f;
  key_id_type field;
  typedef void result_type;
  
  transform_field_wrapper(Functor* f, FieldType field):
        f(f),field(get_id_from_name(field)) { }

  void operator()(internal_graph_type::vertex_type& vtx) {
    var& v = vtx.data()(field);
    v = (*f)(v);
  }
};


template <typename FieldType, typename Functor>
struct map_field_wrapper {
  Functor* f;
  key_id_type field;
  typedef var result_type;

  map_field_wrapper(Functor* f, FieldType field):
        f(f),field(get_id_from_name(field)) { }

  var operator()(internal_graph_type::vertex_type& vtx) {
    var& v = vtx.data()(field);
    return (*f)(v);
  }
};


// the combiner base type
template <typename CombinerType>
struct generic_combiner: public combiner_functor {
  CombinerType* ct;
  void operator()(var& a, const var& b) {
    (*ct)(a, b);
  }
};

template <typename GatherSelectType>
struct generic_gather_select : public gather_select_functor {
  GatherSelectType* gt;
  edge_dir_type operator()(const vars& center) {
    return (*gt)(center);
  }
};


template <typename GatherType>
struct generic_gather : public gather_functor {
  GatherType* gt;
  var operator()(const vars& center, 
                 vars& edge, 
                 const vars& other,
                 edge_direction direction) {
    return (*gt)(center, edge, other, direction);
  }
};

template <typename ApplyType>
struct generic_apply : public apply_functor {
  ApplyType* at;
  bool operator()(vars& center, 
                  const var& gather_result) {
    return (*at)(center, gather_result);
  }
};

template <typename ScatterSelectType>
struct generic_scatter_select : public scatter_select_functor {
  ScatterSelectType* st;
  edge_dir_type operator()(const vars& center) {
    return (*st)(center);
  }
};


template <typename ScatterType>
struct generic_scatter : public scatter_functor {
  ScatterType* st;
  bool operator()(const vars& center, 
                 vars& edge, 
                 const vars& other,
                 edge_direction direction) {
    return (*st)(center, edge, other, direction);
  }
};


} // extension
} // graphlab
#endif
