#ifndef GRAPHLAB_EXTENSION_GAS_HPP
#define GRAPHLAB_EXTENSION_GAS_HPP

#include <vector>
#include <graphlab/vertex_program/ivertex_program.hpp>
#include "extension_data.hpp"
#include "extension_gas_base_types.hpp"

namespace graphlab {
namespace extension{

  /// A collection of all the user ops together
struct gas_op_descriptor{
  gather_functor* gather_op;
  gather_select_functor* gather_select_op;
  combiner_functor* combiner_op;
  apply_functor* apply_op;
  scatter_functor* scatter_op;
  scatter_select_functor* scatter_select_op;
};
// the active set of GAS sets to run
extern std::vector<gas_op_descriptor> descriptor_set;
extern lockfree_push_back<std::vector<gas_op_descriptor> > descriptor_access;
typedef uint32_t descriptor_id_type;



// A wrapper around the gather operation
struct gather_var {
  var v;
  descriptor_id_type descriptor_id;
  combiner_functor* combiner_op;

  gather_var():descriptor_id(-1), combiner_op(NULL) { }
  var& operator+=(const var& other) {
    (*combiner_op)(v, other);
    return v;
  }

  var& operator+=(const gather_var& other) {
    (*combiner_op)(v, other.v);
    return v;
  }


  inline void save(graphlab::oarchive& oarc) const {
    oarc << v << descriptor_id;
  }

  inline void load(graphlab::iarchive& iarc) {
    iarc >> v >> descriptor_id;
    gas_op_descriptor gas;
    bool ret = descriptor_access.query_unsafe(descriptor_id, gas);
    if (ret) combiner_op = gas.combiner_op;
  }
};



// a wrapper around the message
struct message_type: public graphlab::IS_POD_TYPE {
  descriptor_id_type descriptor;
  message_type(): descriptor(-1){ }
  message_type(descriptor_id_type d): descriptor(d) { };
  message_type& operator+=(const message_type& other) {
    return *this;
  }
};

struct extension_update_functor: 
    public graphlab::ivertex_program<internal_graph_type, gather_var, message_type> {
                                                                                  
public:
    typedef graphlab::ivertex_program<internal_graph_type, var> parent_type;

    descriptor_id_type descriptor_id;
//    gas_op_descriptor gas;

    extension_update_functor():descriptor_id(-1) {
    }
    
    extension_update_functor(size_t id):descriptor_id(id) { }

    inline void init(icontext_type& context,
                     const vertex_type& vertex, 
                     const message_type& msg) { 
      descriptor_id = msg.descriptor;
    }


    edge_dir_type gather_edges(icontext_type& context,
                               const vertex_type& vertex) const {
      gas_op_descriptor* gas = 
          descriptor_access.query_unsafe(descriptor_id);
      ASSERT_TRUE(gas != NULL);
      if (gas->gather_select_op) {
        return (*gas->gather_select_op)(vertex.data());
      } else {
        return ALL_EDGES;
      }
    }

    inline gather_var gather(icontext_type& context, 
                      const vertex_type& vertex,
                      edge_type& edge) const {
      gas_op_descriptor* gas = 
          descriptor_access.query_unsafe(descriptor_id);
      ASSERT_TRUE(gas != NULL);
      vertex_type other_vertex = edge.source().id() == vertex.id() ? 
                                    edge.target() : edge.source();
      gather_var ret;
      ret.v  = (*gas->gather_op)(vertex.data(), 
                                edge.data(), 
                                other_vertex.data(),
                                edge.source().id() == vertex.id() ? 
                                      OUT_EDGE : IN_EDGE);
      ret.descriptor_id = descriptor_id;
      ret.combiner_op = gas->combiner_op;
      return ret;
    }

    inline void apply(icontext_type& context, vertex_type& vertex,
                      const gather_type& total) {
      gas_op_descriptor* gas = 
          descriptor_access.query_unsafe(descriptor_id);
      ASSERT_TRUE(gas != NULL);
      bool sched = (*gas->apply_op)(vertex.data(), total.v);
      if (sched) context.signal(vertex, descriptor_id);
    }

    edge_dir_type scatter_edges(icontext_type& context,
                               const vertex_type& vertex) const {
      gas_op_descriptor* gas = 
          descriptor_access.query_unsafe(descriptor_id);
      ASSERT_TRUE(gas != NULL);
      if (gas->scatter_select_op) {
        return (*gas->scatter_select_op)(vertex.data());
      } else {
        return ALL_EDGES;
      }
    }

    inline void scatter(icontext_type& context, const vertex_type& vertex,
                 edge_type& edge) const {
      gas_op_descriptor* gas = 
          descriptor_access.query_unsafe(descriptor_id);
      ASSERT_TRUE(gas != NULL);
      vertex_type other_vertex = edge.source().id() == vertex.id() ? 
          edge.target() : edge.source();

      bool ret = (*gas->scatter_op)(vertex.data(), 
                                   edge.data(), 
                                   other_vertex.data(),
                                   edge.source().id() == vertex.id() ? 
                                   OUT_EDGE : IN_EDGE);
      if (ret) {
        context.signal(other_vertex, descriptor_id);
      }
    }

    inline void save(graphlab::oarchive& oarc) const {
      oarc << descriptor_id;
    }

    inline void load(graphlab::iarchive& iarc) {
      iarc >> descriptor_id;
    }
};





} // namespace extension
} // namespace graphlab

#endif
