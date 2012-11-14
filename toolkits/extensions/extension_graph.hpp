#ifndef GRAPHLAB_EXTENSION_GRAPH_HPP
#define GRAPHLAB_EXTENSION_GRAPH_HPP

#include "extension_data.hpp"
#include "extension_gas.hpp"
#include "extension_gas_lambda_wrapper.hpp"
#include "extension_main.hpp"
namespace graphlab {
namespace dc_impl {
extern distributed_control* get_last_dc();
}
namespace extension {

struct extension_graph_writer{
  std::string field;
  extension_graph_writer(std::string field):field(field) { }
  std::string save_vertex(internal_graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data()(field) << "\n";
    return strm.str();
  }

  std::string save_edge(internal_graph_type::edge_type e) {
    std::stringstream strm;
    strm << e.source().id() << "\t" << e.target().id() << e.data()(field) << "\n";
    return strm.str();
  }
}; 


class extension_graph {
 public:
  dc_dist_object<extension_graph> rmi;
  internal_graph_type internal_graph;
  mutex lock;
  bool finalized;

  extension_graph()
     :rmi(*dc_impl::get_last_dc(), this), 
     internal_graph(*dc_impl::get_last_dc(), __glopts),finalized(false) { }


  extension_graph(distributed_control& dc, 
                  const graphlab_options& opts = graphlab_options() ) 
     :rmi(dc, this), internal_graph(dc, opts),finalized(false) { }

  template <typename FieldType, typename TransformType>
  void transform_field(FieldType field,
                          TransformType transform_functor) {
    finalize();
    lock.lock();
    transform_field_wrapper<FieldType, TransformType> fw(&transform_functor, field);
    internal_graph.transform_vertices(fw);
    lock.unlock();
  }

/*
  template <typename FieldType, typename MapFunctor>
  var map_reduce_field(FieldType field,
                       MapFunctor map_functor) {
    finalize();
    lock.lock();
    map_field_wrapper<FieldType, MapFunctor> mw(&map_functor, field);
    var ret = internal_graph.map_reduce_vertices<var>(mw);
    lock.unlock();
    return ret;
  }
*/

  void load_structure(std::string prefix, std::string format) {
    lock.lock();
    internal_graph.load_format(prefix, format);
    lock.unlock();
  }

  void save_vertices(std::string prefix, std::string field) {
    internal_graph.save(prefix, extension_graph_writer(field),
                        false,    // do not gzip
                        true,     // save vertices
                        false);   // do not save edges
  }

  internal_graph_type& graph() {
    return internal_graph;
  }

  void finalize() {
    lock.lock();
    if (!finalized) {
      internal_graph.finalize();
      internal_graph.transform_vertices([] (internal_graph_type::vertex_type& v) {
                                        v.data().field("in_degree") = (double)v.num_in_edges();
                                        v.data().field("out_degree") = (double)v.num_out_edges();
                                        });
    }
    lock.unlock();
  }

  void synchronous_dispatch_new_engine(size_t desc_id);

  /// GAS which defaults to all out and all in edges
  template <typename GatherType,
           typename CombinerType,
           typename ApplyType,
           typename ScatterType>
    void GAS(GatherType gather,
             CombinerType combiner,
             ApplyType apply,
             ScatterType scatter,
             size_t iterations = 0) {
      finalize();
      generic_gather<GatherType> g;
      g.gt = &gather;

      generic_combiner<CombinerType> c;
      c.ct = &combiner;

      generic_apply<ApplyType> a;
      a.at = &apply;

      generic_scatter<ScatterType> s;
      s.st = &scatter;

      gas_op_descriptor gd;
      gd.gather_select_op = NULL;
      gd.gather_op = &g;
      gd.combiner_op = &c;
      gd.apply_op = &a;
      gd.scatter_select_op = NULL;
      gd.scatter_op = &s;

      lock.lock();
      descriptor_id_type descid = descriptor_access.push_back(gd);
      lock.unlock();

      synchronous_dispatch_new_engine(descid);
    }

  /// Regular GAS
  template <typename GatherSelectType,
           typename GatherType,
           typename CombinerType,
           typename ApplyType,
           typename ScatterSelectType,
           typename ScatterType>
    void GAS(GatherSelectType gatherselect,
             GatherType gather,
             CombinerType combiner,
             ApplyType apply,
             ScatterSelectType scatterselect,
             ScatterType scatter,
             size_t iterations = 0) {
      finalize();
      generic_gather_select<GatherSelectType> gs;
      gs.gt = &gatherselect;

      generic_gather<GatherType> g;
      g.gt = &gather;

      generic_combiner<CombinerType> c;
      c.ct = &combiner;

      generic_apply<ApplyType> a;
      a.at = &apply;

      generic_scatter_select<ScatterSelectType> ss;
      ss.st = &scatterselect;

      generic_scatter<ScatterType> s;
      s.st = &scatter;

      gas_op_descriptor gd;
      gd.gather_select_op = &gs;
      gd.gather_op = &g;
      gd.combiner_op = &c;
      gd.apply_op = &a;
      gd.scatter_select_op = &ss;
      gd.scatter_op = &s;

      lock.lock();
      descriptor_id_type descid = descriptor_access.push_back(gd);
      lock.unlock();

      synchronous_dispatch_new_engine(descid);
    }
};




} // namespace extension
} // namespace graphlab

#endif
