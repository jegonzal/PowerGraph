#ifndef GRAPHLAB_WARP_PARFOR_ALL_VERTICES_HPP
#define GRAPHLAB_WARP_PARFOR_ALL_VERTICES_HPP
#include <boost/function.hpp>
#include <graphlab/parallel/fiber_group.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/graph/vertex_set.hpp>
#include <graphlab/rpc/dc.hpp>
namespace graphlab {
namespace warp {

namespace warp_impl {


/*
 * Actual Parfor implementation.
 * Holds a reference to all the arguments.
 * Each fiber increments the atomic counter and runs the fn on it)
 */
template <typename GraphType>
struct parfor_all_vertices_impl{

  GraphType& graph; 
  boost::function<void(typename GraphType::vertex_type)> fn;
  vertex_set& vset;
  atomic<size_t> ctr;

  parfor_all_vertices_impl(GraphType& graph,
                           boost::function<void(typename GraphType::vertex_type)> fn,
                           vertex_set& vset): graph(graph),fn(fn),vset(vset),ctr(0) { }

  void run_fiber() {
    while (1) {
      size_t lvid = ctr.inc_ret_last();
      if (lvid >= graph.num_local_vertices()) break;
      typename GraphType::local_vertex_type l_vertex = graph.l_vertex(lvid);
      if (l_vertex.owned()) {
        typename GraphType::vertex_type vertex(l_vertex);
        fn(vertex);
      }
    } 
  }
};

} // namespace warp_impl

template <typename GraphType, typename FunctionType>
void parfor_all_vertices(GraphType& graph,
                         FunctionType fn,
                         vertex_set vset = GraphType::complete_set(),
                         size_t nfibers = 2000) {
  distributed_control::get_instance()->barrier();
  fiber_group group;
  warp_impl::parfor_all_vertices_impl<GraphType> parfor(graph, fn, vset);
  
  for (size_t i = 0;i < nfibers; ++i) {
    group.launch(boost::bind(&warp_impl::parfor_all_vertices_impl<GraphType>::run_fiber, &parfor));
  }
  group.join();
  distributed_control::get_instance()->barrier();
  graph.synchronize(vset);
}

} // namespace warp
} // namespace graphlab
#endif
