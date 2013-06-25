#ifndef GRAPHLAB_WARP_GRAPH_MAP_REDUCE_HPP
#define GRAPHLAB_WARP_GRAPH_MAP_REDUCE_HPP

#include <boost/bind.hpp>
#include <graphlab/util/generics/conditional_combiner_wrapper.hpp>
#include <graphlab/parallel/fiber_group.hpp>
#include <graphlab/parallel/fiber_control.hpp>
#include <graphlab/parallel/fiber_remote_request.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

namespace warp {

namespace warp_impl {

/**
 * \internal
 * The default combiner used for combining mapped results from
 * map_reduce_neighborhood; merges self with other using operator +=. 
 */
template <typename T>
void default_combiner(T& self, const T& other) {
  self += other;
}


template <typename T, typename ExtraArgs>
void extended_default_combiner(T& self, const T& other, const ExtraArgs& unused) {
  self += other;
}


template <typename RetType, typename GraphType>
struct map_reduce_neighborhood_impl {

  typedef typename GraphType::vertex_type vertex_type;
  typedef typename GraphType::edge_type edge_type;
  typedef typename GraphType::local_vertex_type local_vertex_type;
  typedef typename GraphType::local_edge_type local_edge_type;
  typedef typename GraphType::vertex_record vertex_record;


/**************************************************************************/
/*                                                                        */
/*              Basic MapReduce Neighborhood Implementation               */
/*                                                                        */
/**************************************************************************/
/*
 * The master calls basic_mapreduce_neighborhood.
 * Which then issues calls to basic_local_mapper on each machine with a replica.
 */

  static conditional_combiner_wrapper<RetType> basic_local_mapper(GraphType& graph,
                                                           edge_dir_type edge_direction,
                                                           RetType (*mapper)(edge_type edge, vertex_type other),
                                                           void (*combiner)(RetType&, const RetType&),
                                                           vertex_id_type vid) {
    lvid_type lvid = graph.local_vid(vid);
    local_vertex_type local_vertex(graph.l_vertex(lvid));
    
    conditional_combiner_wrapper<RetType> accum(combiner);
    if(edge_direction == IN_EDGES || edge_direction == ALL_EDGES) {
      foreach(local_edge_type local_edge, local_vertex.in_edges()) {
        edge_type edge(local_edge);
        vertex_type other(local_edge.source());
        lvid_type a = edge.source().local_id(), b = edge.target().local_id();
        //vertexlocks[std::min(a,b)].lock();
        //vertexlocks[std::max(a,b)].lock();
        accum += mapper(edge, other);
        //vertexlocks[a].unlock();
        //vertexlocks[b].unlock();
      }
    } 
    // do out edges
    if(edge_direction == OUT_EDGES || edge_direction == ALL_EDGES) {
      foreach(local_edge_type local_edge, local_vertex.out_edges()) {
        edge_type edge(local_edge);
        vertex_type other(local_edge.target());
        lvid_type a = edge.source().local_id(), b = edge.target().local_id();
        //vertexlocks[std::min(a,b)].lock();
        //vertexlocks[std::max(a,b)].lock();
        accum += mapper(edge, other);
        //vertexlocks[a].unlock();
        //vertexlocks[b].unlock();
      }
    } 
    return accum;
  }


  static conditional_combiner_wrapper<RetType> basic_local_mapper_from_remote(size_t objid,
                                                           edge_dir_type edge_direction,
                                                           size_t mapper_ptr,
                                                           size_t combiner_ptr,
                                                           vertex_id_type vid) {
    // cast the mappers and combiners back into their pointer types
    RetType (*mapper)(edge_type edge, vertex_type other) = 
        reinterpret_cast<RetType(*)(edge_type, vertex_type)>(mapper_ptr);
    void (*combiner)(RetType&, const RetType&) = 
        reinterpret_cast<void (*)(RetType&, const RetType&)>(combiner_ptr);
    return basic_local_mapper(
        *reinterpret_cast<GraphType*>(distributed_control::get_instance()->get_registered_object(objid)),
        edge_direction,
        mapper,
        combiner,
        vid);
  }

  static RetType basic_map_reduce_neighborhood(typename GraphType::vertex_type current,
                                               edge_dir_type edge_direction,
                                               RetType (*mapper)(edge_type edge,
                                                                 vertex_type other),
                                               void (*combiner)(RetType& self, 
                                                                const RetType& other)) {
    // get a reference to the graph
    GraphType& graph = current.graph_ref;
    // get the object ID of the graph
    size_t objid = graph.get_rpc_obj_id();
    typename GraphType::vertex_record vrecord = graph.l_get_vertex_record(current.local_id());

    // make sure we are running on a master vertex
    ASSERT_EQ(vrecord.owner, distributed_control::get_instance_procid());
    
    // create num-mirrors worth of requests
    std::vector<request_future<conditional_combiner_wrapper<RetType> > > requests(vrecord.num_mirrors());
    
    size_t ctr = 0;
    foreach(procid_t proc, vrecord.mirrors()) {
      // issue the communication
        requests[ctr] = fiber_remote_request(proc, 
                                             map_reduce_neighborhood_impl<RetType, GraphType>::basic_local_mapper_from_remote,
                                             objid,
                                             edge_direction,
                                             reinterpret_cast<size_t>(mapper),
                                             reinterpret_cast<size_t>(combiner),
                                             current.id());
        ++ctr;
    }
    // compute the local tasks
    conditional_combiner_wrapper<RetType> accum = basic_local_mapper(graph, 
                                                                     edge_direction, 
                                                                     mapper, 
                                                                     combiner,
                                                                     current.id());
    accum.set_combiner(combiner);
    // now, wait for everyone
    for (size_t i = 0;i < requests.size(); ++i) {
      accum += requests[i]();
    }
    return accum.value;
  }



/**************************************************************************/
/*                                                                        */
/*           Extended MapReduce Neighborhood Implementation               */
/*                                                                        */
/**************************************************************************/
/*
 * The master calls extended_mapreduce_neighborhood.
 * Which then issues calls to extended_local_mapper on each machine with a replica.
 * The extended mapreduce neighborhood allows the mapper and combiner to take
 * an optional argument
 */




  template <typename ExtraArg>
  static conditional_combiner_wrapper<RetType> extended_local_mapper(GraphType& graph,
                                                              edge_dir_type edge_direction,
                                                              RetType (*mapper)(edge_type edge, vertex_type other, const ExtraArg&),
                                                              void (*combiner)(RetType&, const RetType&, const ExtraArg&),
                                                              vertex_id_type vid,
                                                              const ExtraArg& extra) {

    lvid_type lvid = graph.local_vid(vid);
    local_vertex_type local_vertex(graph.l_vertex(lvid));
    
    conditional_combiner_wrapper<RetType> accum(boost::bind(combiner, _1, _2, extra));
    if(edge_direction == IN_EDGES || edge_direction == ALL_EDGES) {
      foreach(local_edge_type local_edge, local_vertex.in_edges()) {
        edge_type edge(local_edge);
        vertex_type other(local_edge.source());
        lvid_type a = edge.source().local_id(), b = edge.target().local_id();
        //vertexlocks[std::min(a,b)].lock();
        //vertexlocks[std::max(a,b)].lock();
        accum += mapper(edge, other, extra);
        //vertexlocks[a].unlock();
        //vertexlocks[b].unlock();
      }
    } 
    // do out edges
    if(edge_direction == OUT_EDGES || edge_direction == ALL_EDGES) {
      foreach(local_edge_type local_edge, local_vertex.out_edges()) {
        edge_type edge(local_edge);
        vertex_type other(local_edge.target());
        lvid_type a = edge.source().local_id(), b = edge.target().local_id();
        //vertexlocks[std::min(a,b)].lock();
        //vertexlocks[std::max(a,b)].lock();
        accum += mapper(edge, other, extra);
        //vertexlocks[a].unlock();
        //vertexlocks[b].unlock();
      }
    } 
    return accum;
  }

  template <typename ExtraArg>
  static conditional_combiner_wrapper<RetType> extended_local_mapper_from_remote(size_t objid,
                                                              edge_dir_type edge_direction,
                                                              size_t mapper_ptr,
                                                              size_t combiner_ptr,
                                                              vertex_id_type vid,
                                                              const ExtraArg& extra) {
    // cast the mappers and combiners back into their pointer types
    RetType (*mapper)(edge_type edge, vertex_type other, const ExtraArg&) = 
        reinterpret_cast<RetType(edge_type, vertex_type)>(mapper_ptr);
    void (*combiner)(RetType&, const RetType&, const ExtraArg&) = 
        reinterpret_cast<void (RetType&, const RetType&)>(combiner_ptr);
    return extended_local_mapper<ExtraArg>(
        *reinterpret_cast<GraphType*>(distributed_control::get_instance()->get_registered_object(objid)),
        edge_direction,
        mapper,
        combiner,
        vid,
        extra);
  }

  template <typename ExtraArg>
  static RetType extended_map_reduce_neighborhood(typename GraphType::vertex_type current,
                                                  edge_dir_type edge_direction,
                                                  const ExtraArg& extra,
                                                  RetType (*mapper)(edge_type edge,
                                                                    vertex_type other,
                                                                    const ExtraArg& extra),
                                                  void (*combiner)(RetType& self, 
                                                                   const RetType& other,
                                                                   const ExtraArg& extra)) {
    // get a reference to the graph
    GraphType& graph = current.graph_ref;
    typename GraphType::vertex_record vrecord = graph.l_get_vertex_record(current.local_id());
    size_t objid = graph.get_rpc_obj_id();

    // make sure we are running on a master vertex
    ASSERT_EQ(vrecord.owner, distributed_control::get_instance_procid());
    
    // create num-mirrors worth of requests
    std::vector<request_future<conditional_combiner_wrapper<RetType> > > requests(vrecord.num_mirrors());
    
    size_t ctr = 0;
    foreach(procid_t proc, vrecord.mirrors()) {
      // issue the communication
      requests[ctr] = fiber_remote_request(proc, 
                                           map_reduce_neighborhood_impl<RetType, GraphType>::extended_local_mapper_from_remote<ExtraArg>,
                                           objid,
                                           edge_direction,
                                           reinterpret_cast<size_t>(mapper),
                                           reinterpret_cast<size_t>(combiner),
                                           current.id());
        ++ctr;
    }
    // compute the local tasks
    conditional_combiner_wrapper<RetType> accum = 
        extended_local_mapper<ExtraArg>(graph, edge_direction, mapper, 
                                        combiner, current.id(), extra);

    accum.set_combiner(boost::bind(combiner, _1, _2, boost::ref(extra)));
    // now, wait for everyone
    for (size_t i = 0;i < requests.size(); ++i) {
      accum += requests[i]();
    }
    return accum.value;
  }

};

} // namespace warp::warp_impl


template <typename RetType, typename VertexType>
RetType map_reduce_neighborhood(VertexType current,
                                edge_dir_type edge_direction,
                                RetType (*mapper)(typename VertexType::graph_type::edge_type edge,
                                                  VertexType other),
                                void (*combiner)(RetType& self, 
                                                 const RetType& other) = warp_impl::default_combiner<RetType>) {
  return warp_impl::
      map_reduce_neighborhood_impl<RetType, 
                                  typename VertexType::graph_type>::
                                      basic_map_reduce_neighborhood(current, edge_direction, 
                                                                    mapper, combiner);
}

template <typename RetType, typename ExtraArg, typename VertexType>
RetType map_reduce_neighborhood(VertexType current,
                                edge_dir_type edge_direction,
                                const ExtraArg& extra,
                                RetType (*mapper)(typename VertexType::graph_type::edge_type edge,
                                                  VertexType other,
                                                  const ExtraArg& extra),
                                void (*combiner)(RetType& self, 
                                                 const RetType& other,
                                                 const ExtraArg& extra) = warp_impl::extended_default_combiner<RetType, ExtraArg>) {
  return warp_impl::
      map_reduce_neighborhood_impl<RetType, 
                                  typename VertexType::graph_type>::template
                                      extended_map_reduce_neighborhood<ExtraArg>(current, edge_direction, 
                                                                                 extra, 
                                                                                 mapper, combiner);
}




} // namespace warp

} // namespace graphlab

#include <graphlab/macros_undef.hpp>
#endif
