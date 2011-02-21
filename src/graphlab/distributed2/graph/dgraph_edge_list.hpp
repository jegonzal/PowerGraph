#ifndef GRAPHLAB_DISTRIBUTED_DGRAPH_EDGE_LIST_HPP
#define GRAPHLAB_DISTRIBUTED_DGRAPH_EDGE_LIST_HPP
#include <vector>
#include <algorithm>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/unordered_map.hpp>
#include <graphlab/graph/graph.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

namespace dgraph_elist_impl {
  
edge_id_t eid_identity(edge_id_t eid) {
  return eid;
}


edge_id_t eid_transform(const std::vector<edge_id_t> *eid_adapter,
                        edge_id_t eid) {
  return (*eid_adapter)[eid];
}

}


  /** This class defines a set of edges */
class dgraph_edge_list {
public:
  typedef boost::function<edge_id_t (edge_id_t)> TransformType;
  typedef boost::transform_iterator<TransformType, const edge_id_t*> iterator;
  typedef boost::transform_iterator<TransformType, const edge_id_t*> const_iterator;
  typedef edge_id_t value_type;
private:
  // this class can wrap two methods of describing an edge list
  // method 1: a regular edge_list which requires some conversion
  // method 2: a vector

  // method 1:
  edge_list elist;
  const std::vector<edge_id_t> *eid_adapter_ptr;
  // method 2:
  std::vector<edge_id_t> edgeidvec;
  
  size_t numel;
public:


  // an edge list which wraps a regular elist. Method 1.
  inline dgraph_edge_list(edge_list elist,
                          const std::vector<edge_id_t> &eid_adapter)
                              :elist(elist), eid_adapter_ptr(&eid_adapter), numel(elist.size()){ }

  // an edge list which wraps a vector. Method 2.
  inline dgraph_edge_list(const std::vector<edge_id_t> &edgeidvec) :
                          eid_adapter_ptr(NULL), edgeidvec(edgeidvec), numel(edgeidvec.size()) { }

  /** \brief Get the size of the edge list */
  inline size_t size() const { return numel; }

  /** \brief Get the ith edge in the edge list */
  inline edge_id_t operator[](size_t i) const {
    assert(i < size());
    return *(begin() + i);
  }

  /** \brief Returns a pointer to the start of the edge list */
  inline iterator begin() const {
    if (eid_adapter_ptr != NULL) {
      return boost::make_transform_iterator(elist.begin(),
                                               boost::bind(dgraph_elist_impl::eid_transform,
                                                            eid_adapter_ptr, _1));

    }
    else {
      return boost::make_transform_iterator(&(edgeidvec[0]),
                                      dgraph_elist_impl::eid_identity);
    }
  }

  /** \brief Returns a pointer to the end of the edge list */
  inline iterator end() const {
    if (eid_adapter_ptr != NULL) {
      return boost::make_transform_iterator(elist.end(),
                                             boost::bind(dgraph_elist_impl::eid_transform,
                                                          eid_adapter_ptr, _1));

    }
    else {
      return boost::make_transform_iterator(&(edgeidvec[edgeidvec.size()]),
                                             dgraph_elist_impl::eid_identity);

    }
  }

  /** \brief Fill a vector with edge id list */
  inline void fill_vector(std::vector<edge_id_t>& lvalue) const {
    lvalue.clear();
    foreach(edge_id_t eid, *this) lvalue.push_back(eid);
  }

  /** \brief test if the edge list is empty */
  inline bool empty() const { return size() == 0; }
};

}
#include <graphlab/macros_undef.hpp>
#endif