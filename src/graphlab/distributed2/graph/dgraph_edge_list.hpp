/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


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
  
edge_id_t eid_identity(edge_id_t eid);

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
  // method 2:
  std::vector<edge_id_t> edgeidvec;
  size_t numel;
  bool useelist;
  
public:


  // an edge list which wraps a regular elist. Method 1.
  inline dgraph_edge_list(edge_list elist)
                              :elist(elist), numel(elist.size()), useelist(true){ }

  // an edge list which wraps a vector. Method 2.
  inline dgraph_edge_list(const std::vector<edge_id_t> &edgeidvec) :
                          edgeidvec(edgeidvec), numel(edgeidvec.size()), useelist(false) { }

  /** \brief Get the size of the edge list */
  inline size_t size() const { return numel; }

  /** \brief Get the ith edge in the edge list */
  inline edge_id_t operator[](size_t i) const {
    assert(i < size());
    return *(begin() + i);
  }

  /** \brief Returns a pointer to the start of the edge list */
  inline iterator begin() const {
    if (useelist) {
      return boost::make_transform_iterator(elist.begin(),
                                            dgraph_elist_impl::eid_identity);

    }
    else {
      return boost::make_transform_iterator(&(edgeidvec[0]),
                                      dgraph_elist_impl::eid_identity);
    }
  }

  /** \brief Returns a pointer to the end of the edge list */
  inline iterator end() const {
    if (useelist) {
      return boost::make_transform_iterator(elist.end(),
                                             dgraph_elist_impl::eid_identity);

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

