/**  
 * Copyright (c) 2011 Carnegie Mellon University. 
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

/* *
 * Author: Haijie Gu (haijieg@cs.cmu.edu)
 * Date: 11/10/2011
 *
 * CSR+CSC implementation of a graph storage.
 * */

#ifndef GRAPHLAB_GRAPH_STORAGE_HPP
#define GRAPHLAB_GRAPH_STORAGE_HPP


#ifndef __NO_OPENMP__
#include <omp.h>
#endif


#include <cmath>

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>

#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>

#include <boost/version.hpp>
#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>
#include <boost/iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/util/generics/shuffle.hpp>
#include <graphlab/graph/graph_basic_types.hpp>


#include <graphlab/parallel/atomic.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {


 
  template<typename VertexData, typename EdgeData>
  class json_parser;



  template<typename VertexData, typename EdgeData>
  class graph_storage {
  public:
    typedef graphlab::lvid_type lvid_type;
    typedef graphlab::edge_id_type edge_id_type;

    /** The type of the edge data stored in the graph. */
    typedef EdgeData edge_data_type;

    /** The type of the vertex data stored in the graph. */
    typedef VertexData vertex_data_type;

    /** \internal
     * \brief The type representing a range of edges. Edges
     * have consecutive internal ids.
     * */
    typedef std::pair<size_t, size_t>  edge_range_type;

    friend class json_parser<VertexData, EdgeData>;


    /* ----------------------------------------------------------------------------- */
    /* helper data field and structures: edge_data_list, class edge, class edge_list */
    /* ----------------------------------------------------------------------------- */
  public:
    // Edge class for temporary storage. Will be finalized into the CSR+CSC form.
    class edge_info {
    public:
      std::vector<EdgeData> data;
      std::vector<lvid_type> source_arr;
      std::vector<lvid_type> target_arr;
    public:
      edge_info () {}
      void reserve_edge_space(size_t n) {
        data.reserve(n);
        source_arr.reserve(n);
        target_arr.reserve(n);
      }
      // \brief Add an edge to the temporary storage.
      void add_edge(lvid_type source, lvid_type target, EdgeData _data) {
        data.push_back(_data);
        source_arr.push_back(source);
        target_arr.push_back(target);
      }
      // \brief Add edges in block to the temporary storage.
      void add_block_edges(const std::vector<lvid_type>& src_arr, 
                           const std::vector<lvid_type>& dst_arr, 
                           const std::vector<EdgeData>& edata_arr) {
        data.insert(data.end(), edata_arr.begin(), edata_arr.end());
        source_arr.insert(source_arr.end(), src_arr.begin(), src_arr.end());
        target_arr.insert(target_arr.end(), dst_arr.begin(), dst_arr.end());
      }
      // \brief Remove all contents in the storage. 
      void clear() {
        std::vector<EdgeData>().swap(data);
        std::vector<lvid_type>().swap(source_arr);
        std::vector<lvid_type>().swap(target_arr);
      }
      // \brief Return the size of the storage.
      size_t size() const {
        return source_arr.size();
      }
      // \brief Return the estimated memory footprint used.
      size_t estimate_sizeof() const {
        return data.capacity()*sizeof(EdgeData) + 
          source_arr.capacity()*sizeof(lvid_type)*2 + 
          sizeof(data) + sizeof(source_arr)*2 + sizeof(edge_info);
      }
    }; // end of class edge_info.

    // A class of edge information. Used as value type of the edge_list.
    class edge_type {
    public:
      /** \brief Creates an empty edge type. */
      edge_type () : _source(-1), _target(-1), _edge_id(-1), 
                     _dir(NO_EDGES), _empty(true) { }
      /** \brief Creates an edge of given source id, target id, edge
       * id and direction enum.  \internal edge id is used to locate
       * the edge data in the edge_data array.  edge_dir type is
       * defined in graph_basic.hpp. **/
      edge_type (const lvid_type _source, 
                 const lvid_type _target, 
                 const edge_id_type _eid, edge_dir_type _dir) :
        _source(_source), _target(_target), _edge_id(_eid), 
        _dir(_dir), _empty(false) {
          if (_dir != OUT_EDGES) std::swap(this->_source, this->_target);
        }
    public:
      /** \brief Returns the source vertex id of the edge. */
      inline lvid_type source() const {
        return _source;
      }
      /** \brief Returns the target vertex id of the edge. */
      inline lvid_type target() const { 
        return _target;
      }
      /** \brief Returns the direction of the edge. */
      inline edge_dir_type get_dir() const {
        return _dir;
      }
      /** \brief Returns whether this is an empty edge. */
      inline bool empty() const { return _empty; }
      // Data fields. 
    private:
      lvid_type _source;
      lvid_type _target;
      edge_id_type _edge_id;
      edge_dir_type _dir;
      bool _empty;

      friend class graph_storage;
    }; // end of class edge_type.

    // Internal iterator on edge_types.
    class edge_iterator  {
    public:
      typedef std::random_access_iterator_tag iterator_category;
      typedef edge_type    value_type;
      typedef ssize_t      difference_type;
      typedef edge_type*   pointer;
      typedef edge_type   reference;

    public:
      // Cosntructors
      /** \brief Creates an empty iterator. */
      edge_iterator () : offset(-1), empty(true) { }
      /** \brief Creates an iterator at a specific edge.
       * The edge location is defined by the follows: 
       * A center vertex id,  an offset to the center, the direction and the
       * pointer to the start of edge id array. */
      edge_iterator (lvid_type _center, size_t _offset, 
                     edge_dir_type _itype, const lvid_type* _vid_arr) :
        center(_center), offset(_offset), itype(_itype), vid_arr(_vid_arr), 
        empty(false) { }
      /** \brief Returns the value of the iterator. An empty iterator always returns empty edge type*/ 
      inline edge_type operator*() const  {
        //  ASSERT_TRUE(!empty);
        return make_value();
      }
#if BOOST_VERSION < 105000
      typedef boost::detail::
      operator_arrow_result<edge_type, edge_type, edge_type*> arrow_type;
      inline typename arrow_type::type operator->() const {
        return arrow_type::make(make_value());
      }
#else
      typedef typename boost::detail::
      operator_arrow_dispatch<edge_type, edge_type*>::result_type arrow_type;
      inline arrow_type operator->() const {
        return arrow_type(make_value());
      }
#endif

      /** \brief Returns if two iterators point to the same edge. */
      inline bool operator==(const edge_iterator& it) const {
        return (empty && it.empty) || 
          (empty == it.empty && itype == it.itype && center == it.center && 
           offset == it.offset);
      }

      /** \brief Returns if two iterators don't point to the same edge. */
      inline bool operator!=(const edge_iterator& it) const { 
        return !(*this == it);
      }

      /** \brief Increases the iterator. */
      inline edge_iterator& operator++() {
        //ASSERT_TRUE(!empty);
        ++offset;
        return *this;
      }

      /** \brief Increases the iterator. */
      inline edge_iterator operator++(int) {
        //ASSERT_TRUE(!empty);
        const edge_iterator copy(*this);
        operator++();
        return copy;
      }

      /** \brief Computes the difference of two iterators. */
      inline ssize_t operator-(const edge_iterator& it) const {
        return offset - it.offset;
      }

      /** \brief Returns a new iterator whose value is increased by i difference units. */
      inline edge_iterator operator+(difference_type i) const {
        return edge_iterator(center, offset+i, itype, vid_arr);
      }

      /** \brief Increases the iterator by i difference units. */
      inline edge_iterator& operator+=(difference_type i) {
        offset+=i;
        return *this;
      }

      /** \brief Generate the return value of the iterator. */
      inline edge_type make_value() const {
        return empty ? edge_type() : edge_type(center, vid_arr[offset], offset, itype);
      }

    private:
      lvid_type center;
      size_t offset;
      edge_dir_type itype;
      const lvid_type* vid_arr;
      bool empty;
    }; // end of class edge_iterator.

    /** Represents an iteratable list of edge_types. */
    class edge_list {
    public:
      typedef edge_iterator iterator;
      typedef edge_iterator const_iterator;
      typedef edge_type value_type;
    private:
      edge_iterator begin_iter, end_iter;
    public:
      /** Cosntructs an edge_list with begin and end.  */
      edge_list(const edge_iterator begin_iter = edge_iterator(), 
                const edge_iterator end_iter = edge_iterator()) : 
        begin_iter(begin_iter), end_iter(end_iter) { }
      inline size_t size() const { return end_iter - begin_iter;}            
      inline edge_type operator[](size_t i) const {return *(begin_iter + i);}
      iterator begin() const { return begin_iter; }
      iterator end() const { return end_iter; }
      bool empty() const { return size() == 0; }
    }; // end of class edge_list.


  public:
    // CONSTRUCTORS ============================================================>
    graph_storage() : use_skip_list(false) {  }

    // METHODS =================================================================>
   
    /** \internal \brief  Set graph storage to use skip list.
     * Skip list is used to jump between ...
     * */ 
    void set_use_skip_list (bool x) { use_skip_list = x;}

    /** \brief Returns the number of edges in the graph. */
    size_t edge_size() const { return num_edges; }

    /** \brief Returns the number of vertices in the graph. */
    size_t vertices_size() const { return num_vertices; }

    /** \brief Returns the number of in edges of the vertex. */
    size_t num_in_edges (const lvid_type v) const {
      if (v >= num_vertices)
        return 0;

      ASSERT_LT(v, num_vertices);

      size_t begin = CSC_dst[v];
      if (begin >= num_edges) return 0;
      // Search is the next valid src vertex after v.
      size_t search = use_skip_list ? 
        nextValid(CSC_dst_skip, v, true) : nextValid(CSC_dst, v, false);
      size_t end = (search >= num_vertices) ? num_edges: CSC_dst[search];
      return (end-begin);
    }

    /** \brief Returns the number of out edges of the vertex. */
    size_t num_out_edges (const lvid_type v) const {
      if (v >= num_vertices)
        return 0;

      ASSERT_LT(v, num_vertices);
      size_t begin = CSR_src[v];
      if (begin >= num_edges) return 0;

      size_t search = use_skip_list ? nextValid(CSR_src_skip, v, true) : nextValid(CSR_src, v, false);
      size_t end = (search >= num_vertices) ? num_edges: CSR_src[search];
      return (end-begin);
    }

    /** \brief Returns the edge id of the edge. 
     * Edges are assigned with consecutive int ids,
     * ordered first by source and then by target.
     * */
    edge_id_type edge_id(const edge_type& edge) const {
      ASSERT_FALSE(edge.empty());
      return edge.get_dir() == OUT_EDGES ?
        edge._edge_id :
        c2r_map[edge._edge_id];
    }

    /** \brief Returns the reference of edge data of an edge. */
    edge_data_type& edge_data(lvid_type source, lvid_type target) {
      ASSERT_LT(source, num_vertices);
      ASSERT_LT(target, num_vertices);
      edge_type ans = find(source, target);
      return edge_data(ans);
    }

    /** \brief Returns the constant reference of edge data of an edge. */
    const edge_data_type& edge_data(lvid_type source, 
                                    lvid_type target) const {
      ASSERT_LT(source, num_vertices);
      ASSERT_LT(target, num_vertices);
      edge_type ans = find(source, target);
      return edge_data(ans);
    }

    /** \brief Returns the reference of edge data of an edge. */
    edge_data_type& edge_data(edge_type edge) {
      ASSERT_FALSE(edge.empty());
      return edge_data_list[edge.get_dir() == OUT_EDGES ?
                            edge._edge_id :
                            c2r_map[edge._edge_id]];
    }

    /** \brief Returns the constant reference of edge data of an edge. */
    const edge_data_type& edge_data(edge_type edge) const {
      ASSERT_FALSE(edge.empty());
      return edge_data_list[edge.get_dir() == OUT_EDGES ?
                            edge._edge_id :
                            c2r_map[edge._edge_id]];
    }

    /** \brief Returns a list of in edges of a vertex. */
    edge_list in_edges(const lvid_type v) const {
      if (v >= num_vertices)
        return edge_list();

      std::pair<bool, edge_range_type> rangePair = inEdgeRange(v);
      if (rangePair.first) {
        edge_range_type range = rangePair.second;
        edge_dir_type dir = IN_EDGES;

        edge_iterator begin (v, range.first, dir, &(CSC_src[0]));
        edge_iterator end (v, range.second+1, dir, &(CSC_src[0]));
        // std::cout << "in range (" << range.first << "," <<
        // range.second << ")" << std::endl; std::cout << "in edge
        // size: " << end-begin << std::endl;
        return edge_list(begin, end);
      } else { return edge_list(); }
    }

    /** \brief Returns a list of out edges of a vertex. */
    edge_list out_edges(const lvid_type v) const {
      if (v >= num_vertices)
        return edge_list();

      std::pair<bool, edge_range_type> rangePair = outEdgeRange(v);
      if (rangePair.first) {
        edge_range_type range = rangePair.second;
        edge_iterator begin (v, range.first, OUT_EDGES, &(CSR_dst[0]));
        edge_iterator end (v, range.second+1, OUT_EDGES, &(CSR_dst[0]));
        // std::cout << "out range (" << range.first << "," <<
        // range.second << ")" << std::endl; std::cout << "out_edge
        // size: " << end-begin << std::endl;
        return edge_list(begin, end);
      } else { return edge_list(); }
    }

    /** \brief Returns an edge type of a given source target pair. */
    edge_type find (const lvid_type src, 
                    const lvid_type dst) const {
      // DEBUG printf("Find: %u, %u \n", src, dst);
      
      // Get the out edge range of the src, as well as the in edge
      // range of the dst.
      // Search CSR or CSC, whichever has less candidates.
      std::pair<bool, edge_range_type> dstRangePair = outEdgeRange(src);
      std::pair<bool, edge_range_type> srcRangePair = inEdgeRange(dst);
      if( srcRangePair.first && dstRangePair.first) {
        // The edge may exist. 
        edge_range_type srcRange =  srcRangePair.second;
        edge_range_type dstRange = dstRangePair.second;

        if ((srcRange.second - srcRange.first) < 
            (dstRange.second - dstRange.first)) {
          // Out edge candidate size is smaller, search CSC.
          size_t efind =  binary_search(CSC_src, srcRange.first, 
                                        srcRange.second, src);
          return efind >= num_edges ? 
            edge_type() : edge_type(dst, src, efind, IN_EDGES);
        } else {
          // In edge candidate size is smaller, search CSR.
          size_t efind = binary_search(CSR_dst, dstRange.first, 
                                       dstRange.second, dst);
          return efind >= num_edges ? 
            edge_type() : edge_type(src, dst, efind, OUT_EDGES);
        }
      } else {
        return edge_type();
      }
    } // end of find.

     /** \brief Finalize the graph storage. 
      * Construct the CSC, CSR, by sorting edges to maximize the
      * efficiency of graphlab.  
      * This function takes O(|V|log(degree)) time and will 
      * fail if there are any duplicate edges.
      *
      * Assumption: 
      * _num_of_v == 1 + max(max_element(edges.source_arr), max_element(edges.target_arr))
      */
    void finalize(size_t _num_of_v, edge_info &edges) {
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize starts." << std::endl;
#endif
      num_vertices = _num_of_v;
      num_edges = edges.size();

      // Permute_index, alias of c2r_map. Confusing but efficient.
      std::vector<edge_id_type>& permute_index = c2r_map;
      permute_index.reserve(num_edges); 

      // Counter_index.
      std::vector<atomic<int> > counter_array(num_vertices+1, 0);
      permute_index.assign(num_edges, 0);

      // Sort edges by source;
      // Begin of counting sort.
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize: Sort by source vertex" << std::endl;
#endif
      counting_sort(edges.source_arr, counter_array, permute_index); 

      // Parallel sort target for each source= x interval: counter_array[x] - counter_array[x+1];
#ifdef _OPENMP
#pragma omp parallel for
#else
      logstream(LOG_DEBUG) << "Graph2 finalize: Parallel sort is disabled." << std::endl;
#endif

      for (ssize_t j = 0; j < ssize_t(num_vertices); ++j) {
        if (counter_array[j] < counter_array[j+1]) {
          std::sort(permute_index.begin()+counter_array[j], 
                    permute_index.begin()+counter_array[j+1],
                    cmp_by_any_functor<lvid_type> (edges.target_arr)); 
        }
      }
      // End of counting sort.

      // Inplace permute of edge_data, edge_src, edge_target array.
      // Modified from src/graphlab/util/generics/shuffle.hpp.
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize: Inplace permute by source vertex" << std::endl;
#endif
      lvid_type swap_src; lvid_type swap_target;
      for (size_t i = 0; i < permute_index.size(); ++i) {
        if (i != permute_index[i]) {
          // Reserve the ith entry;
          size_t j = i;
          EdgeData swap_data = edges.data[i];
          swap_src = edges.source_arr[i];
          swap_target = edges.target_arr[i];
          // Begin swap cycle:
          while (j != permute_index[j]) {
            size_t next = permute_index[j];
            if (next != i) {
              edges.data[j] = edges.data[next];
              edges.source_arr[j] = edges.source_arr[next];
              edges.target_arr[j] = edges.target_arr[next];
              permute_index[j] = j;
              j = next;
            } else {
              // end of cycle
              edges.data[j] = swap_data;
              edges.source_arr[j] = swap_src;
              edges.target_arr[j] = swap_target;
              permute_index[j] = j;
              break;
            }
          }
        }
      }

      bool duplicate_edge_warn = false;

      // Construct CSR_src:
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG)<< "Graph2 finalize: build CSR_src..." << std::endl;
#endif
      CSR_src.reserve(num_vertices);
      if (use_skip_list) {
        CSR_src_skip.reserve(num_vertices);
      }
      size_t lastSrc = -1;
      lvid_type old_src = -1;
      lvid_type old_dst = -1;
      // Iterate over the edges. 
      for (size_t it = 0; it < num_edges; ++it) {
        lvid_type src = edges.source_arr[it];
        lvid_type dst = edges.target_arr[it];
        // Check duplicate edge.
        if (src == old_src && dst == old_dst) {
          if (!duplicate_edge_warn)
            logstream(LOG_WARNING)
              << "Duplicate edge "
              << it << ":(" << src << ", " << dst << ") "
              << "found! Graphlab does not support graphs "
              << "with duplicate edges. This error will be reported only once." << std::endl;
          duplicate_edge_warn = true;
          continue;
        } else {
          old_src = src;
          old_dst = dst;
        }
        // Fill in CSR_src and CSR_src_skip. 
        if (src != lastSrc) {
          for (size_t j = (lastSrc+1); j < src; ++j) {
            CSR_src.push_back(-1);
            if (use_skip_list) 
              CSR_src_skip.push_back(src-lastSrc-1);
          }
          CSR_src.push_back(it);
          if (use_skip_list) 
            CSR_src_skip.push_back(0);
          lastSrc = src;
        }
      }
      // Fill in the remaining row index list.
      for( size_t j = (lastSrc +1); j < num_vertices; ++j) {
        CSR_src.push_back(-1);
        if (use_skip_list) 
          CSR_src_skip.push_back(num_vertices-lastSrc-1);
      }
      ASSERT_EQ(CSR_src.size(), num_vertices);
      if (use_skip_list)
        ASSERT_EQ(CSR_src_skip.size(), num_vertices);

      // End of building CSR


      // Begin building CSC
      // Directed graph need both CSC and CSR
      // Construct c2r_map, sort the ids according to column first order.
      // Begin of counting sort.
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize: Sort by source vertex" << std::endl;
#endif
      counting_sort(edges.target_arr, counter_array, permute_index); 
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ssize_t i = 0; i < ssize_t(num_vertices); ++i) {
        if (counter_array[i] < counter_array[i+1]) {
          std::sort(permute_index.begin()+counter_array[i],
                    permute_index.begin() + counter_array[i+1],
                    cmp_by_any_functor<lvid_type>(edges.source_arr)); 
        }
      }
      // End of counting sort.

#ifdef AVOID_OUTOFPLACE_PERMUTE
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize: Inplace permute by target vertex" << std::endl;
#endif
      inplace_shuffle(edges.source_arr.begin(), edges.source_arr.end(), permute_index);
      /// YUCHENG to JAY: why is this here? It looks like a copy and paste from above
//       counting_sort(edges.target_arr, counter_array, permute_index); 
//       for (ssize_t i = 0; i < ssize_t(num_vertices); ++i) {
//         if (counter_array[i] < counter_array[i+1]) {
//           std::sort(permute_index.begin()+counter_array[i],
//                     permute_index.begin() + counter_array[i+1],
//                     cmp_by_any_functor<lvid_type>(edges.source_arr)); 
//         }
//       }
#else
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) << "Graph2 finalize: Outofplace permute by target vertex" << std::endl;
#endif
      outofplace_shuffle(edges.source_arr, permute_index);
#endif

      /* DEBUG
         printf("c2r_map: \n");
         foreach(edge_id_type e, c2r_map)
         std::cout << e << " ";
         std::cout << std::endl;
      */


      // Construct CSC_dst:
      CSC_dst.reserve(num_vertices);
      if (use_skip_list) {
        CSC_dst_skip.reserve(num_vertices);
      }
      size_t lastDst = -1;
#ifdef DEBUG_GRAPH
      logstream(LOG_DEBUG) <<"Graph2 finalize: Build CSC_dst..." << std::endl;
#endif
      // Iterate over the edges. 
      for (size_t it = 0; it < num_edges; ++it) {
        lvid_type dst = edges.target_arr[c2r_map[it]];

        // Fill in CSC_dst and CSR_src_skip. 
        if (dst != lastDst) {
          for (size_t j = (lastDst + 1); j < dst; ++j) {
            CSC_dst.push_back(-1);
            if (use_skip_list) 
              CSC_dst_skip.push_back(dst-lastDst-1);
          }
          CSC_dst.push_back(it);
          if (use_skip_list) 
            CSC_dst_skip.push_back(0);
          lastDst = dst;
        }
      }
      // Fill in the remaining row index list.
      for( size_t j = (lastDst +1); j < num_vertices; ++j) {
        CSC_dst.push_back(-1);
        if (use_skip_list) 
          CSC_dst_skip.push_back(num_vertices-lastDst-1);
      }
      ASSERT_EQ(CSC_dst.size(), num_vertices);
      if (use_skip_list)
        ASSERT_EQ(CSC_dst_skip.size(), num_vertices);

      // Swap edges.source with CSC_src
      CSC_src.swap(edges.source_arr);
      // End of building CSC

      
      // Swap edges.target with CSR_dst
      CSR_dst.swap(edges.target_arr);
      // Swap edge data and perserve c2r_map.
      edge_data_list.swap(edges.data);
#ifdef DEBGU_GRAPH
      logstream(LOG_DEBUG) << "End of finalize." << std::endl;
#endif

      /* DEBUG */
      // printf("CSR dst:\n");
      // foreach(lvid_type i, CSR_dst)
      //   std::cout << i << " ";
      // std::cout << std::endl;
      // printf("CSR src:\n");
      // foreach(size_t i, CSR_src)
      //   std::cout << i << " "; 
      // std::cout << std::endl;
      // printf("CSR src skip:\n");
      // foreach(size_t i, CSR_src_skip)
      //   std::cout << i << " "; 
      // std::cout << std::endl;

      /* DEBUG */
      // printf("CSC dst:\n");
      // foreach(lvid_type i, CSC_dst)
      //   std::cout << i << " ";
      // std::cout << std::endl;
      // printf("CSC src:\n");
      // foreach(size_t i, CSC_src)
      //   std::cout << i << " "; 
      // std::cout << std::endl;
      // printf("CSC dst skip:\n");
      // foreach(size_t i, CSC_dst_skip)
      //   std::cout << i << " "; 
      // std::cout << std::endl;
     
    } // end of finalize.

    /** \brief Reset the storage. */
    void clear() {
      CSR_src.clear();
      CSR_dst.clear();
      CSC_src.clear();
      CSC_dst.clear();
      CSR_src_skip.clear();
      CSC_dst_skip.clear();
      c2r_map.clear();
      edge_data_list.clear();
    }

    /** \brief Reset the storage and free the reserved memory. */
    void clear_reserve() {
      std::vector<edge_id_type>().swap(CSR_src);
      std::vector<lvid_type>().swap(CSR_dst);
      std::vector<lvid_type>().swap(CSC_src);
      std::vector<edge_id_type>().swap(CSC_dst);
      std::vector<edge_id_type>().swap(c2r_map);
      std::vector<EdgeData>().swap(edge_data_list);
      std::vector<edge_id_type>().swap(CSR_src_skip);
      std::vector<edge_id_type>().swap(CSC_dst_skip);
    }

    size_t estimate_sizeof() const {
      // const size_t word_size = sizeof(size_t);
      const size_t vid_size = sizeof(lvid_type);
      const size_t eid_size = sizeof(edge_id_type);
      // Actual content size;
      const size_t CSR_size = eid_size * CSR_src.capacity() + 
        vid_size * CSR_dst.capacity();
      const size_t CSC_size = eid_size *CSC_dst.capacity() + 
        vid_size * CSC_src.capacity() + eid_size * c2r_map.capacity();
      const size_t edata_size = sizeof(EdgeData) * edge_data_list.capacity();

      // Container size;
      const size_t container_size = sizeof(CSR_src) + sizeof(CSR_dst) + 
        sizeof(CSC_src) + sizeof(CSC_dst) + sizeof(c2r_map) + 
        sizeof(edge_data_list);
      // Skip list size:
      const size_t skip_list_size = sizeof(CSR_src_skip) + 
        sizeof(CSC_dst_skip) + CSR_src_skip.capacity() * vid_size + 
        CSC_dst_skip.capacity() * vid_size;

      logstream(LOG_DEBUG) << "CSR size: " 
                << (double)CSR_size/(1024*1024)
                << " CSC size: " 
                << (double)CSC_size/(1024*1024) 
                << " edata size: "
                << (double)edata_size/(1024*1024)
                << " skiplist size: " 
                << (double)(skip_list_size)/(1024*1024)
                << " container size: " 
                << (double)container_size/(1024*1024) 
                << " \n Total size: " 
                << double(CSR_size + CSC_size + container_size + skip_list_size) << std::endl;

      return CSR_size + CSC_size + edata_size + container_size + 
        skip_list_size;
    } // end of estimate_sizeof

    /** To be deprecated. */
    // This is a log(V) operation.
    // Do not use operations on edge id.
    // Use edge_list instead.
    lvid_type source(edge_id_type eid) const {
      ASSERT_LT(eid, num_edges);
      // Helper function: binary search the CSR_row;
      return lookup_source(eid);
    }

    // This is a log(V) operation.
    // Do not use operations on edge id.
    // Use edge_list instead.
    lvid_type target(edge_id_type eid) const {
      ASSERT_LT(eid, num_edges);
      return CSR_dst[eid];
    }



    // ------------- Private data storage ----------------
  private:
    /** Number of vertices in the storage (not counting singletons)*/
    size_t num_vertices;
    /** Number of edges in the storage. */
    size_t num_edges;

    /** Array of edge data sorted by source vid. */
    std::vector<EdgeData> edge_data_list;

    /** \internal 
     * Row index of CSR, corresponding to the source vertices. */
    std::vector<edge_id_type> CSR_src;

    /** 
     * \internal 
     * Suppose CSR_src is: 1 x x 3 x x x x 5
     * where x means no out edges.
     *     CSR_src_skip =  0 2 2 0 4 0 0 4 0
     * is used to jump to the prev/next valid vertex in CSR_src.  
     * Optional.
     */
    std::vector<edge_id_type> CSR_src_skip;

    /** \internal 
     * Col index of CSR, corresponding to the target vertices. */
    std::vector<lvid_type> CSR_dst;

    /** \internal 
     * Map the sort-by-col edge id to sort-by-row edge id */
    std::vector<edge_id_type> c2r_map;

    /** \internal
     * Row index of CSC, corresponding to the target vertices. */
    std::vector<edge_id_type> CSC_dst;

    /* 
     * \internal 
     * Suppose CSC_dst is: 1 x x 3 x x x x 5
     * where x means no in edges.
     *     CSC_dst_skip =  0 2 2 0 4 0 0 4 0
     * is used to jump to the prev/next valid vertex in CSC_dst.  
     * Optional.
     */
    std::vector<edge_id_type> CSC_dst_skip;
    /** \internal
     * Col index of CSC, corresponding to the source vertices. */
    std::vector<lvid_type> CSC_src;

    /** Graph storage traits. */
    bool use_skip_list;


 /****************************************************************************
 *                       Internal Functions                                 *
 *                     ----------------------                               *
 * These functions functions and types provide internal access to the       *
 * underlying graph representation. They should not be used unless you      *
 * *really* know what you are doing.                                        *
 ****************************************************************************/
  private:
    /** \internal
     *  Returns the begin and end index of the in edge of vertex v. */
    inline std::pair<bool, edge_range_type> inEdgeRange(lvid_type v) const {
      ASSERT_LT(v, num_vertices);

      size_t col_start = CSC_dst[v];
      if (col_start >= num_edges) {
        // No inbound edges.
        return std::make_pair(false, std::make_pair(0,0));
      } else {

        // Find the start column of the next vertex.
        lvid_type nextV = use_skip_list ? nextValid(CSC_dst_skip, v, true) : 
          nextValid(CSC_dst, v, false);
        size_t col_end = (nextV < num_vertices) ? CSC_dst[nextV] : num_edges;
        return std::make_pair(true, std::make_pair(col_start, col_end-1));
      }
    } // End of inEdgeRange;

    /** \internal
     *  Returns the begin and end index of the out edge of vertex v. */
    inline std::pair<bool, edge_range_type> outEdgeRange(lvid_type v) const {
      ASSERT_LT(v, num_vertices);
      size_t row_start = CSR_src[v];
      if (row_start >= num_edges) {
        // No outbound edges.
        return std::make_pair(false, std::make_pair(0,0));;
      } else {
        // Find the start column of the next vertex.
        lvid_type nextV = use_skip_list ? nextValid(CSR_src_skip, v, true) : 
          nextValid(CSR_src, v, false);
        size_t row_end = (nextV < num_vertices) ? CSR_src[nextV] :num_edges; 

        return std::make_pair(true, std::make_pair(row_start, row_end-1));
      }
    } // End of outEdgeRange;


    //-------------Private Helper functions------------
    /** \internal
     *  Compare functor of any type*/
    template <typename anyvalue>
    struct cmp_by_any_functor {
      const std::vector<anyvalue>& vec;
      cmp_by_any_functor(const std::vector<anyvalue>& _vec) : vec(_vec) { }
      bool operator()(size_t me, size_t other) const {
        return (vec[me] < vec[other]);
      }
    };

    /** \internal
     *  Counting sort vector in ascending order and 
     *  fill the counting array and permute index array. */
    template <typename valuetype>
    void counting_sort(const std::vector<valuetype>& value_array, std::vector< atomic<int> >& counter_array, std::vector<edge_id_type>& permute_index) {
      counter_array.assign(counter_array.size(), 0);
      permute_index.assign(permute_index.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ssize_t i = 0; i < ssize_t(value_array.size()); ++i) {
        size_t val = value_array[i];
        counter_array[val].inc();
      }

      for (size_t i = 1; i < counter_array.size(); ++i) {
        counter_array[i] += counter_array[i-1];
      }
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (ssize_t i = 0; i < ssize_t(value_array.size()); ++i) {
        size_t val = value_array[i];
        permute_index[counter_array[val].dec()] = i;
      }
    }

    /** \internal
     *  Binary search vfind in a vector of lvid_type 
     *  within range [start, end]. Returns (size_t)(-1) if not found. */
    size_t binary_search(const std::vector<lvid_type>& vec, 
                         size_t start, size_t end, 
                         lvid_type vfind) const {
      ASSERT_LT(vfind, num_vertices);
      while(start <= end) {
        size_t mid = (start+end)/2;
        lvid_type vpoke = vec[mid];
        if(vpoke == vfind) {
          return mid;
        } else if (vpoke > vfind) {
          end = mid - 1;
        } else {
          start = mid + 1;
        }
      }
      // Not found;
      return -1;
    }// End of binary_search

    /** \internal
     * To be deprecated. */
    // This is a log(V) operation.
    // Do not use operations on edge id.
    // Use edge_list instead.
    lvid_type
    lookup_source(edge_id_type eid) const {
      // Binary search: find i such that CSR_src[i] <= eid, CSR_src[i+1] > eid;
      ASSERT_LT(eid, num_edges);
      size_t start = 0;
      size_t end = num_vertices-1;
      if (CSR_src[start] >= num_edges) start = use_skip_list ? nextValid(CSR_src_skip, start, true) : nextValid(CSR_src, start, false);
      if (CSR_src[end] >= num_edges) end = use_skip_list ? prevValid(CSR_src_skip, end, true) : prevValid(CSR_src, end, false);

      // Keep the invariant that CSR_src[start] <= eid, CSR_src[end] > eid.
      while (start <= end) {
        if (CSR_src[end] <= eid)
          return end;
        if (CSR_src[start] == eid)
          return start;
        if (CSR_src[start] > eid)
          {
            ASSERT_LT(0, start);
            lvid_type ans = use_skip_list ? prevValid(CSR_src_skip, start, true) : prevValid(CSR_src, start, false);
            //lvid_type ans = start -1;
            ASSERT_LT(ans, num_vertices);
            return ans;
          }

        size_t mid = (start + end)/2;
        // mid may fall in to an invalid grid 
        if (CSR_src[mid] >= num_edges) mid = prevValid(CSR_src, mid, false);

        if (CSR_src[mid] == eid)
          return mid;

        if (CSR_src[mid] > eid) {
          end = use_skip_list ? prevValid(CSR_src_skip, mid, true) : prevValid(CSR_src, mid, false); 
          //end = mid-1;
        } else {
          //start = mid+1;
          start = use_skip_list ? nextValid(CSR_src_skip, mid, true) : nextValid(CSR_src, mid, false);
        }
      } 
            
      ASSERT_TRUE(false);
      return start;
    }
    
    /** \internal
     * Returns the next valid vertex to the current vertex id in vertex_array.
     * Return num_vertices if there is no valid vertex next to the curent vertex.
     *
     * A vertex is valid if it has out edges (or in edges).
     *
     * If use_skip_list = false 
     * vertex_array can be: CSR_src, or CSC_dst
     *
     * If use_skip_list = true
     * vertex_array can be: CSR_src_skip, or CSC_dst_skip
     * Uses skip array to find the next valid vertex in O(1).
     *
     * Can be applied on invalid vertex. 
     * This function is useful in binary search where the middle is not
     * assumed to be valid. */
    inline lvid_type 
    nextValid(const std::vector<edge_id_type>& vertex_array, 
              lvid_type curv, bool use_skip_list) const {

      if (curv == num_vertices-1) return num_vertices;

      if (use_skip_list) {
        return (curv + 1 + vertex_array[curv+1]);
      } else {
        lvid_type search = curv+1;
        while (search < num_vertices && vertex_array[search] >= num_edges) 
          ++search;
        return search;
      }
    }

    /** \internal
     * Symmetric function of nextValid.
     * Returns the previous valid vertex to the current vertex id in vertex_array.
     * Return num_vertices if there is no valid vertex previous to the curent vertex.
     */
    inline lvid_type 
    prevValid(const std::vector<lvid_type>& vertex_array, 
              lvid_type curv, bool use_skip_list) const {
      if (curv == 0) return -1;

      if (use_skip_list) {
        return (curv - 1 - vertex_array[curv-1]);
      } else {
        lvid_type search = curv-1;
        while (search >= 0 && vertex_array[search] >= num_edges)
          --search;
        return search;
      }
    }

  public:

    /** \internal
     * Returns a reference of CSR_src.*/
    const std::vector<lvid_type>& get_csr_src() const {
      return CSR_src;
    }
    /** \internal
     * Returns a reference of CSR_dst.*/
    const std::vector<edge_id_type>& get_csr_dst() const {
      return CSR_dst;
    }
    /** \internal
     * Returns a reference of CSC_src.*/
    const std::vector<edge_id_type>& get_csc_src() const {
      return CSC_src;
    }
    /** \internal
     * Returns a reference of CSC_dst.*/
    const std::vector<lvid_type>& get_csc_dst() const {
      return CSC_dst;
    }
    /** \internal
     * Returns a reference of edge_data_list.*/
    const std::vector<EdgeData>& get_edge_data() const {
      return edge_data_list;
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();
      arc >> use_skip_list
          >> num_vertices
          >> num_edges
          >> edge_data_list
          >> CSR_src
          >> CSR_dst
          >> CSC_src
          >> CSC_dst
          >> c2r_map
          >> CSR_src_skip
          >> CSC_dst_skip;
    }

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      arc << use_skip_list
          << num_vertices
          << num_edges
          << edge_data_list
          << CSR_src
          << CSR_dst
          << CSC_src
          << CSC_dst
          << c2r_map
          << CSR_src_skip
          << CSC_dst_skip;
    }

    /** swap two graph storage*/
    void swap(graph_storage& other) {
      std::swap(use_skip_list, other.use_skip_list);
      std::swap(num_vertices, other.num_vertices);
      std::swap(num_edges, other.num_edges);
      std::swap(edge_data_list, other.edge_data_list);
      std::swap(CSR_src, other.CSR_src);
      std::swap(CSR_dst, other.CSR_dst);
      std::swap(CSC_src, other.CSC_src);
      std::swap(CSC_dst, other.CSC_dst);
      std::swap(c2r_map, other.c2r_map);
      std::swap(CSR_src_skip, other.CSR_src_skip);
      std::swap(CSC_dst_skip, other.CSC_dst_skip);
    }

  };// End of graph store;
}// End of namespace;

namespace std {
  template<typename VertexData, typename EdgeData>
  inline void swap(graphlab::graph_storage<VertexData,EdgeData>& a, 
                   graphlab::graph_storage<VertexData,EdgeData>& b) {
    a.swap(b);
  } // end of swap
}; // end of std namespace


#include <graphlab/macros_undef.hpp>
#endif
