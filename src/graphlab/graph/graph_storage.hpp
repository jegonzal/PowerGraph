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
#include <omp.h>
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
#include <graphlab/macros_def.hpp>

namespace graphlab {

  template<typename VertexData, typename EdgeData>
  class graph_storage {
  public:
    typedef graphlab::vertex_id_type vertex_id_type;
    typedef graphlab::edge_id_type edge_id_type;

    typedef EdgeData edge_data_type;
    typedef VertexData vertex_data_type;
    typedef std::pair<size_t, size_t>  edge_range_type;

    /* ----------------------------------------------------------------------------- */
    /* Helper data field and structures: edge_data_list, class edge, class edge_list */
    /* ----------------------------------------------------------------------------- */
  public:
    // Edge class for temporary storage. Will be finalized into the CSR+CSC form.
    class edge_info {
      public:
        std::vector<EdgeData> data;
        std::vector<vertex_id_type> source_arr;
        std::vector<vertex_id_type> target_arr;
      public:
        edge_info () {}
        void add_edge(vertex_id_type source, vertex_id_type target, EdgeData _data) {
          data.push_back(_data);
          source_arr.push_back(source);
          target_arr.push_back(target);
        }
        void add_block_edges(const std::vector<vertex_id_type>& src_arr, 
            const std::vector<vertex_id_type>& dst_arr, 
            const std::vector<EdgeData>& edata_arr) {
          data.insert(data.end(), edata_arr.begin(), edata_arr.end());
          source_arr.insert(source_arr.end(), src_arr.begin(), src_arr.end());
          target_arr.insert(target_arr.end(), dst_arr.begin(), dst_arr.end());
        }
        void clear() {
          std::vector<EdgeData>().swap(data);
          std::vector<vertex_id_type>().swap(source_arr);
          std::vector<vertex_id_type>().swap(target_arr);
        }
        size_t size() const {
          return source_arr.size();
        }
        size_t estimate_sizeof() const {
          return data.capacity()*sizeof(EdgeData) + source_arr.capacity()*sizeof(vertex_id_type)*2 + sizeof(data) + sizeof(source_arr)*2 + sizeof(edge_info);
        }
    }; // end of class edge_info.

    // A class of edge information. Used as value type of the edge_list.
    class edge_type {
      public:
        enum edge_dir{OUTEDGE, INEDGE, NONE};
      public:
        edge_type () : _source(-1), _target(-1), _edge_id(-1), _dir(NONE), _empty(true) { }
        edge_type (const vertex_id_type _source, const vertex_id_type _target, 
            const edge_id_type _eid, edge_dir _dir) :
          _source(_source), _target(_target), _edge_id(_eid), _dir(_dir), _empty(false) { }
      public:
        inline vertex_id_type source() const {
          ASSERT_FALSE(empty()); 
          return _dir==OUTEDGE ? _source : _target; 
        }

        inline vertex_id_type target() const { 
          ASSERT_FALSE(empty());
          return _dir==OUTEDGE ? _target : _source; 
        }

        inline edge_dir get_dir() const {
          return _dir;
        }

        // Do not use.
        inline uint offset() const {
          return _edge_id;
        }

        // Do not use.
        edge_id_type edge_id() const {
          return _edge_id;
        }

        inline bool empty() const { return _empty; }
        // Data fields. 
      private:
        vertex_id_type _source;
        vertex_id_type _target;
        edge_id_type _edge_id;
        edge_dir _dir;
        bool _empty;

        friend class graph_storage;
    }; // end of class edge_type.

    // Internal iterator on edge_types.
    class edge_iterator : 
      public std::iterator<std::forward_iterator_tag, edge_type> {
    public:
      typedef typename edge_type::edge_dir iterator_type;
      typedef edge_type reference;
    public:
      // Cosntructors
      edge_iterator () : offset(-1), empty(true) { }
     
      edge_iterator (vertex_id_type _center, size_t _offset, 
                     iterator_type _itype, const edge_id_type* _vid_arr) : 
        center(_center), offset(_offset), itype(_itype), vid_arr(_vid_arr), empty(false) { }
      
      edge_iterator (const edge_iterator& it) :
        center(it.center), offset(it.offset), itype(it.itype), vid_arr(it.vid_arr),
       empty(it.empty) { }
  
      inline edge_type operator*() const  {
        ASSERT_TRUE(!empty);
        return make_value();
      }

      typedef boost::detail::
      operator_arrow_result<edge_type, edge_type, edge_type*> arrow_type;
      inline typename arrow_type::type operator->() const {
        return arrow_type::make(make_value());
      }


      inline bool operator==(const edge_iterator& it) const {
        return (empty && it.empty) || 
          (empty == it.empty && itype == it.itype && center == it.center && 
           offset == it.offset);
      }

      inline bool operator!=(const edge_iterator& it) const { 
        return !(*this == it);
      }

      inline edge_iterator& operator++() {
        ASSERT_TRUE(!empty);
        ++offset;
        return *this;
      }

      inline edge_iterator operator++(int) {
        ASSERT_TRUE(!empty);
        const edge_iterator copy(*this);
        operator++();
        return copy;
      }


      inline int operator-(const edge_iterator& it) const {
        ASSERT_TRUE(!empty && itype == it.itype && center == it.center);
        return offset - it.offset;
      }

      inline edge_iterator operator+(size_t i) const {
        edge_iterator retval(center, offset+i, itype, vid_arr);
        return retval;
      }

      // Generate the ret value of the iterator.
      inline edge_type make_value() const {
          return edge_type(center, vid_arr[offset], offset, itype);
      }

    private:
      vertex_id_type center;
      size_t offset;
      iterator_type itype;
      const edge_id_type* vid_arr;
      bool empty;
    }; // end of class edge_iterator.

    // Represents an iteratable list of edge_types.
    class edge_list {
      // Type interface for boost foreach.
    public:
      typedef edge_iterator iterator;
      typedef edge_iterator const_iterator;
      typedef edge_type value_type;

    public:
      // Construct an empty edge list
      edge_list() : list_size(0) { }
      // Cosntruct an edge_list with begin and end. 
      edge_list(edge_iterator begin, edge_iterator end) : 
        begin_ptr(begin), end_ptr(end) { 
        list_size = (size_t)(end_ptr-begin_ptr);
      }
      // Copy constructor
      edge_list(const edge_list& other) : 
        begin_ptr(other.begin_ptr), end_ptr(other.end_ptr), 
        list_size(other.list_size) { }

      inline size_t size() const { return list_size;}
            
      inline edge_type operator[](size_t i) const {
        ASSERT_LT(i, list_size);
        return *(begin_ptr + i);
      }

      iterator begin() const { return begin_ptr; }
      iterator end() const { return end_ptr; }
      bool empty() const { return size() == 0; }

    private:
      edge_iterator begin_ptr;
      edge_iterator end_ptr;
      size_t list_size;
    }; // end of class edge_list.

  public:
    graph_storage() : use_skip_list(false) { 
    }

    void set_is_directed (bool x) {
      logstream(LOG_INFO) << "Using graph storage that only support directed graph. Undirect setting will be ignored." << std::endl;
    }
    void set_use_skip_list (bool x) { use_skip_list = x;}
    void get_is_directed () const { return true; }

    size_t edge_size() const { return num_edges; }

    size_t vertices_size() const { return num_vertices; }

    // Return the size of in_neighbours.
    size_t num_in_edges (const vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);

      size_t begin = CSC_dst[v];
      if (begin >= num_edges) return 0;
      // Search is the next valid src vertex after v.
      size_t search = use_skip_list ? nextValid(CSC_dst_skip, v, true) : nextValid(CSC_dst, v, false);
      size_t end = (search >= num_vertices) ? num_edges: CSC_dst[search];
      return (end-begin);
    }

    // Return the size of out_neighbours.
    size_t num_out_edges (const vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);
      size_t begin = CSR_src[v];
      if (begin >= num_edges) return 0;

      size_t search = use_skip_list ? nextValid(CSR_src_skip, v, true) : nextValid(CSR_src, v, false);
      size_t end = (search >= num_vertices) ? num_edges: CSR_src[search];
      return (end-begin);
    }


    edge_id_type edge_id (const edge_type& e) const {
      ASSERT_FALSE(e.empty());
      if (e.get_dir() == edge_type::OUT_EDGE)
      {
        return e.edge_id();
      } else {
        return e.edge_id() + num_edges;
      }
    }


    edge_data_type& edge_data(vertex_id_type source, vertex_id_type target) {
      ASSERT_LT(source, num_vertices);
      ASSERT_LT(target, num_vertices);
      edge_type ans = find(source, target);
      return edge_data(ans);
    }

    const edge_data_type& edge_data(vertex_id_type source, 
                                    vertex_id_type target) const {
      ASSERT_LT(source, num_vertices);
      ASSERT_LT(target, num_vertices);
      edge_type ans = find(source, target);
      return edge_data(ans);
    }

    edge_data_type& edge_data(edge_type edge) {
      ASSERT_FALSE(edge.empty());
      return edge_data_list[edge.get_dir() == edge_type::OUTEDGE ? 
        edge._edge_id :
        c2r_map[edge._edge_id]];
    }

    const edge_data_type& edge_data(edge_type edge) const {
      ASSERT_FALSE(edge.empty());
      return edge_data_list[edge.get_dir() == edge_type::OUTEDGE ? 
        edge._edge_id :
        c2r_map[edge._edge_id]];
    }


    // Return in edge list of a vertex.
    edge_list in_edges(const vertex_id_type v) const {
      std::pair<bool, edge_range_type> rangePair = inEdgeRange(v);
      if (rangePair.first) {
        edge_range_type range = rangePair.second;

        typedef typename edge_type::edge_dir edge_dir;

        edge_dir dir = edge_type::INEDGE;

        edge_iterator begin (v, range.first, dir, &(CSC_src[0]));

        edge_iterator end (v, range.second+1, dir, &(CSC_src[0]));
        // std::cout << "in range (" << range.first << "," <<
        // range.second << ")" << std::endl; std::cout << "in edge
        // size: " << end-begin << std::endl;
        return edge_list(begin, end);
      } else { return edge_list(); }
    }

    // Return out edge list of a vertex.
    edge_list out_edges(const vertex_id_type v) const {
      std::pair<bool, edge_range_type> rangePair = outEdgeRange(v);
      if (rangePair.first) {
        edge_range_type range = rangePair.second;
        edge_iterator begin (v, range.first, edge_type::OUTEDGE, &(CSR_dst[0]));
        edge_iterator end (v, range.second+1, edge_type::OUTEDGE, &(CSR_dst[0]));
        // std::cout << "out range (" << range.first << "," <<
        // range.second << ")" << std::endl; std::cout << "out_edge
        // size: " << end-begin << std::endl;
        return edge_list(begin, end);
      } else { return edge_list(); }
    }

    // Return the edge id given source and destination.
    edge_type find (const vertex_id_type src, 
                    const vertex_id_type dst) const {
      /* DEBUG
         printf("Find: %u, %u \n", src, dst);
      */
      // Get the out edge range of the src, as well as the in edge
      // range of the dst.
        // Directed graph, search CSR or CSC, whichever has less candidates.
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
            if (efind >= num_edges) {
              return edge_type();
            } else {
                return edge_type(dst, src, efind, edge_type::INEDGE);
            }
          } else {
            // In edge candidate size is smaller, search CSR.
            size_t efind = binary_search(CSR_dst, dstRange.first, 
                                         dstRange.second, dst);
            if (efind >= num_edges) {
              return edge_type();
            } else {
              return edge_type(src, dst, efind, edge_type::OUTEDGE);
            }
          }
        } else {
          return edge_type();
        }
    } // end of find.

    // Finalize the graph storage. Construct CSC, CSRs.
    void finalize(size_t _num_of_v, edge_info &edges) {
      #ifdef DEBUG_GRAPH
      std::cout << "Graph finalize..." << std::endl;
      #endif
      num_vertices = _num_of_v;
      num_edges = edges.size();

      // Permute_index.
      std::vector<edge_id_type>& permute_index = c2r_map;
      permute_index.reserve(num_edges); // aliased as permute index.
      // Counter_index.
      std::vector<atomic<int> > counter_array(num_vertices+1, 0);
      permute_index.assign(num_edges, 0);

      // Sort edges by source;
      // Begin of counting sort.
      #ifdef DEBUG_GRAPH
      std::cout << "Sort by src..." << std::endl;
      #endif
      counting_sort(edges.source_arr, counter_array, permute_index); 
      #ifdef DEBUG_GRAPH 
      std::cout << "finish counting sort." << std::endl;
      #endif
      // Parallel sort target for each source= x interval: counter_array[x] - counter_array[x+1];
#ifndef AVOID_PARALLEL_SORT
#pragma omp parallel for
#else
      logstream(LOG_INFO) << "Parallel sort is disabled." << std::endl;
#endif

      for (ssize_t j = 0; j < ssize_t(num_vertices); ++j) {
        if (counter_array[j] < counter_array[j+1]) {
          std::sort(permute_index.begin()+counter_array[j], 
                    permute_index.begin()+counter_array[j+1],
                    cmp_by_any_functor<vertex_id_type> (edges.target_arr)); 
        }
      }
      // End of counting sort.

      // Inplace permute of edge_data, edge_src, edge_target array.
      // Modified from src/graphlab/util/generics/shuffle.hpp.
      #ifdef DEBUG_GRAPH
      std::cout << "Inplace permute by src..." << std::endl;
      #endif
      vertex_id_type swap_src; vertex_id_type swap_target;
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
      std::cout << "Build CSR_src..." << std::endl;
#endif
      CSR_src.reserve(num_vertices);
      if (use_skip_list) {
        CSR_src_skip.reserve(num_vertices);
      }
      size_t lastSrc = -1;
      vertex_id_type old_src = -1;
      vertex_id_type old_dst = -1;
      // Iterate over the edges. 
      for (size_t it = 0; it < num_edges; ++it) {
        vertex_id_type src = edges.source_arr[it];
        vertex_id_type dst = edges.target_arr[it];
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
        std::cout << "Sort by dst..." << std::endl;
#endif
        counting_sort(edges.target_arr, counter_array, permute_index); 
#ifdef DEBUG_GRAPH
        std::cout << "finish counting sort." << std::endl;
#endif
#ifndef AVOID_PARALLEL_SORT
#pragma omp parallel for
#endif
        for (ssize_t i = 0; i < ssize_t(num_vertices); ++i) {
          if (counter_array[i] < counter_array[i+1]) {
            std::sort(permute_index.begin()+counter_array[i],
                      permute_index.begin() + counter_array[i+1],
                      cmp_by_any_functor<vertex_id_type>(edges.source_arr)); 
          }
        }
        // End of counting sort.

#ifdef DEBUG_GRAPH
        std::cout << "Outplace permute by dst..." << std::endl;
#endif
        outofplace_shuffle(edges.source_arr, permute_index);
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
        std::cout <<"Build CSC_dst..." << std::endl;
#endif
        // Iterate over the edges. 
        for (size_t it = 0; it < num_edges; ++it) {
          vertex_id_type dst = edges.target_arr[c2r_map[it]];

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
      std::cout << "End of finalize." << std::endl;
#endif

      /* DEBUG */
      // printf("CSR dst:\n");
      // foreach(vertex_id_type i, CSR_dst)
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
      // foreach(vertex_id_type i, CSC_dst)
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

    void clear_reserve() {
      std::vector<edge_id_type>().swap(CSR_src);
      std::vector<vertex_id_type>().swap(CSR_dst);
      std::vector<vertex_id_type>().swap(CSC_src);
      std::vector<edge_id_type>().swap(CSC_dst);
      std::vector<edge_id_type>().swap(c2r_map);
      std::vector<EdgeData>().swap(edge_data_list);
      std::vector<vertex_id_type>().swap(CSR_src_skip);
      std::vector<vertex_id_type>().swap(CSC_dst_skip);
    }

    size_t estimate_sizeof() const {
      // const size_t word_size = sizeof(size_t);
      const size_t vid_size = sizeof(vertex_id_type);
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

      std::cout << "CSR size: " 
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
    vertex_id_type source(edge_id_type eid) const {
      ASSERT_LT(eid, num_edges);
      // Helper function: binary search the CSR_row;
      return lookup_source(eid);
    }

    // This is a log(V) operation.
    // Do not use operations on edge id.
    // Use edge_list instead.
    vertex_id_type target(edge_id_type eid) const {
      ASSERT_LT(eid, num_edges);
      return CSR_dst[eid];
    }



    // ------------- Private data storage ----------------
  private:
    size_t num_vertices;
    size_t num_edges;

    // Array for storing edge data, sorted by source vid.
    std::vector<EdgeData> edge_data_list;

    /** Row of CSR */
    std::vector<edge_id_type> CSR_src;

    /* Suppose CSR_src is: 1 x x 3 x x x x 5
     * where x means no out edges.
     *     CSR_src_skip =  0 2 2 0 4 0 0 4 0
     * is used to jump to the prev/next valid vertex in CSR_src.  
     * Optional.
     */
    std::vector<vertex_id_type> CSR_src_skip;

    /** Col of CSR */
    std::vector<vertex_id_type> CSR_dst;

    /** Map the sort by col edge id to sort by row edge id */
    std::vector<edge_id_type> c2r_map;

    /** Col of CSC */
    std::vector<edge_id_type> CSC_dst;

    /* Suppose CSC_dst is: 1 x x 3 x x x x 5
     * where x means no out edges.
     *     CSC_dst_skip =  0 2 2 0 4 0 0 4 0
     * is used to jump to the prev/next valid vertex in CSC_dst.  
     * Optional.
     */
    std::vector<vertex_id_type> CSC_dst_skip;
    /** Row of CSC */
    std::vector<vertex_id_type> CSC_src;

    // graph traits
    bool use_skip_list;

  private:
    // Get the start, end index of the inbound edge of vertex v.
    inline std::pair<bool, edge_range_type> inEdgeRange(vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);

      size_t col_start = CSC_dst[v];
      if (col_start >= num_edges) {
        // No inbound edges.
        return std::make_pair(false, std::make_pair(0,0));
      } else {

        // Find the start column of the next vertex.
        vertex_id_type nextV = use_skip_list ? nextValid(CSC_dst_skip, v, true) : 
          nextValid(CSC_dst, v, false);
        size_t col_end = (nextV < num_vertices) ? CSC_dst[nextV] : num_edges;
        return std::make_pair(true, std::make_pair(col_start, col_end-1));
      }
    } // End of inEdgeRange;

    // Get the start, end index of the outbound edge of vertex v.
    inline std::pair<bool, edge_range_type> outEdgeRange(vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);
      size_t row_start = CSR_src[v];
      if (row_start >= num_edges) {
        // No outbound edges.
        return std::make_pair(false, std::make_pair(0,0));;
      } else {
        // Find the start column of the next vertex.
        vertex_id_type nextV = use_skip_list ? nextValid(CSR_src_skip, v, true) : 
          nextValid(CSR_src, v, false);
        size_t row_end = (nextV < num_vertices) ? CSR_src[nextV] :num_edges; 

        return std::make_pair(true, std::make_pair(row_start, row_end-1));
      }
    } // End of outEdgeRange;


    //-------------Private Helper functions------------
    // Sort by value of vec in the ascending order.
    template <typename anyvalue>
    struct cmp_by_any_functor {
      const std::vector<anyvalue>& vec;
      cmp_by_any_functor(const std::vector<anyvalue>& _vec) : vec(_vec) { }
      bool operator()(size_t me, size_t other) const {
        return (vec[me] < vec[other]);
      }
    };

    template <typename valuetype>
    void counting_sort(const std::vector<valuetype>& value_array, std::vector< atomic<int> >& counter_array, std::vector<edge_id_type>& permute_index) {
      counter_array.assign(counter_array.size(), 0);
      permute_index.assign(permute_index.size(), 0);
#ifndef AVOID_PARALLEL_SORT
#pragma omp parallel for
#endif
      for (ssize_t i = 0; i < ssize_t(value_array.size()); ++i) {
        size_t val = value_array[i];
        counter_array[val].inc();
      }

      for (size_t i = 1; i < counter_array.size(); ++i) {
        counter_array[i] += counter_array[i-1];
      }
#ifndef AVOID_PARALLEL_SORT
#pragma omp parallel for
#endif
      for (ssize_t i = 0; i < ssize_t(value_array.size()); ++i) {
        size_t val = value_array[i];
        permute_index[counter_array[val].dec()] = i;
      }
    }


    size_t binary_search(const std::vector<vertex_id_type>& vec, 
                         size_t start, size_t end, 
                         vertex_id_type vfind) const {
      ASSERT_LT(vfind, num_vertices);
      while(start <= end) {
        size_t mid = (start+end)/2;
        vertex_id_type vpoke = vec[mid];
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

    vertex_id_type
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
            vertex_id_type ans = use_skip_list ? prevValid(CSR_src_skip, start, true) : prevValid(CSR_src, start, false);
            //vertex_id_type ans = start -1;
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
    
    // If use_skip_list, use the Skip array to find the next valid vertex. O(1)
    //
    // Else use the vertex_id (CSR_row or CSC_dst) to find the next
    // vertex. O(V).  Can be applied on invalid vertex in CSR_row or
    // CSC_dst. Useful in binary search where the middle is not
    // assumed to be valid.
    inline vertex_id_type 
    nextValid(const std::vector<vertex_id_type>& vertex_array, 
              vertex_id_type curv, bool use_skip_list) const {

      if (curv == num_vertices-1) return num_vertices;

      if (use_skip_list) {
        return (curv + 1 + vertex_array[curv+1]);
      } else {
        vertex_id_type search = curv+1;
        while (search < num_vertices && vertex_array[search] >= num_edges) 
          ++search;
        return search;
      }
    }

    // Use the Skip array to find the previous valid vertex. O(1)
    inline vertex_id_type 
    prevValid(const std::vector<vertex_id_type>& vertex_array, 
              vertex_id_type curv, bool use_skip_list) const {
      if (curv == 0) return -1;

      if (use_skip_list) {
        return (curv - 1 - vertex_array[curv-1]);
      } else {
        vertex_id_type search = curv-1;
        while (search >= 0 && vertex_array[search] >= num_edges)
          --search;
        return search;
      }
    }

  public:

    const std::vector<vertex_id_type>& get_csr_src() const {
      return CSR_src;
    }
    const std::vector<edge_id_type>& get_csr_dst() const {
      return CSR_dst;
    }
    const std::vector<edge_id_type>& get_csc_src() const {
      return CSC_src;
    }
    const std::vector<vertex_id_type>& get_csc_dst() const {
      return CSC_dst;
    }
    const std::vector<EdgeData>& get_edge_data() const {
      return edge_data_list;
    }

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

  };// End of graph store;
}// End of namespace;
#include <graphlab/macros_undef.hpp>
#endif
