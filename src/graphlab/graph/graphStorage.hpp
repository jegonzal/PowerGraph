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

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {

  template<typename VertexData, typename EdgeData>
  class graphStorage {
  public:
    typedef uint32_t vertex_id_type;
    typedef uint32_t edge_id_type;
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
      edge_info() : _source(-1), _target(-1) { }
      edge_info(const edge_info& other) :
        _source(other.source()), _target(other.target()),
        _data(other.data()) { }
      edge_info(vertex_id_type source, vertex_id_type target) :
        _source(source), _target(target)  { }
      edge_info(vertex_id_type source, vertex_id_type target, EdgeData data) : 
        _source(source), _target(target), _data(data) {}
    public:
      inline vertex_id_type source() const { return _source; }
      inline vertex_id_type target() const { return _target; }   
      inline EdgeData& data() { return _data; }
      inline const EdgeData& data() const { return _data; }
    private:
      vertex_id_type _source;
      vertex_id_type _target;
      EdgeData _data;
    }; // end of class edge.

    // A class of edge information. Used as value type of the edge_list.
    class edge_type {
    public:
      edge_type () : _source(-1), _target(-1), _edge_id(-1), _empty(true) { }
      edge_type (const vertex_id_type _source, const vertex_id_type _target, 
                 const edge_id_type _eid) :
        _source(_source), _target(_target), _edge_id(_eid), _empty(false) { }
    public:
      inline vertex_id_type source() const {
        ASSERT_FALSE(empty()); 
        return _source; 
      }

      inline vertex_id_type target() const { 
        ASSERT_FALSE(empty());
        return _target; 
      }
      inline bool empty() const { return _empty; }
      // Data fields. 
    private:
      vertex_id_type _source;
      vertex_id_type _target;
      edge_id_type _edge_id;
      bool _empty;

      friend class graphStorage;
    }; // end of class edge_type.

    // Internal iterator on edge_types.
    class edge_iterator : public std::iterator<std::forward_iterator_tag, edge_type> {
    public:
      enum iterator_type{INEDGE, OUTEDGE}; 
      typedef edge_type reference;
    public:
      // Cosntructors
      edge_iterator () : offset(-1), empty(true) { }
      edge_iterator (vertex_id_type _center, size_t _offset, iterator_type _itype, const graphStorage* _gstore_ptr) : 
        center(_center), offset(_offset), itype(_itype), empty(false), gstore_ptr(_gstore_ptr) { }
      edge_iterator (const edge_iterator& it) :
        center(it.center), offset(it.offset), itype(it.itype), empty(it.empty), gstore_ptr(it.gstore_ptr) { }
    public:

      inline edge_type operator*() const {
        ASSERT_TRUE(!empty);
        return makeValue();
      }

      inline bool operator==(const edge_iterator& it) const {
        if (empty && it.empty)
          return true;
        if (empty != it.empty)
          return false;
        return (itype == it.itype && center == it.center && offset == it.offset);
      }

      inline edge_iterator operator++() {
        ASSERT_TRUE(!empty);
        ++offset;
        return *this;
      }

      inline int operator-(const edge_iterator& it) const {
        ASSERT_TRUE(!empty && itype == it.itype && center == it.center);
        return offset - it.offset;
      }

      inline edge_iterator operator+(size_t i) const {
        edge_iterator retval(center, offset+i, itype, gstore_ptr);
        return retval;
      }

    private:
      // Generate the ret value of the iterator.
      inline edge_type makeValue() const {
        edge_type ret;
        if (itype == INEDGE) {
          edge_type rvalue(gstore_ptr->CSC_src[offset], center, gstore_ptr->c2r_map[offset]);
          ret = rvalue;
        } else if (itype == OUTEDGE) {
          edge_type rvalue(center, gstore_ptr->CSR_dst[offset], offset);
          ret = rvalue;
        } else {
          logstream(LOG_FATAL) << "Edge iterator type is invalid." << std::endl;
        }
        return ret;
      }
    private:
      vertex_id_type center;
      size_t offset;
      iterator_type itype;
      bool empty;
      const graphStorage* gstore_ptr;
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
      edge_list(edge_iterator begin, edge_iterator end) : begin_ptr(begin), end_ptr(end) { 
        list_size = (size_t)(end_ptr-begin_ptr);
      }
      // Copy constructor
      edge_list(const edge_list& other) : begin_ptr(other.begin_ptr), end_ptr(other.end_ptr), list_size(other.list_size) { }

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
    graphStorage() { }

    size_t edge_size() const {
      return num_edges;
    }

    size_t vertices_size() const {
      return num_vertices;
    }

    // Return the size of in_neighbours.
    size_t num_in_edges (const vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);
      size_t begin = CSC_dst[v];
      if (begin >= num_edges) return 0;

      // Search is the next valid src vertex after v.
      size_t search = nextValid(CSC_dst_skip, v);
      size_t end = (search >= num_vertices) ? num_edges: CSC_dst[search];
      return (end-begin);
    }

    // Return the size of out_neighbours.
    size_t num_out_edges (const vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);
      size_t begin = CSR_src[v];
      if (begin >= num_edges) return 0;

      size_t search = nextValid(CSR_src_skip, v);
      size_t end = (search >= num_vertices) ? num_edges: CSR_src[search];
      return (end-begin);
    }

    edge_data_type& edge_data(vertex_id_type source, vertex_id_type target) {
      ASSERT_LT(source, num_vertices);
      ASSERT_LT(target, num_vertices);
      edge_type ans = find(source, target);
      return edge_data(ans);
    }

    const edge_data_type& edge_data(vertex_id_type source, vertex_id_type target) const {
      ASSERT_LT(source, num_vertices);
      ASSERT_LT(target, num_vertices);
      edge_type ans = find(source, target);
      return edge_data(ans);
    }

    edge_data_type& edge_data(edge_type edge) {
      ASSERT_FALSE(edge.empty());
      return edge_data_list[edge._edge_id];
    }

    const edge_data_type& edge_data(edge_type edge) const {
      ASSERT_FALSE(edge.empty());
      return edge_data_list[edge._edge_id];
    }


    // Return in edge list of a vertex.
    edge_list in_edges(const vertex_id_type v) const {
      std::pair<bool, edge_range_type> rangePair = inEdgeRange(v);
      if (rangePair.first) {
        edge_range_type range = rangePair.second;
        edge_iterator begin (v, range.first, edge_iterator::INEDGE, this);
        edge_iterator end (v, range.second+1, edge_iterator::INEDGE, this);
        // std::cout << "in range (" << range.first << "," << range.second << ")" << std::endl;
        // std::cout << "in edge size: " << end-begin << std::endl;
        return edge_list(begin, end);
      } else {
        return edge_list();
      }
    }

    // Return out edge list of a vertex.
    edge_list out_edges(const vertex_id_type v) const {
      std::pair<bool, edge_range_type> rangePair = outEdgeRange(v);
      if (rangePair.first) {
        edge_range_type range = rangePair.second;
        edge_iterator begin (v, range.first, edge_iterator::OUTEDGE, this);
        edge_iterator end (v, range.second+1, edge_iterator::OUTEDGE, this);

        // std::cout << "out range (" << range.first << "," << range.second << ")" << std::endl;
        // std::cout << "out_edge size: " << end-begin << std::endl;
        return edge_list(begin, end);
      } else {
        return edge_list();
      }
    }

    // Return the edge id given source and destination.
    edge_type find (const vertex_id_type src, 
                    const vertex_id_type dst) const {
      /* DEBUG
         printf("Find: %u, %u \n", src, dst);
      */
      // Get the out edge range of the src, as well as the in edge range of the dst. 
      std::pair<bool, edge_range_type> srcRangePair = inEdgeRange(dst);
      std::pair<bool, edge_range_type> dstRangePair = outEdgeRange(src);

      if( srcRangePair.first && dstRangePair.first) {
        // The edge may exist. 
        edge_range_type srcRange =  srcRangePair.second;
        edge_range_type dstRange = dstRangePair.second;

        if ((srcRange.second - srcRange.first) < (dstRange.second-dstRange.first)) {
          // Out edge candidate size is smaller, search CSR.
          size_t efind =  binary_search(CSC_src, srcRange.first, srcRange.second, src);
          if (efind >= num_edges) {
            return edge_type();
          } else {
            return edge_type(src, dst, c2r_map[efind]);
          }
        } else {
          // In edge candidate size is smaller, search CSC.
          size_t efind = binary_search(CSR_dst, dstRange.first, dstRange.second, dst);
          if (efind >= num_edges) {
            return edge_type();
          } else {
            return edge_type(src, dst, efind);
          }
        }
      } else {
        return edge_type();
      }
    } // end of find.

    // Finalize the graph storage. Construct CSC, CSRs.
    void finalize(size_t _num_of_v, std::vector<edge_info> &edges_tmp) {
      num_vertices = _num_of_v;
      num_edges = edges_tmp.size();

      // We sort the edges_tmp in reverse order by row, because it is easy to use pop_back 
      // to transfer the edgedata in to the compact list at the end of this function.
      std::sort(edges_tmp.begin(), edges_tmp.end(), cmp_by_src_functor());

      /* DEBUG
         foreach(edge_info e, edges_tmp)
         std::cout << "(" << e.source() << "," << e.target() << ") ";
         std::cout << std::endl;
      */

      // Construct CSR:
      CSR_src.reserve(num_vertices);
      CSR_dst.reserve(num_edges);

      size_t colIdx = 0;
      size_t lastSrc = -1;

      vertex_id_type old_src = -1;
      vertex_id_type old_dst = -1;

      // Iterate over the edges_tmp in the reverse order.
      for (size_t it = 0; it < num_edges; ++it) {
        edge_info e = edges_tmp[it];
        vertex_id_type src = e.source();
        vertex_id_type dst = e.target();
        if (src == old_src && dst == old_dst) {
          logstream(LOG_FATAL)
            << "Duplicate edge "
            << it << ":(" << src << ", " << dst << ") "
            << "found! Graphlab does not support graphs "
            << "with duplicate edges." << std::endl;
        } else {
          old_src = src;
          old_dst = dst;
        }
        // Fill the columns index list;
        CSR_dst.push_back(dst);
        // Fill the row index list;
        if (src != lastSrc) {
          for (size_t j = (lastSrc+1); j < src; ++j) {
            CSR_src.push_back(-1);
          }
          CSR_src.push_back(colIdx);
          lastSrc = src;
        }
        ++colIdx;
      }
      // Fill in the remaining row index list.
      for( size_t j = (lastSrc +1); j < num_vertices; ++j) {
        CSR_src.push_back(-1);
      }
      ASSERT_EQ(CSR_src.size(), num_vertices);
      /* DEBUG */
      // printf("CSR dst:\n");
      // foreach(vertex_id_type i, CSR_dst)
      //   std::cout << i << " ";
      // std::cout << std::endl;
      // printf("CSR src:\n");
      // foreach(size_t i, CSR_src)
      //   std::cout << i << " "; 
      // std::cout << std::endl;

      // Construct c2r map, sort the ids according to column first order.
      c2r_map.reserve(num_edges);
      for(size_t i = 0; i < num_edges; ++i) {
        edge_id_type eid = (edge_id_type)(i);
        c2r_map.push_back(eid);
      }
      // Functor used to sort ids.
      cmp_by_dst_functor cmp_by_dst(edges_tmp);
      std::sort(c2r_map.begin(), c2r_map.end(), cmp_by_dst);

      /* DEBUG
         printf("c2r_map: \n");
         foreach(edge_id_type e, c2r_map)
         std::cout << e << " ";
         std::cout << std::endl;
      */

      // Construct CSC:
      CSC_dst.reserve(num_vertices);
      CSC_src.reserve(num_edges);

      size_t rowIdx = 0;
      size_t lastDst = -1;

      // Iterate over the edges_tmp
      for (size_t it = 0; it < num_edges; ++it) {
        edge_info e = edges_tmp[c2r_map[it]];
        vertex_id_type src = e.source();
        vertex_id_type dst = e.target();

        // Fill the columns index list;
        CSC_src.push_back(src);

        // Fill the row index list;
        if (dst != lastDst) {
          for (size_t j = (lastDst+1); j < dst; ++j) {
            CSC_dst.push_back(-1);
          }
          CSC_dst.push_back(rowIdx);
          lastDst = dst;
        }
        ++rowIdx;
      }
      // Fill in the remaining col index list.
      for( size_t j = (lastDst +1); j < num_vertices; ++j) {
        CSC_dst.push_back(-1);
      }
      ASSERT_EQ(CSC_dst.size(), num_vertices);
      /*DEBUG*/
      // printf("CSC src:\n");
      // foreach(vertex_id_type i, CSC_src)
      //   std::cout << i << " ";
      // std::cout << std::endl;
      // printf("CSR dst:\n");
      // foreach(size_t i, CSC_dst)
      //   std::cout << i << " "; 
      // std::cout << std::endl;

      // Construct the CSR_src_skip, CSC_dst_skip array for fast iteration. 
      CSR_src_skip.assign(num_vertices, 0);
      CSC_dst_skip.assign(num_vertices, 0);
      size_t iter = 0;
      while (iter < num_vertices) {
        if (CSR_src[iter] >= num_edges) 
          {
            size_t begin_invalid_block = iter;
            size_t end_invalid_block = iter;
            while (end_invalid_block < num_vertices && CSR_src[end_invalid_block] >= num_edges)
              ++end_invalid_block;
            size_t block_size = end_invalid_block - begin_invalid_block;
            CSR_src_skip[begin_invalid_block] = block_size;
            CSR_src_skip[end_invalid_block-1] = block_size;
            iter = end_invalid_block;
          } else {
          ++iter;
        }
      }
      iter = 0;
      while (iter < num_vertices) {
        if (CSC_dst[iter] >= num_edges) 
          {
            size_t begin_invalid_block = iter;
            size_t end_invalid_block = iter;
            while (end_invalid_block < num_vertices && CSC_dst[end_invalid_block] >= num_edges)
              ++end_invalid_block;
            size_t block_size = end_invalid_block - begin_invalid_block;
            CSC_dst_skip[begin_invalid_block] = block_size;
            CSC_dst_skip[end_invalid_block-1] = block_size;
            iter = end_invalid_block;
          } else {
          ++iter;
        }
      }
      // std::cout << "Skip list size: " << double(CSR_src_skip.capacity() * sizeof(vertex_id_type)* 2) / (1024*1024) << "MB"<<std::endl;
      // printf("CSR skip:\n");
      // foreach(size_t i, CSR_src_skip)
      //   std::cout << i << " ";
      // std::cout << std::endl;
      // printf("CSC skip:\n");
      // foreach(size_t i, CSC_dst_skip)
      //   std::cout << i << " "; 
      // std::cout << std::endl;

          
      // Transfer the edge data to a compact list.
      edge_data_list.reserve(num_edges);
      for (size_t i = 0; i < num_edges; ++i) {
        edge_info& e = edges_tmp[i];
        edge_data_list.push_back(e.data());
      }
      std::vector<edge_info>().swap(edges_tmp);
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
      std::vector<size_t>().swap(CSR_src);
      std::vector<vertex_id_type>().swap(CSR_dst);
      std::vector<vertex_id_type>().swap(CSC_src);
      std::vector<size_t>().swap(CSC_dst);
      std::vector<edge_id_type>().swap(c2r_map);
      std::vector<EdgeData>().swap(edge_data_list);
      std::vector<vertex_id_type>().swap(CSR_src_skip);
      std::vector<vertex_id_type>().swap(CSC_dst_skip);
    }

    size_t estimate_sizeof() const {
      size_t word_size = sizeof(size_t);
      size_t vid_size = sizeof(vertex_id_type);
      size_t eid_size = sizeof(edge_id_type);

      // Actual content size;
      size_t CSR_size = word_size * CSR_src.capacity() + vid_size * CSR_dst.capacity();
      size_t CSC_size = word_size *CSC_dst.capacity() + vid_size * CSC_src.capacity() + eid_size * c2r_map.capacity();
      size_t edata_size = sizeof(EdgeData) * edge_data_list.capacity();

      // Container size;
      size_t container_size = sizeof(CSR_src) + sizeof(CSR_dst) + sizeof(CSC_src) + sizeof(CSC_dst) + sizeof(c2r_map) + sizeof(edge_data_list);

      // Skip list size:
      size_t skip_list_size = sizeof(CSR_src_skip) + sizeof(CSC_dst_skip) + CSR_src_skip.capacity() * vid_size + CSC_dst_skip.capacity() * vid_size;
      return CSR_size + CSC_size + edata_size + container_size + skip_list_size;;
    }

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
    std::vector<size_t> CSR_src;

    /* Suppose CSR_src is: 1 x x 3 x x x x 5
     * where x means no out edges.
     *     CSR_src_skip =  0 2 2 0 4 0 0 4 0
     * is used to jump to the prev/next valid vertex in CSR_src.  
     */
    std::vector<vertex_id_type> CSR_src_skip;

    /** Col of CSR */
    std::vector<vertex_id_type> CSR_dst;

    /** Map the sort by col edge id to sort by row edge id */
    std::vector<edge_id_type> c2r_map;
    /** Col of CSC */
    std::vector<size_t> CSC_dst;

    /* Suppose CSC_dst is: 1 x x 3 x x x x 5
     * where x means no out edges.
     *     CSC_dst_skip =  0 2 2 0 4 0 0 4 0
     * is used to jump to the prev/next valid vertex in CSC_dst.  
     */
    std::vector<vertex_id_type> CSC_dst_skip;
    /** Row of CSC */
    std::vector<vertex_id_type> CSC_src;

  private:
    // Get the start, end index of the inbound edge of vertex v.
    std::pair<bool, edge_range_type> inEdgeRange(vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);

      size_t col_start = CSC_dst[v];
      if (col_start >= num_edges) {
        // No inbound edges.
        return std::make_pair(false, std::make_pair(0,0));
      } else {

        // Find the start column of the next vertex.
        vertex_id_type nextV = nextValid(CSC_dst_skip, v);
        size_t col_end = (nextV < num_vertices) ? CSC_dst[nextV] : num_edges;
        return std::make_pair(true, std::make_pair(col_start, col_end-1));
      }
    } // End of inEdgeRange;

    // Get the start, end index of the outbound edge of vertex v.
    std::pair<bool, edge_range_type> outEdgeRange(vertex_id_type v) const {
      ASSERT_LT(v, num_vertices);
      size_t row_start = CSR_src[v];
      if (row_start >= num_edges) {
        // No outbound edges.
        return std::make_pair(false, std::make_pair(0,0));;
      } else {
        // Find the start column of the next vertex.
        vertex_id_type nextV = nextValid(CSR_src_skip, v);
        size_t row_end = (nextV < num_vertices) ? CSR_src[nextV] :num_edges; 

        return std::make_pair(true, std::make_pair(row_start, row_end-1));
      }
    } // End of outEdgeRange;


    //-------------Private Helper functions------------

    // Sort by src in the ascending order.
    struct cmp_by_src_functor {
      bool operator()(const edge_info& me, const edge_info& other) const {
        bool less = (me.source() < other.source()) || (me.source()== other.source()&& me.target()< other.target());
        /*
          if (reverse) {
          return !less; 
          } else {
          return less;
          }
        */
        return less;
      }
    };

    // Sort by dst in the ascending order.
    struct cmp_by_dst_functor {   
      typedef std::vector<edge_info>& elist_ptr;
      elist_ptr edgelist;
      size_t size;
      cmp_by_dst_functor(elist_ptr edgelist) : edgelist(edgelist) {size = edgelist.size();}

      bool operator()(const size_t me, const size_t other) {
        const edge_info& e1 = edgelist[me];
        const edge_info& e2 = edgelist[other];
        return (e1.target() < e2.target()) || (e1.target() == e2.target() && e1.source() < e2.source());
      }
    };

    size_t binary_search(const std::vector<vertex_id_type>& vec, size_t start, size_t end, vertex_id_type vfind) const {
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
      if (CSR_src[start] >= num_edges) start = nextValid(CSR_src_skip, start);
      if (CSR_src[end] >= num_edges) end = prevValid(CSR_src_skip, end);

      // Keep the invariant that CSR_src[start] <= eid, CSR_src[end] > eid.
      while (start <= end) {
        if (CSR_src[end] <= eid)
          return end;
        if (CSR_src[start] == eid)
          return start;
        if (CSR_src[start] > eid)
          {
            ASSERT_LT(0, start);
            vertex_id_type ans = prevValid(CSR_src_skip, start);
            //vertex_id_type ans = start -1;
            ASSERT_LT(ans, num_vertices);
            return ans;
          }

        size_t mid = (start + end)/2;
        // mid may fall in to an invalid grid 
        if (CSR_src[mid] >= num_edges) mid = prevValid2(CSR_src, mid);

        if (CSR_src[mid] == eid)
          return mid;

        if (CSR_src[mid] > eid) {
          end = prevValid(CSR_src_skip, mid); 
          //end = mid-1;
        } else {
          //start = mid+1;
          start = nextValid(CSR_src_skip, mid);
        }
      } 
            
      ASSERT_TRUE(false);
      return start;
    }

    // Use the Skip array to find the next valid vertex. O(1)
    inline vertex_id_type nextValid(const std::vector<vertex_id_type>& vertex_skip, vertex_id_type curv) const {
      if (curv == num_vertices-1) return num_vertices;
      return (curv + 1 + vertex_skip[curv+1]);
    }

    // Use the Skip array to find the previous valid vertex. O(1)
    inline vertex_id_type prevValid(const std::vector<vertex_id_type>& vertex_skip, vertex_id_type curv) const {
      if (curv == 0) return -1;
      return (curv - 1 - vertex_skip[curv-1]);
    }

    // Use the vertex_id (CSR_row or CSC_dst) to find the next vertex. O(V).
    // Can be applied on invalid vertex in CSR_row or CSC_dst. Useful in binary search where the middle is not
    // assumed to be valid.
    inline vertex_id_type nextValid2(const std::vector<size_t>& vertex_array, vertex_id_type curv) const {
      vertex_id_type search = curv+1;
      while (search < num_vertices && vertex_array[search] >= num_edges) ++search;
      return search;
    }

    // Use the vertex_id (CSR_row or CSC_dst) to find the next vertex. O(V).
    // Similar to nextValid2.
    inline vertex_id_type prevValid2(const std::vector<size_t>& vertex_array, vertex_id_type curv) const {
      vertex_id_type search = curv-1;
      while (search >= 0 && vertex_array[search] >= num_edges) --search;
      return search;
    }

  public:

    void load(iarchive& arc) {
      clear();
      arc >> num_vertices
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
      arc << num_vertices
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
  };// End of graphStore;
}// End of namespace;
#include <graphlab/macros_undef.hpp>
#endif
