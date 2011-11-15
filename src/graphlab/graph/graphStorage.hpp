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
        typedef std::pair<size_t, size_t>  edge_range_type;
        class edge; class edge_list;

      public:
        graphStorage() { }

                // Finalize the graph storage. Construct CSC, CSRs.
        void finalize(size_t _num_of_v, std::vector<edge> &edges_tmp) {
          num_vertices = _num_of_v;
          num_edges = edges_tmp.size();

          // We sort the edges_tmp in reverse order by row, because it is easy to use pop_back 
          // to transfer the edgedata in to the compact list at the end of this function.
          bool reverse = true;
          std::sort(edges_tmp.begin(), edges_tmp.end(), cmp_by_src_functor(reverse));

          /* DEBUG
          foreach(edge e, edges_tmp)
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
            size_t id = num_edges-1-it;
            edge e = edges_tmp[id];
            vertex_id_type src = e.source();
            vertex_id_type dst = e.target();
            if (src == old_src && dst == old_dst) {
              logstream(LOG_FATAL)
                << "Duplicate edge "
                << id << ":(" << src << ", " << dst << ") "
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
          /* DEBUG 
          printf("CSR dst:\n");
          foreach(vertex_id_type i, CSR_dst)
            std::cout << i << " ";
          std::cout << std::endl;
          printf("CSR src:\n");
          foreach(size_t i, CSR_src)
            std::cout << i << " "; 
          std::cout << std::endl;
          */

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
            edge e = edges_tmp[num_edges -1 -c2r_map[it]];
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
          /*
          printf("CSC src:\n");
          foreach(vertex_id_type i, CSC_src)
            std::cout << i << " ";
          std::cout << std::endl;
          printf("CSR dst:\n");
          foreach(size_t i, CSC_dst)
            std::cout << i << " "; 
          std::cout << std::endl;
          */

          // Transfer the edge data to a compact list.
          edge_data_list.reserve(num_edges);
          for (size_t i = 0; i < num_edges; ++i) {
            edge& e = edges_tmp[num_edges-i-1];
            edge_data_list.push_back(e.data());
            edges_tmp.pop_back();
          }
          edges_tmp.clear();
          std::vector<edge>().swap(edges_tmp);
        }

        void clear() {
          CSR_src.clear();
          CSR_dst.clear();
          CSC_src.clear();
          CSC_dst.clear();
          c2r_map.clear();
          edge_data_list.clear();
        }

        void resetMem() {
          std::vector<size_t>().swap(CSR_src);
          std::vector<vertex_id_type>().swap(CSR_dst);
          std::vector<vertex_id_type>().swap(CSC_src);
          std::vector<size_t>().swap(CSC_dst);
          std::vector<edge_id_type>().swap(c2r_map);
          std::vector<EdgeData>().swap(edge_data_list);
        }

        size_t get_storage_size() const {
          size_t word_size = sizeof(size_t);
          size_t vid_size = sizeof(vertex_id_type);
          size_t eid_size = sizeof(edge_id_type);

          // Actual content size;
          size_t CSR_size = word_size * CSR_src.capacity() + vid_size * CSR_dst.capacity();
          size_t CSC_size = word_size *CSC_dst.capacity() + vid_size * CSC_src.capacity() + eid_size * c2r_map.capacity();
          size_t edata_size = sizeof(EdgeData) * edge_data_list.capacity();

          // Container size;
          size_t container_size = sizeof(CSR_src) + sizeof(CSR_dst) + sizeof(CSC_src) + sizeof(CSC_dst) + sizeof(c2r_map) + sizeof(edge_data_list);

//          printf("CSR_src %u, CSC_dst %u, CSR_dst %u, CSC_src %u, c2r_map %u, e_data_list %u \n", CSR_src.capacity(), CSC_dst.capacity(), CSR_dst.capacity(), CSC_src.capacity(), c2r_map.capacity(), edge_data_list.capacity());
          return CSR_size + CSC_size + edata_size + container_size;
        }

        size_t edge_size() const {
          return num_edges;
        }

        size_t vertices_size() const {
          return num_vertices;
        }

        vertex_id_type  source(edge_id_type eid) const {
          ASSERT_LT(eid, num_edges);
          // Helper function: binary search the CSR_row;
          return lookup_source(eid);
        }

        vertex_id_type  target(edge_id_type eid) const {
          ASSERT_LT(eid, num_edges);
          return CSR_dst[eid];
        }

        size_t num_in_neighbors(vertex_id_type v) const {
          ASSERT_LT(v, num_vertices);
          size_t begin = CSC_dst[v];
          size_t end = (v == num_vertices-1) ? num_edges: CSC_dst[v+1];
          return (end-begin);
        }

        size_t num_out_neighbors(vertex_id_type v) const {
          ASSERT_LT(v, num_vertices);
          size_t begin = CSR_src[v];
          size_t end = (v == num_vertices-1) ? num_edges: CSR_src[v+1];
          return (end-begin);
        }

        std::vector<vertex_id_type> in_vertices(vertex_id_type v) const {
          ASSERT_LT(v, num_vertices);
          std::vector<vertex_id_type> ret;
          size_t begin = CSC_dst[v];
          size_t end = (v == num_vertices-1) ? num_edges: CSC_dst[v+1];
          for (size_t i = begin; i < end; ++i) {
            ret.push_back(CSC_src[i]);
          }
          return ret;
        }

        std::vector<vertex_id_type> out_vertices(vertex_id_type v) const {
          ASSERT_LT(v, num_vertices);
          std::vector<vertex_id_type> ret;
          size_t begin = CSR_src[v];
          size_t end = (v == num_vertices-1) ? num_edges: CSR_src[v+1];
          for (size_t i = begin; i < end; ++i) {
            ret.push_back(CSR_dst[i]);
          }
          return ret;
        }

        edge_list in_edge_ids(vertex_id_type v) const {
          ASSERT_LT(v, num_edges); 
          size_t begin = CSC_dst[v];
          size_t end = (v == num_vertices-1) ? num_edges: CSC_dst[v+1];
          return edge_list(&(c2r_map[begin]), (end-begin));
        }

        edge_list out_edge_ids(vertex_id_type v) const {
          ASSERT_LT(v, num_edges);
          edge_id_type begin = CSR_src[v];
          edge_id_type end = (v == num_vertices-1) ? num_edges: CSR_src[v+1];
          return edge_list(begin, end);
        }

        // Return the edge id given source and destination.
        std::pair<bool, edge_id_type> 
          find (vertex_id_type src, vertex_id_type dst) const {
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
                  return std::make_pair(false, 0);
                } else {
                  return std::make_pair(true, c2r_map[efind]);
                }
              } else {
                // In edge candidate size is smaller, search CSC.
                size_t efind = binary_search(CSR_dst, dstRange.first, dstRange.second, dst);
                if (efind >= num_edges) {
                  return std::make_pair(false, 0);
                } else {
                  return std::make_pair(true, efind);
                }
              }
            } else {
              return std::make_pair(false, 0);
            }
          }

        /* Helper data field and structures: edge_data_list, class edge, class edge_list */
      public:
        std::vector<EdgeData> edge_data_list;
        // Edge class for temporary storage. Will be finalized into the CSR+CSC form.
        class edge {
          public:
            edge() : _source(-1), _target(-1) { }
            edge(const edge& other) :
              _source(other.source()), _target(other.target()),
              _data(other.data()) { }
            edge(vertex_id_type source, vertex_id_type target) :
              _source(source), _target(target)  { }
            edge(vertex_id_type source, vertex_id_type target, EdgeData data) : 
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
        }; // end of edge_data

        // Represents an iteratable list of edge_ids.
        class edge_list {
          public:
            // Define types of union iterator: an shared interface of CSR (lazy) and CSC (normal) iterators.
            typedef boost::counting_iterator<edge_id_type> csrIter; 
            typedef const edge_id_type* cscIter;
            struct unionIterType {
              csrIter lazyIter;
              cscIter normalIter;
            };

            class unionIterator : public std::iterator<std::forward_iterator_tag, edge_id_type> {
              public:
                // Constructors
                unionIterator () : is_lazy(false) {}
                unionIterator (const edge_id_type* p) : is_lazy(false) {
                  uiter.normalIter = p;
                }
                unionIterator (edge_id_type i) : is_lazy(true) {
                  uiter.lazyIter = csrIter(i);
                }
                unionIterator(const unionIterator& it) : 
                  uiter(it.uiter), is_lazy(it.is_lazy){}

                // Interface of an iterator.
                const edge_id_type& operator*() const {
                  ret = is_lazy ? *(uiter.lazyIter) : *(uiter.normalIter);
                  return ret;
                }

                edge_id_type& operator*() {
                  ret = is_lazy ? *(uiter.lazyIter) : *(uiter.normalIter);
                  return ret;
                }

                bool operator==(const unionIterator& it) const {
                  if (it.is_lazy == is_lazy) {
                    return (is_lazy && it.uiter.lazyIter == uiter.lazyIter) || (!is_lazy && it.uiter.normalIter == uiter.normalIter);
                  } else {
                    return false;
                  }
                }

                unionIterator operator++() {
                  if (is_lazy) {
                    ++uiter.lazyIter;
                  } else {
                    ++uiter.normalIter;
                  }
                  return *this;
                }

                int operator-(const unionIterator& it) const {
                  return is_lazy ? uiter.lazyIter - it.uiter.lazyIter : uiter.normalIter - it.uiter.normalIter;
                }

                unionIterator operator+(size_t i) const {
                  unionIterator retval(*this);
                  if (is_lazy) {
                    retval.uiter.lazyIter += i;
                  } else {
                    retval.uiter.normalIter += i;
                  }
                  return retval;
                }

                bool isLazy() {return is_lazy;} 
                // Data fields
              public:
                unionIterType uiter;

              private:
                mutable edge_id_type ret;
                bool is_lazy;
            };

          public:
            // Define types of an iteratable interface.
            typedef unionIterator iterator; // Should not be used
            typedef unionIterator const_iterator;
            typedef edge_id_type value_type;

          public:
            /** \brief Construct an empty edge list */
            edge_list() : begin_ptr(NULL), end_ptr(NULL), is_lazy(false) { }

            /** \brief Construct an edge list from an std vector */
            edge_list(const std::vector<edge_id_type>& edges) :
              is_lazy(false),
              begin_ptr(&(*edges.begin())), 
              end_ptr(begin_ptr + edges.size()) {}

            /** \brief Construct an edge list from an in memory array */
            edge_list(const edge_id_type* begin_ptr, size_t len) :
              is_lazy(false), begin_ptr(begin_ptr),  end_ptr(begin_ptr + len) {} 

            edge_list(const edge_id_type begin, const edge_id_type end) :
              is_lazy(true), begin_ptr(begin), end_ptr(end) {}


            /** \brief Get the size of the edge list */
            size_t size() const { return (size_t)(end_ptr - begin_ptr);}

            /** \brief Get the ith edge in the edge list */
            edge_id_type operator[](size_t i) const {
              ASSERT_LT(i,  size());
              return *(begin_ptr + i);
            }

            /** \brief Returns a pointer to the start of the edge list */
            iterator begin() const {
              return begin_ptr;
            }

            /** \brief Returns a pointer to the end of the edge list */
            iterator end() const {
              return end_ptr;
            } 

            /** \brief Fill a vector with edge id list */
            void fill_vector(std::vector<edge_id_type>& lvalue) const {
              lvalue.clear();
              foreach(edge_id_type eid, *this) lvalue.push_back(eid);    
            }
            /** \brief test if the edge list is empty */
            bool empty() const { return size() == 0; }

            bool isLazy() const {return is_lazy;}

          private:
            bool is_lazy;
            iterator begin_ptr; // Points to first element
            iterator end_ptr; // One past end   
        }; // End of edge list


        // ------------- Private data storage ----------------
      private:
        size_t num_vertices;
        size_t num_edges;

        /** Row of CSR */
        std::vector<size_t> CSR_src;
        /** Col of CSR */
        std::vector<vertex_id_type> CSR_dst;

        /** Map the sort by col edge id to sort by row edge id */
        std::vector<edge_id_type> c2r_map;
        /** Col of CSC */
        std::vector<size_t> CSC_dst;
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
            vertex_id_type nextV = v + 1;
            while(nextV < num_vertices && CSC_dst[nextV] >= num_edges) ++nextV; 
            size_t col_end = (nextV < num_vertices) ? CSC_dst[nextV] : num_vertices;
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
            vertex_id_type nextV = v + 1;
            while(nextV < num_vertices && CSR_src[nextV] >= num_edges) ++nextV; 
            size_t row_end = (nextV < num_vertices) ? CSR_src[nextV] :num_vertices; 
            return std::make_pair(true, std::make_pair(row_start, row_end-1));
          }
        } // End of outEdgeRange;


        //-------------Private Helper functions------------

        // Sort by src in the ascending order.
        struct cmp_by_src_functor {
          bool reverse;
          cmp_by_src_functor (bool order) : reverse(order) {}
          bool operator()(const edge& me, const edge& other) const {
            bool less = (me.source() < other.source()) || (me.source()== other.source()&& me.target()< other.target());
            if (reverse) {
              return !less; 
            } else {
              return less;
            }
          }
        };

        // Sort by dst in the ascending order.
        struct cmp_by_dst_functor {   
          typedef std::vector<edge>& elist_ptr;
          elist_ptr edgelist;
          size_t size;
          cmp_by_dst_functor(elist_ptr edgelist) : edgelist(edgelist) {size = edgelist.size();}

          bool operator()(const size_t me, const size_t other) {
            const edge& e1 = edgelist[size - 1 - me];
            const edge& e2 = edgelist[size - 1 - other];
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
            size_t start = 0;
            size_t end = num_vertices-1;
            // Keep the invariant that CSR_src[start] <= eid, CSR_src[end] > eid.
            while (start <= end) {
              if (CSR_src[end] <= eid)
                return end;
              if (CSR_src[start] == eid)
                return start;
              if (CSR_src[start] > eid)
              {
                ASSERT_LT(0, start);
                return start-1;
              }

              size_t mid = (start + end)/2;
              if (CSR_src[mid] == eid)
                return mid;
              if (CSR_src[mid] > eid) {
                end = mid -1;
              } else {
                start = mid + 1;
              }
            } 
            ASSERT_TRUE(false);
            return start;
          }

        void load(iarchive& arc) {
          clear();
          arc >> num_vertices
            >> num_edges
            >> edge_data_list
            >> CSR_src
            >> CSR_dst
            >> CSC_src
            >> CSC_dst
            >> c2r_map;
        }

        void save(oarchive& arc) {
          arc << num_vertices
            << num_edges
            << edge_data_list
            << CSR_src
            << CSR_dst
            << CSC_src
            << CSC_dst
            << c2r_map;
        }

    };// End of graphStore;
}// End of namespace;
#include <graphlab/macros_undef.hpp>
#endif
