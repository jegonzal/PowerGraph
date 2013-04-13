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
  class graph_storage {
  public:
    typedef graphlab::lvid_type lvid_type;
    typedef graphlab::edge_id_type edge_id_type;

    /** The type of the edge data stored in the graph. */
    typedef EdgeData edge_data_type;

    /** The type of the vertex data stored in the graph. */
    typedef VertexData vertex_data_type;


 public:


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
      
      edge_iterator (lvid_type _center, size_t _offset, 
                     edge_dir_type _itype, const edge_id_type* _vid_arr) :
        center(_center), offset(_offset), itype(_itype), vid_arr(_vid_arr), 
        empty(false) { }
     
      /** \brief Returns the value of the iterator. An empty iterator always returns empty edge type*/ 
      inline edge_type operator*() const  {
      }

      inline typename arrow_type::type operator->() const {
        return arrow_type::make(make_value());
      }

      operator_arrow_dispatch<edge_type, edge_type*>::result_type arrow_type;
      inline arrow_type operator->() const {
        return arrow_type(make_value());
      }
      /** \brief Returns if two iterators point to the same edge. */
      inline bool operator==(const edge_iterator& it) const {
            }

      /** \brief Returns if two iterators don't point to the same edge. */
      inline bool operator!=(const edge_iterator& it) const { 
        return !(*this == it);
      }

      /** \brief Increases the iterator. */
      inline edge_iterator& operator++() {
      }

      /** \brief Increases the iterator. */
      inline edge_iterator operator++(int) {
      }

      /** \brief Computes the difference of two iterators. */
      inline ssize_t operator-(const edge_iterator& it) const {
      }

      /** \brief Returns a new iterator whose value is increased by i difference units. */
      inline edge_iterator operator+(difference_type i) const {
      }

      /** \brief Increases the iterator by i difference units. */
      inline edge_iterator& operator+=(difference_type i) {
        return *this;
      }

      /** \brief Generate the return value of the iterator. */
      inline edge_type make_value() const {
      }

    private:
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
   

    /** \brief Returns the number of edges in the graph. */
    size_t edge_size() const { return num_edges; }

    /** \brief Returns the number of vertices in the graph. */
    size_t vertices_size() const { return num_vertices; }

    /** \brief Returns the number of in edges of the vertex. */
    size_t num_in_edges (const lvid_type v) const {
    }

    /** \brief Returns the number of out edges of the vertex. */
    size_t num_out_edges (const lvid_type v) const {
    }

    /** \brief Returns the edge id of the edge. 
     * Edges are assigned with consecutive int ids,
     * ordered first by source and then by target.
     * */
    edge_id_type edge_id(const edge_type& edge) const {
    }

    /** \brief Returns the reference of edge data of an edge. */
    edge_data_type& edge_data(lvid_type source, lvid_type target) {
    }

    /** \brief Returns the constant reference of edge data of an edge. */
    const edge_data_type& edge_data(lvid_type source, 
                                    lvid_type target) const {
    }

    /** \brief Returns the reference of edge data of an edge. */
    edge_data_type& edge_data(edge_type edge) {
    }

    /** \brief Returns the constant reference of edge data of an edge. */
    const edge_data_type& edge_data(edge_type edge) const {
        }

    /** \brief Returns a list of in edges of a vertex. */
    edge_list in_edges(const lvid_type v) const {
    }

    /** \brief Returns a list of out edges of a vertex. */
    edge_list out_edges(const lvid_type v) const {
    }

    /** \brief Returns an edge type of a given source target pair. */
    edge_type find (const lvid_type src, 
                    const lvid_type dst) const {
    } // end of find.

     /** \brief Finalize the graph storage. 
      */
    void finalize(size_t _num_of_v, edge_info &edges){  
    } // end of finalize.

    /** \brief Reset the storage. */
    void clear() {
    }


    // ------------- Private data storage ----------------
  private:
    /** Number of vertices in the storage (not counting singletons)*/
    size_t num_vertices;
    /** Number of edges in the storage. */
    size_t num_edges;

    struct row {
    };

    std::vector<std::vector> CSC_dst;



  public:
    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();
    }

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
    }

    /** swap two graph storage*/
    void swap(graph_storage& other) {
    }

  };// End of graph store;
}// End of namespace;

#include <graphlab/macros_undef.hpp>
#endif
