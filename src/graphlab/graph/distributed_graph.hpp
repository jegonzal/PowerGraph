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

#ifndef GRAPHLAB_DISTRIBUTED_GRAPH_HPP
#define GRAPHLAB_DISTRIBUTED_GRAPH_HPP

#ifndef __NO_OPENMP__
#include <omp.h>
#endif

#include <cmath>

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <graphlab/util/dense_bitset.hpp>


#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>
#include <sstream>

#include <boost/functional.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>


#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/branch_hints.hpp>
#include <graphlab/util/generics/conditional_addition_wrapper.hpp>

#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/graph/local_graph.hpp>
#include <graphlab/graph/ingress/idistributed_ingress.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/ingress/distributed_batch_ingress2.hpp>
#include <graphlab/graph/ingress/distributed_oblivious_ingress.hpp>
#include <graphlab/graph/ingress/distributed_random_ingress.hpp>
#include <graphlab/graph/ingress/distributed_identity_ingress.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <graphlab/util/cuckoo_map_pow2.hpp>

#include <graphlab/util/fs_util.hpp>
#include <graphlab/util/hdfs.hpp>


#include <graphlab/graph/builtin_parsers.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab { 

  /**
   * test
   */
  template<typename VertexData, typename EdgeData>
  class distributed_graph {
  public:

    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    /**
       The line parse is any function (or functor) that has the form:
     
       <code>
        bool line_parser(distributed_graph& graph, const std::string& filename,
                         const std::string& textline);
       </code>

       the line parser returns true if the line is parsed successfully and
       calls graph.add_vertex(...) or graph.add_edge(...)
     */
    typedef boost::function<bool(distributed_graph&, const std::string&,
                                 const std::string&)> line_parser_type;


    typedef fixed_dense_bitset<RPC_MAX_N_PROCS> mirror_type;

    /// The type of the local graph used to store the graph data 
    typedef graphlab::local_graph<VertexData, EdgeData> local_graph_type;

    friend class distributed_ingress_base<VertexData, EdgeData>;
    friend class distributed_random_ingress<VertexData, EdgeData>;
    friend class distributed_identity_ingress<VertexData, EdgeData>;
    friend class distributed_batch_ingress<VertexData, EdgeData>;
    friend class distributed_oblivious_ingress<VertexData, EdgeData>;

    typedef graphlab::vertex_id_type vertex_id_type;
    typedef graphlab::lvid_type lvid_type;
    typedef graphlab::edge_id_type edge_id_type;
    
    struct vertex_type;
    typedef bool edge_list_type;  
    class edge_type;
    
    struct local_vertex_type;
    struct local_edge_list_type;
    class local_edge_type;
    
    /** Vertex object which provides access to the vertex data
     * and information about it.
     */
    struct vertex_type {
      distributed_graph& graph_ref;
      lvid_type lvid;
      vertex_type(distributed_graph& graph_ref, lvid_type lvid):
            graph_ref(graph_ref), lvid(lvid) { }

      bool operator==(vertex_type& v) const {
        return lvid == v.lvid;
      }
      
      /// \brief Returns a constant reference to the data on the vertex
      const vertex_data_type& data() const {
        return graph_ref.get_local_graph().vertex_data(lvid);
      }

      /// \brief Returns a reference to the data on the vertex
      vertex_data_type& data() {
        return graph_ref.get_local_graph().vertex_data(lvid);
      }

      /// \brief Returns the number of in edges of the vertex
      size_t num_in_edges() const {
        return graph_ref.l_get_vertex_record(lvid).num_in_edges;
      }

      /// \brief Returns the number of out edges of the vertex
      size_t num_out_edges() const {
        return graph_ref.l_get_vertex_record(lvid).num_out_edges;
      }
      
      /// \brief Returns the vertex ID of the vertex       
      vertex_id_type id() const {
        return graph_ref.global_vid(lvid);
      }
 
      /// \cond GRAPHLAB_INTERNAL     
      /// \brief Returns a list of in edges (not implemented) 
      edge_list_type in_edges() __attribute__ ((noreturn)) {
        ASSERT_TRUE(false);
      }

      /// \brief Returns a list of out edges (not implemented) 
      edge_list_type out_edges() __attribute__ ((noreturn)) {
        ASSERT_TRUE(false);
      }
      /// \endcond 
      
      /** 
       *  \brief Returns the local ID of the vertex
       */
      lvid_type local_id() const {
        return lvid;
      }
      
    };

    
    /** Edge object which provides access to a single edge
        on the graph */
    class edge_type {
    private:
      distributed_graph& graph_ref;
      typename local_graph_type::edge_type e;
    public:
      edge_type(distributed_graph& graph_ref,
                typename local_graph_type::edge_type e):
                                          graph_ref(graph_ref), e(e) { }

      /// \brief Returns the source vertex of the edge
      vertex_type source() { return vertex_type(graph_ref, e.source().id()); }
      /// \brief Returns the target vertex of the edge
      vertex_type target() { return vertex_type(graph_ref, e.target().id()); }
      
      /// \brief Returns a constant reference to the data on the edge 
      const edge_data_type& data() const { return e.data(); }
      
      /// \brief Returns a reference to the data on the edge 
      edge_data_type& data() { return e.data(); }
    }; 

    // CONSTRUCTORS ==========================================================>
    distributed_graph(distributed_control& dc, 
                      const graphlab_options& opts = graphlab_options() ) : 
      rpc(dc, this), finalized(false), vid2lvid(-1),
      nverts(0), nedges(0), local_own_nverts(0), nreplicas(0),
      ingress_ptr(NULL) {
      rpc.barrier();

      set_options(opts);
    }


    void set_options(const graphlab_options& opts) {
      size_t bufsize = 50000;
      bool usehash = false;
      bool userecent = false;
      std::string ingress_method = "random";
      std::vector<std::string> keys = opts.get_graph_args().get_option_keys();
      foreach(std::string opt, keys) {
        if (opt == "ingress") {
          opts.get_graph_args().get_option("ingress", ingress_method);
        } else if (opt == "bufsize") {
          opts.get_graph_args().get_option("bufsize", bufsize);
        } else if (opt == "usehash") {
          opts.get_graph_args().get_option("usehash", usehash);
        } else if (opt == "userecent") {
          opts.get_graph_args().get_option("userecent", userecent);
        } else {
          logstream(LOG_ERROR) << "Unexpected Graph Option: " << opt << std::endl;
        }
      }
      set_ingress_method(ingress_method, bufsize, usehash, userecent);
    }

    // METHODS ===============================================================>
    /**
     * Finalize is used to complete graph ingress by resolving vertex
     * ownship and completing local data structures.
     */
    void finalize() {
      if (finalized) return;
      ASSERT_NE(ingress_ptr, NULL);
      logstream(LOG_INFO) << "Distributed graph: enter finalize" << std::endl;
      ingress_ptr->finalize();
      rpc.barrier(); delete ingress_ptr; ingress_ptr = NULL;
      finalized = true;
    }
    
    bool is_finalized() {
      return finalized;
    }
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const { return nverts; }

    /** \brief Get the number of edges */
    size_t num_edges() const { return nedges; }

    /// \brief converts a vertex ID to a vertex object
    vertex_type vertex(vertex_id_type vid) {
      return vertex_type(*this, local_vid(vid));
    }


    /** \brief Get a list of all in edges of a given vertex ID. Not Implemented */
    edge_list_type in_edges(const vertex_id_type vid) const 
      __attribute__((noreturn)) {
      // Not implemented.
      logstream(LOG_WARNING) << "in_edges not implemented. " << std::endl;
      ASSERT_TRUE(false);
    }


    /**
     * \brief Returns the number of in edges of a given vertex ID.
     * 
     * Returns the number of in edges of a given vertex ID.
     * Equivalent to vertex(vid).num_in_edges()
     */
    size_t num_in_edges(const vertex_id_type vid) const {
      return get_vertex_record(vid).num_in_edges;
    }

    /** Get a list of all out edges of a given vertex ID. Not Implemented */
    edge_list_type out_edges(const vertex_id_type vid) const 
      __attribute__((noreturn)) {
            // Not implemented.
      logstream(LOG_WARNING) << "in_edges not implemented. " << std::endl;
      ASSERT_TRUE(false);
    }

    /**
     * \brief Returns the number of out edges of a given vertex ID.
     * 
     * Returns the number of in edges of a given vertex ID.
     * Equivalent to vertex(vid).num_out_edges()
     */
    size_t num_out_edges(const vertex_id_type vid) const {
      return get_vertex_record(vid).num_out_edges;
    }
   
    /** 
     * \brief Creates a vertex containing the vertex data
     */
    void add_vertex(const vertex_id_type& vid, 
                    const VertexData& vdata = VertexData() ) {
      ASSERT_NE(ingress_ptr, NULL);
      if(vid == vertex_id_type(-1)) {
        logstream(LOG_FATAL)
          << "\n\tAdding a vertex with id -1 is not allowed."
          << "\n\tThe -1 vertex id is reserved for internal use."
          << std::endl;
      }
      ingress_ptr->add_vertex(vid, vdata);
    }

    /**
     * \brief Creates an edge connecting vertex source to vertex target.  
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const EdgeData& edata = EdgeData()) {
      ASSERT_NE(ingress_ptr, NULL);
      if(source == vertex_id_type(-1)) {
        logstream(LOG_FATAL)
          << "\n\tThe source vertex with id vertex_id_type(-1)\n"
          << "\tor unsigned value " << vertex_id_type(-1) << " in edge \n"
          << "\t(" << source << "->" << target << ") is not allowed.\n" 
          << "\tThe -1 vertex id is reserved for internal use."
          << std::endl;
      }
      if(target == vertex_id_type(-1)) {
        logstream(LOG_FATAL)
          << "\n\tThe target vertex with id vertex_id_type(-1)\n"
          << "\tor unsigned value " << vertex_id_type(-1) << " in edge \n"
          << "\t(" << source << "->" << target << ") is not allowed.\n" 
          << "\tThe -1 vertex id is reserved for internal use."
          << std::endl;
      }
      if(source == target) {
        logstream(LOG_FATAL)
          << "\n\tTrying to add self edge (" << source << "->" << target << ")."
          << "\n\tSelf edges are not allowed."
          << std::endl;
      }
      ASSERT_NE(source, target);
      ingress_ptr->add_edge(source, target, edata);
    }


   /**
    * \brief Performs a map-reduce operation on each vertex in the 
    * graph returning the result.
    * 
    * Given a map function, map_reduce_vertices() call the map function on all
    * vertices in the graph. The return values are then summed together and the
    * final result returned. The map function should only read the vertex data
    * and should not make any modifications. map_reduce_vertices() must be
    * called on all machines simultaneously.
    *
    * ### Basic Usage 
    * For instance, if the graph has float vertex data, and float edge data:
    * \code
    *   typedef graphlab::distributed_graph<float, float> graph_type;
    * \endcode
    * 
    * To compute an absolute sum over all the vertex data, we would write
    * a function which reads in each a vertex, and returns the absolute
    * value of the data on the vertex.
    * \code
    * float absolute_vertex_data(graph_type::vertex_type vertex) {
    *   return std::fabs(vertex.data());
    * }
    * \endcode
    * After which calling:
    * \code
    * float sum = graph.map_reduce_vertices<float>(absolute_vertex_data);
    * \endcode
    * will call the <code>absolute_vertex_data()</code> function
    * on each vertex in the graph. <code>absolute_vertex_data()</code>
    * reads the value of the vertex and returns the absolute result.
    * This return values are then summed together and returned. 
    * All machines see the same result.
    *
    * The template argument <code><float></code> is needed to inform
    * the compiler regarding the return type of the mapfunction.
    *
    * ### Relations
    * This function is similar to 
    * graphlab::iengine::map_reduce_vertices()
    * with the difference that this does not take a context
    * and thus cannot influence engine signalling.
    * transform_vertices() can be used to perform a similar
    * but may also make modifications to graph data.
    *
    * \tparam ResultType The output of the map function. Must have
    *                    operator+= defined, and must be \ref Serializable.
    * \tparam VertexMapperType The type of the map function. 
    *                          Not generally needed.
    *                          Can be inferred by the compiler.
    * \param mapfunction The map function to use. Must take 
    *                   a \ref vertex_type, or a reference to a 
    *                   \ref vertex_type as its only argument.
    *                   Returns a ResultType which must be summable
    *                   and \ref Serializable .
    */
    template <typename ResultType, typename MapFunctionType>
    ResultType map_reduce_vertices(MapFunctionType mapfunction) {
      ASSERT_TRUE(finalized);
      rpc.barrier();
      bool global_result_set = false;
      ResultType global_result = ResultType();
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        bool result_set = false;
        ResultType result = ResultType();
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int i = 0; i < (int)local_graph.num_vertices(); ++i) {
          if (lvid2record[i].owner == rpc.procid()) {
            if (!result_set) {
              vertex_type vtx(l_vertex(i));
              result = mapfunction(vtx);
              result_set = true;
            }
            else if (result_set){
              vertex_type vtx(l_vertex(i));
              result += mapfunction(vtx);
            }
          }
        }
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
          if (!global_result_set) {
            global_result = result;
            global_result_set = true;
          }
          else {
            global_result += result;
          }
        }
      }
      conditional_addition_wrapper<ResultType> wrapper(global_result, global_result_set);
      rpc.all_reduce(wrapper);
      return wrapper.value;
    }

   /**
    * \brief Performs a map-reduce operation on each edge in the 
    * graph returning the result.
    * 
    * Given a map function, map_reduce_edges() call the map function on all
    * edges in the graph. The return values are then summed together and the
    * final result returned. The map function should only read data
    * and should not make any modifications. map_reduce_edges() must be
    * called on all machines simultaneously.
    *
    * ### Basic Usage 
    * For instance, if the graph has float vertex data, and float edge data:
    * \code
    *   typedef graphlab::distributed_graph<float, float> graph_type;
    * \endcode
    * 
    * To compute an absolute sum over all the edge data, we would write
    * a function which reads in each a edge, and returns the absolute
    * value of the data on the edge.
    * \code
    * float absolute_edge_data(graph_type::edge_type edge) {
    *   return std::fabs(edge.data());
    * }
    * \endcode
    * After which calling:
    * \code
    * float sum = engine.map_reduce_edges<float>(absolute_edge_data);
    * \endcode
    * will call the <code>absolute_edge_data()</code> function
    * on each edge in the graph. <code>absolute_edge_data()</code>
    * reads the value of the edge and returns the absolute result.
    * This return values are then summed together and returned. 
    * All machines see the same result.
    *
    * The template argument <code><float></code> is needed to inform
    * the compiler regarding the return type of the mapfunction.
    *
    * ### Relations
    * This function similar to 
    * graphlab::distributed_graph::map_reduce_edges()
    * with the difference that this does not take a context
    * and thus cannot influence engine signalling.
    * Finally transform_edges() can be used to perform a similar
    * but may also make modifications to graph data.
    *
    * \tparam ResultType The output of the map function. Must have
    *                    operator+= defined, and must be \ref Serializable.
    * \tparam EdgeMapperType The type of the map function. 
    *                          Not generally needed.
    *                          Can be inferred by the compiler.
    * \param mapfunction The map function to use. Must take 
    *                   a \ref edge_type, or a reference to a 
    *                   \ref edge_type as its only argument.
    *                   Returns a ResultType which must be summable
    *                   and \ref Serializable .
    */
   template <typename ResultType, typename MapFunctionType>
    ResultType map_reduce_edges(MapFunctionType mapfunction) {
      ASSERT_TRUE(finalized);
      rpc.barrier();
      bool global_result_set = false;
      ResultType global_result = ResultType();
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        bool result_set = false;
        ResultType result = ResultType();
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int i = 0; i < (int)local_graph.num_vertices(); ++i) {
          foreach(const local_edge_type& e, l_vertex(i).in_edges()) {
            if (!result_set) {
              edge_type edge(e);
              result = mapfunction(edge);
              result_set = true;
            }
            else if (result_set){
              edge_type edge(e);
              result += mapfunction(edge);
            }
          }
        }
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
         if (!global_result_set) {
            global_result = result;
            global_result_set = true;
          }
          else {
            global_result += result;
          }
        }
      }

      conditional_addition_wrapper<ResultType> wrapper(global_result, global_result_set);
      rpc.all_reduce(wrapper);
      return wrapper.value;
    }

    /**
     * \brief Performs a transformation operation on each vertex in the graph.
     *
     * Given a mapfunction, transform_vertices() calls mapfunction on 
     * every vertex in graph. The map function may make modifications
     * to the data on the vertex. transform_vertices() must be called by all
     * machines simultaneously.
     *
     * ### Basic Usage 
     * For instance, if the graph has integer vertex data, and integer edge
     * data: 
     * \code
     *   typedef graphlab::distributed_graph<size_t, size_t> graph_type;
     * \endcode
     * 
     * To set each vertex value to be the number of out-going edges,
     * we may write the following function:     
     * \code
     * void set_vertex_value(graph_type::vertex_type vertex) {
     *   vertex.data() = vertex.num_out_edges();
     * }
     * \endcode
     *
     * Calling transform_vertices():
     * \code
     *   engine.transform_vertices(set_vertex_value);
     * \endcode
     * will run the <code>set_vertex_value()</code> function
     * on each vertex in the graph, setting its new value. 
     *
     * ### Relations
     * map_reduce_vertices() provide similar signalling functionality, 
     * but should not make modifications to graph data. 
     * graphlab::iengine::transform_vertices() provide
     * the same graph modification capabilities, but with a context
     * and thus can perform signalling.
     *
     * \tparam VertexMapperType The type of the map function. 
     *                          Not generally needed.
     *                          Can be inferred by the compiler.
     * \param mapfunction The map function to use. Must take an
     *                   \ref icontext_type& as its first argument, and
     *                   a \ref vertex_type, or a reference to a 
     *                   \ref vertex_type as its second argument.
     *                   Returns void.
     */ 
    template <typename TransformType>
    void transform_vertices(TransformType transform_functor) {
      ASSERT_TRUE(finalized);
      rpc.barrier();
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (int i = 0; i < (int)local_graph.num_vertices(); ++i) {
        if (lvid2record[i].owner == rpc.procid()) {
          vertex_type vtx(l_vertex(i));
          transform_functor(vtx);
        }
      }
      rpc.barrier();
      synchronize();
    }

    /**
     * \brief Performs a transformation operation on each edge in the graph.
     *
     * Given a mapfunction, transform_edges() calls mapfunction on 
     * every edge in graph. The map function may make modifications
     * to the data on the edge. transform_edges() must be called on
     * all machines simultaneously.
     *
     * ### Basic Usage 
     * For instance, if the graph has integer vertex data, and integer edge
     * data: 
     * \code
     *   typedef graphlab::distributed_graph<size_t, size_t> graph_type;
     * \endcode
     * 
     * To set each edge value to be the number of out-going edges
     * of the target vertex, we may write the following:
     * \code
     * void set_edge_value(graph_type::edge_type edge) {
     *   edge.data() = edge.target().num_out_edges();
     * }
     * \endcode
     *
     * Calling transform_edges():
     * \code
     *   engine.transform_edges(set_edge_value);
     * \endcode
     * will run the <code>set_edge_value()</code> function
     * on each edge in the graph, setting its new value. 
     *
     * ### Relations
     * map_reduce_edges() provide similar signalling functionality, 
     * but should not make modifications to graph data. 
     * graphlab::iengine::transform_edges() provide
     * the same graph modification capabilities, but with a context
     * and thus can perform signalling.
     *
     * \tparam EdgeMapperType The type of the map function. 
     *                          Not generally needed.
     *                          Can be inferred by the compiler.
     * \param mapfunction The map function to use. Must take an
     *                   \ref icontext_type& as its first argument, and
     *                   a \ref edge_type, or a reference to a 
     *                   \ref edge_type as its second argument.
     *                   Returns void.
     */ 
    template <typename TransformType>
    void transform_edges(TransformType transform_functor) {
      ASSERT_TRUE(finalized);
      rpc.barrier();
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (int i = 0; i < (int)local_graph.num_vertices(); ++i) {
        foreach(const local_edge_type& e, l_vertex(i).in_edges()) {
          edge_type edge(e);
          transform_functor(edge);
        }
      }
      rpc.barrier();
    }

    // disable documentation for parallel_for stuff. These are difficult
    // to use properly by the user
    /// \cond GRAPHLAB_INTERNAL
    /**
     * \internal 
     * parallel_for_vertices will partition the set of vertices among the
     * vector of accfunctions. Each accfunction is then executed sequentially
     * on the set of vertices it was assigned.
     *
      * \param accfunction must be a void function which takes a single
      * vertex_type argument. It may be a functor and contain state.
      * The function need not be reentrant as it is only called sequentially
     */
    template <typename VertexFunctorType>
    void parallel_for_vertices(
        std::vector<VertexFunctorType>& accfunction) {
      ASSERT_TRUE(finalized);
      rpc.barrier();
      int numaccfunctions = (int)accfunction.size();
      ASSERT_GE(numaccfunctions, 1);
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (int i = 0; i < (int)accfunction.size(); ++i) {
        for (int j = i;j < (int)local_graph.num_vertices(); j+=numaccfunctions) {
          if (lvid2record[j].owner == rpc.procid()) {
            accfunction[i](vertex_type(l_vertex(j)));
          }
        }
      }
      rpc.barrier();
    }

    
    /**
     * \internal
     * parallel_for_edges will partition the set of edges among the
     * vector of accfunctions. Each accfunction is then executed sequentially
     * on the set of edges it was assigned.
     *
     * \param accfunction must be a void function which takes a single
     * edge_type argument. It may be a functor and contain state.
     * The function need not be reentrant as it is only called sequentially
     */
    template <typename EdgeFunctorType>
    void parallel_for_edges(
        std::vector<EdgeFunctorType>& accfunction) {
      ASSERT_TRUE(finalized);
      rpc.barrier();
      int numaccfunctions = (int)accfunction.size();
      ASSERT_GE(numaccfunctions, 1);
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (int i = 0; i < (int)accfunction.size(); ++i) {
        for (int j = i;j < (int)local_graph.num_vertices(); j+=numaccfunctions) {
          foreach(const local_edge_type& e, l_vertex(j).in_edges()) {
            accfunction[i](edge_type(e));
          }
        }
      }
      rpc.barrier();
    }

    /// \endcond
    
    /// \brief Clears the graph. 
    void clear () { 
      foreach (vertex_record& vrec, lvid2record)
        vrec.clear();
      lvid2record.clear();
      vid2lvid.clear();
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      // read the vertices 
      arc >> nverts 
          >> nedges 
          >> local_own_nverts 
          >> nreplicas
          >> begin_eid
          >> vid2lvid
          >> lvid2record
          >> local_graph;
      finalized = true;
      // check the graph condition
    } // end of load


    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      ASSERT_TRUE(finalized);
      // Write the number of edges and vertices
      arc << nverts 
          << nedges 
          << local_own_nverts 
          << nreplicas 
          << begin_eid
          << vid2lvid
          << lvid2record
          << local_graph;
    } // end of save


    /** \brief Load part of the distributed graph from a path*/
    void load_binary(const std::string& prefix) {
      rpc.full_barrier();
      std::string fname = prefix + tostr(rpc.procid()) + ".bin";

      logstream(LOG_INFO) << "Load graph from " << fname << std::endl;
      if(boost::starts_with(fname, "hdfs://")) {
        graphlab::hdfs hdfs;
        graphlab::hdfs::fstream in_file(hdfs, fname);
        boost::iostreams::filtering_stream<boost::iostreams::input> fin;
        fin.push(boost::iostreams::gzip_decompressor());
        fin.push(in_file);
        
        if(!fin.good()) {
          logstream(LOG_FATAL) << "Error opening file: " << fname << std::endl;
          exit(-1);
        }
        iarchive iarc(fin);
        iarc >> *this;
        fin.pop();
        fin.pop();
        in_file.close();
      } else {
        std::ifstream in_file(fname.c_str(),
                              std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_stream<boost::iostreams::input> fin;
        fin.push(boost::iostreams::gzip_decompressor());
        fin.push(in_file);
        iarchive iarc(fin);
        iarc >> *this;
        fin.pop();
        fin.pop();
        in_file.close();
      }
      logstream(LOG_INFO) << "Finish loading graph from " << fname << std::endl;
      rpc.full_barrier();
    } // end of load


    /** \brief Load part of the distributed graph from a path*/
    void save_binary(const std::string& prefix) {
      rpc.full_barrier();
      timer savetime;  savetime.start();
      std::string fname = prefix + tostr(rpc.procid()) + ".bin";
      logstream(LOG_INFO) << "Save graph to " << fname << std::endl;
      if(boost::starts_with(fname, "hdfs://")) {
        graphlab::hdfs hdfs;
        graphlab::hdfs::fstream out_file(hdfs, fname, true);
        boost::iostreams::filtering_stream<boost::iostreams::output> fout;
        fout.push(boost::iostreams::gzip_compressor());        
        fout.push(out_file);
        if (!fout.good()) {
          logstream(LOG_FATAL) << "Error opening file: " << fname << std::endl;
          exit(-1);
        }
        oarchive oarc(fout);
        oarc << *this;
        fout.pop();
        fout.pop();
        out_file.close();
      } else {
        std::ofstream out_file(fname.c_str(),
                               std::ios_base::out | std::ios_base::binary);
        if (!out_file.good()) {
          logstream(LOG_FATAL) << "Error opening file: " << fname << std::endl;
          exit(-1);
        }
        boost::iostreams::filtering_stream<boost::iostreams::output> fout;
        fout.push(boost::iostreams::gzip_compressor());        
        fout.push(out_file);
        oarchive oarc(fout);
        oarc << *this;
        fout.pop();
        fout.pop();
        out_file.close();
      }
      logstream(LOG_INFO) << "Finish saving graph to " << fname << std::endl
                          << "Finished saving binary graph: " 
                          << savetime.current_time() << std::endl;
      rpc.full_barrier();
    } // end of save





    template<typename Writer>
    void save_to_posixfs(const std::string& prefix, Writer writer,
                         bool gzip = true,
                         bool save_vertex = true,
                         bool save_edge = true,
                         size_t files_per_machine = 4) {
      typedef boost::function<void(vertex_type)> vertex_function_type;
      typedef boost::function<void(edge_type)> edge_function_type;
      typedef std::ofstream base_fstream_type;
      typedef boost::iostreams::filtering_stream<boost::iostreams::output>
        boost_fstream_type;
      rpc.full_barrier();
      // figure out the filenames
      std::vector<std::string> graph_files;
      std::vector<base_fstream_type*> outstreams;
      std::vector<boost_fstream_type*> booststreams;
      graph_files.resize(files_per_machine);
      for(size_t i = 0; i < files_per_machine; ++i) {
        graph_files[i] = prefix + "." + tostr(1 + i + rpc.procid() * files_per_machine)
          + "_of_" + tostr(rpc.numprocs() * files_per_machine);
        if (gzip) graph_files[i] += ".gz";
      }

      // create the vector of callbacks
      std::vector<vertex_function_type> vertex_callbacks(graph_files.size());
      std::vector<edge_function_type> edge_callbacks(graph_files.size());
    
      for(size_t i = 0; i < graph_files.size(); ++i) {
        logstream(LOG_INFO) << "Saving to file: " << graph_files[i] << std::endl;
        // open the stream
        base_fstream_type* out_file = 
          new base_fstream_type(graph_files[i].c_str(),
                                std::ios_base::out | std::ios_base::binary);
        // attach gzip if the file is gzip
        boost_fstream_type* fout = new boost_fstream_type;
        // Using gzip filter
        if (gzip) fout->push(boost::iostreams::gzip_compressor());
        fout->push(*out_file);

        outstreams.push_back(out_file);
        booststreams.push_back(fout);
        // construct the callback for the parallel for
        typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;
        vertex_callbacks[i] = 
          boost::bind(&graph_type::template save_vertex_to_stream<boost_fstream_type, Writer>,
                      this, _1, boost::ref(*fout), boost::ref(writer));
        edge_callbacks[i] =
          boost::bind(&graph_type::template save_edge_to_stream<boost_fstream_type, Writer>,
                      this, _1, boost::ref(*fout), boost::ref(writer));
      }

      if (save_vertex) parallel_for_vertices(vertex_callbacks);
      if (save_edge) parallel_for_edges(edge_callbacks);

      // cleanup
      for(size_t i = 0; i < graph_files.size(); ++i) {    
        booststreams[i]->pop();
        if (gzip) booststreams[i]->pop();
        delete booststreams[i];
        delete outstreams[i];
      }
      vertex_callbacks.clear();
      edge_callbacks.clear();
      outstreams.clear();
      booststreams.clear();
      rpc.full_barrier();
    } // end of save to posixfs




    template<typename Writer>
    void save_to_hdfs(const std::string& prefix, Writer writer,
                      bool gzip = true,
                      bool save_vertex = true,
                      bool save_edge = true,
                      size_t files_per_machine = 4) {
      typedef boost::function<void(vertex_type)> vertex_function_type;
      typedef boost::function<void(edge_type)> edge_function_type;
      typedef graphlab::hdfs::fstream base_fstream_type;
      typedef boost::iostreams::filtering_stream<boost::iostreams::output>
        boost_fstream_type;
      rpc.full_barrier();
      // figure out the filenames
      std::vector<std::string> graph_files;
      std::vector<base_fstream_type*> outstreams;
      std::vector<boost_fstream_type*> booststreams;
      graph_files.resize(files_per_machine);
      for(size_t i = 0; i < files_per_machine; ++i) {
        graph_files[i] = prefix + "." + tostr(1 + i + rpc.procid() * files_per_machine)
          + "_of_" + tostr(rpc.numprocs() * files_per_machine);
        if (gzip) graph_files[i] += ".gz";
      }

      ASSERT_TRUE(hdfs::has_hadoop());
      hdfs& hdfs = hdfs::get_hdfs();
    
      // create the vector of callbacks

      std::vector<vertex_function_type> vertex_callbacks(graph_files.size());
      std::vector<edge_function_type> edge_callbacks(graph_files.size());

      for(size_t i = 0; i < graph_files.size(); ++i) {
        logstream(LOG_INFO) << "Saving to file: " << graph_files[i] << std::endl;
        // open the stream
        base_fstream_type* out_file = new base_fstream_type(hdfs,
                                                            graph_files[i],
                                                            true);
        // attach gzip if the file is gzip
        boost_fstream_type* fout = new boost_fstream_type;
        // Using gzip filter
        if (gzip) fout->push(boost::iostreams::gzip_compressor());
        fout->push(*out_file);

        outstreams.push_back(out_file);
        booststreams.push_back(fout);
        // construct the callback for the parallel for
        typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;
        vertex_callbacks[i] = 
          boost::bind(&graph_type::template save_vertex_to_stream<boost_fstream_type, Writer>,
                      this, _1, boost::ref(*fout), writer);
        edge_callbacks[i] =
          boost::bind(&graph_type::template save_edge_to_stream<boost_fstream_type, Writer>,
                      this, _1, boost::ref(*fout), writer);
      }

      if (save_vertex) parallel_for_vertices(vertex_callbacks);
      if (save_edge) parallel_for_edges(edge_callbacks);

      // cleanup
      for(size_t i = 0; i < graph_files.size(); ++i) {
        booststreams[i]->pop();
        if (gzip) booststreams[i]->pop();
        delete booststreams[i];
        delete outstreams[i];
      }
      vertex_callbacks.clear();
      edge_callbacks.clear();
      outstreams.clear();
      booststreams.clear();
      rpc.full_barrier();
    } // end of save to hdfs



    template<typename Writer>
    void save(const std::string& prefix, Writer writer,
              bool gzip = true, bool save_vertex = true, bool save_edge = true,
              size_t files_per_machine = 4) {
      if(boost::starts_with(prefix, "hdfs://")) {
        save_to_hdfs(prefix, writer, gzip, save_vertex, save_edge, files_per_machine);
      } else {
        save_to_posixfs(prefix, writer, gzip, save_vertex, save_edge, files_per_machine);
      }
    } // end of save




    void save_format(const std::string& prefix, const std::string& format,
                        bool gzip = true, size_t files_per_machine = 4) {
      if (format == "snap" || format == "tsv") {
        save(prefix, builtin_parsers::tsv_writer<distributed_graph>(),
             gzip, false, true, files_per_machine);
      } else if (format == "graphjrl") {
         save(prefix, builtin_parsers::graphjrl_writer<distributed_graph>(),
             gzip, true, true, files_per_machine);
      } else if (format == "bin") {
         save_binary(prefix);
      } else {
        logstream(LOG_ERROR)
          << "Unrecognized Format \"" << format << "\"!" << std::endl;
        return;
      }
    } // end of save structure




    /**
       Load a graph from a collection of files in stored in path using
       the user defined line parser.
     */
    void load_from_posixfs(const std::string& original_path, line_parser_type line_parser) {
      std::string directory_name; std::string prefix;
      boost::filesystem::path path(original_path);
      if (boost::filesystem::is_directory(path)) {
        // if this is a directory
        // force a "/" at the end of the path
        // make sure to check that the path is non-empty. (you do not
        // want to make the empty path "" the root path "/" )
        directory_name = path.native();
      }
      else {
        directory_name = path.parent_path().native();
        prefix = path.filename().native();
        directory_name = (directory_name.empty() ? "." : directory_name);
      }
      std::vector<std::string> graph_files;
      fs_util::list_files_with_prefix(directory_name, prefix, graph_files);
      if (graph_files.size() == 0) {
        logstream(LOG_WARNING) << "No files found matching " << original_path << std::endl;
      }
      for(size_t i = 0; i < graph_files.size(); ++i) {
        if (i % rpc.numprocs() == rpc.procid()) {
          logstream(LOG_EMPH) << "Loading graph from file: " << graph_files[i] << std::endl;
          // is it a gzip file ?
          const bool gzip = boost::ends_with(graph_files[i], ".gz");
          // open the stream
          std::ifstream in_file(graph_files[i].c_str(), 
                                std::ios_base::in | std::ios_base::binary);
          // attach gzip if the file is gzip
          boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
          // Using gzip filter
          if (gzip) fin.push(boost::iostreams::gzip_decompressor());
          fin.push(in_file);
          const bool success = load_from_stream(graph_files[i], fin, line_parser);
          if(!success) {
            logstream(LOG_FATAL) 
              << "Error parsing file: " << graph_files[i] << std::endl;
          }
          fin.pop();
          if (gzip) fin.pop();
        }
      }
      rpc.full_barrier();
    } // end of load from posixfs


    /**
       Load a graph from a collection of files in stored in path using
       the user defined line parser.
    */
    void load_from_hdfs(const std::string& original_path, line_parser_type line_parser) {
      // force a "/" at the end of the path
      // make sure to check that the path is non-empty. (you do not
      // want to make the empty path "" the root path "/" )
      std::string path = original_path;
      if (path.length() > 0 && path[path.length() - 1] != '/') path = path + "/";
      ASSERT_TRUE(hdfs::has_hadoop());
      hdfs& hdfs = hdfs::get_hdfs();    
      std::vector<std::string> graph_files;
      graph_files = hdfs.list_files(path);
      if (graph_files.size() == 0) {
        logstream(LOG_WARNING) << "No files found matching " << original_path << std::endl;
      }
      for(size_t i = 0; i < graph_files.size(); ++i) {
        if (i % rpc.numprocs() == rpc.procid()) {
          logstream(LOG_EMPH) << "Loading graph from file: " << graph_files[i] << std::endl;
          // is it a gzip file ?
          const bool gzip = boost::ends_with(graph_files[i], ".gz");
          // open the stream
          graphlab::hdfs::fstream in_file(hdfs, graph_files[i]);
          boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
          if(gzip) fin.push(boost::iostreams::gzip_decompressor());
          fin.push(in_file);      
          const bool success = load_from_stream(graph_files[i], fin, line_parser);
          if(!success) {
            logstream(LOG_FATAL) 
              << "Error parsing file: " << graph_files[i] << std::endl;
          }
          fin.pop();
          if (gzip) fin.pop();
        }
      }
      rpc.full_barrier();
    } // end of load from hdfs


    /**
       Load a the graph from a given path using the user defined line parser.
     */  
    void load(const std::string& path, line_parser_type line_parser) {
      rpc.full_barrier();
      if(boost::starts_with(path, "hdfs://")) {
        load_from_hdfs(path, line_parser);
      } else {
        load_from_posixfs(path, line_parser);
      }
      rpc.full_barrier();
    } // end of load
  

    // Synthetic Generators ===================================================>
    void load_synthetic_powerlaw(size_t nverts, bool in_degree = false,
                                 double alpha = 2.1, size_t truncate = (size_t)(-1)) {
      rpc.full_barrier(); 
      std::vector<double> prob(std::min(nverts, truncate), 0);
      logstream(LOG_INFO) << "constructing pdf" << std::endl;
      for(size_t i = 0; i < prob.size(); ++i)
        prob[i] = std::pow(double(i+1), -alpha);
      logstream(LOG_INFO) << "constructing cdf" << std::endl;
      random::pdf2cdf(prob);
      logstream(LOG_INFO) << "Building graph" << std::endl;
      size_t target_index = rpc.procid();
      size_t addedvtx = 0;

      for(size_t source = rpc.procid(); source < nverts;
          source += rpc.numprocs()) {
        const size_t out_degree = random::sample(prob) + 1;
        for(size_t i = 0; i < out_degree; ++i) {
          target_index = (target_index + 2654435761)  % nverts;
          while (source == target_index) {
            target_index = (target_index + 2654435761)  % nverts;
          }
          if(in_degree) add_edge(target_index, source);
          else add_edge(source, target_index);
        }
        ++addedvtx;
        if (addedvtx % 10000000 == 0) {
          logstream(LOG_EMPH) << addedvtx << " inserted\n";
        }
      }
      rpc.full_barrier(); 
    } // end of load random powerlaw


    /**
       load a graph with a standard format
       \todo: finish documentation of formats
     */
    void load_format(const std::string& path, const std::string& format) {
      line_parser_type line_parser;
      if (format == "snap") {
        line_parser = builtin_parsers::snap_parser<distributed_graph>;
      } else if (format == "adj") {
        line_parser = builtin_parsers::adj_parser<distributed_graph>;
      } else if (format == "tsv") {
        line_parser = builtin_parsers::tsv_parser<distributed_graph>;
      } else if (format == "graphjrl") {
        line_parser = builtin_parsers::graphjrl_parser<distributed_graph>;
      } else if (format == "bin") {
         load_binary(path);
      } else {
        logstream(LOG_ERROR)
          << "Unrecognized Format \"" << format << "\"!" << std::endl;
        return;
      }
      load(path, line_parser);
    } // end of load












/****************************************************************************
 *                       Internal Functions                                 *
 *                     ----------------------                               *
 * These functions functions and types provide internal access to the       *
 * underlying graph representation. They should not be used unless you      *
 * *really* know what you are doing.                                        *
 ****************************************************************************/


    /**
     * \internal
     * The vertex record stores information associated with each
     * vertex on this proc
     */
    struct vertex_record {
      /// The official owning processor for this vertex
      procid_t owner; 
      /// The local vid of this vertex on this proc
      vertex_id_type gvid;
      /// The number of in edges
      vertex_id_type num_in_edges, num_out_edges;
      /** The set of proc that mirror this vertex.  The owner should
          NOT be in this set.*/
      mirror_type _mirrors;
      vertex_record() : 
        owner(-1), gvid(-1), num_in_edges(0), num_out_edges(0) { }
      vertex_record(const vertex_id_type& vid) : 
        owner(-1), gvid(vid), num_in_edges(0), num_out_edges(0) { }
      procid_t get_owner () const { return owner; }
      const mirror_type& mirrors() const { return _mirrors; }
      size_t num_mirrors() const { return _mirrors.popcount(); }

      void clear() {
        _mirrors.clear();
      }

      void load(iarchive& arc) {
        clear();
        arc >> owner
            >> gvid
            >> num_in_edges
            >> num_out_edges
            >> _mirrors;
      }

      void save(oarchive& arc) const {
        arc << owner
            << gvid
            << num_in_edges
            << num_out_edges
            << _mirrors;
      } // end of save
    }; // end of vertex_record




    /** \internal
     * \brief converts a local vertex ID to a local vertex object
     */
    local_vertex_type l_vertex(lvid_type vid) {
      return local_vertex_type(*this, vid);
    }
    
    /** \internal
     *\brief Get the Total number of vertex replicas in the graph */
    size_t num_replicas() const { return nreplicas; }

    /** \internal
     *\brief Get the number of vertices local to this proc */
    size_t num_local_vertices() const { return local_graph.num_vertices(); }

    /** \internal
     *\brief Get the number of edges local to this proc */
    size_t num_local_edges() const { return local_graph.num_edges(); }

    /** \internal
     *\brief Get the number of vertices owned by this proc */
    size_t num_local_own_vertices() const { return local_own_nverts; }

    /** \internal
     *\brief Convert a global vid to a local vid */
    lvid_type local_vid (const vertex_id_type vid) const {
      // typename boost::unordered_map<vertex_id_type, lvid_type>::
      //   const_iterator iter = vid2lvid.find(vid);
      typename cuckoo_map_type::const_iterator iter = vid2lvid.find(vid);
      ASSERT_TRUE(iter != vid2lvid.end());
      return iter->second;
    } // end of local_vertex_id

    /** \internal
     *\brief Convert a local vid to a global vid */
    vertex_id_type global_vid(const lvid_type lvid) const { 
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid].gvid;
    } // end of global_vertex_id



    /**
     * \internal
     * \brief Returns an edge list of all in edges of a local vertex ID
     *        on the local graph
     * 
     * Equivalent to l_vertex(lvid).in_edges()
     */    
    local_edge_list_type l_in_edges(const lvid_type lvid) {
      return local_edge_list_type(*this, local_graph.in_edges(lvid));
    }

    /**
     * \internal
     * \brief Returns the number of in edges of a local vertex ID
     *        on the local graph
     * 
     * Equivalent to l_vertex(lvid).num in_edges()
     */    
    size_t l_num_in_edges(const lvid_type lvid) const { 
      return local_graph.num_in_edges(lvid);
    }

    /**
     * \internal
     * \brief Returns an edge list of all out edges of a local vertex ID
     *        on the local graph
     * 
     * Equivalent to l_vertex(lvid).out_edges()
     */    
    local_edge_list_type l_out_edges(const lvid_type lvid) {
      return local_edge_list_type(*this, local_graph.out_edges(lvid));
    }

    /**
     * \internal
     * \brief Returns the number of out edges of a local vertex ID
     *        on the local graph
     * 
     * Equivalent to l_vertex(lvid).num out_edges()
     */    
    size_t l_num_out_edges(const lvid_type lvid) const {
      return local_graph.num_out_edges(lvid);
    }
    
    procid_t procid() const {
      return rpc.procid();
    }


    procid_t numprocs() const {
      return rpc.numprocs();
    }

    distributed_control& dc() {
      return rpc.dc();
    }



    /** \internal
     * \brief Returns the internal vertex record of a given global vertex ID
     */
    const vertex_record& get_vertex_record(vertex_id_type vid) const {
      // typename boost::unordered_map<vertex_id_type, lvid_type>::
      //   const_iterator iter = vid2lvid.find(vid);
      typename cuckoo_map_type::const_iterator iter = vid2lvid.find(vid);
      ASSERT_TRUE(iter != vid2lvid.end());
      return lvid2record[iter->second];
    }

    /** \internal
     * \brief Returns the internal vertex record of a given local vertex ID
     */
    vertex_record& l_get_vertex_record(lvid_type lvid) {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid];
    }

    /** \internal
     * \brief Returns the internal vertex record of a given local vertex ID
     */
    const vertex_record& l_get_vertex_record(lvid_type lvid) const {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid];
    }

    /** \internal
     * \brief Returns true if the provided global vertex ID is a 
     *        master vertex on this machine and false otherwise.
     */
    bool is_master(vertex_id_type vid) const {
      typename cuckoo_map_type::const_iterator iter = vid2lvid.find(vid);
      return (iter != vid2lvid.end()) && l_is_master(iter->second);
    }
    /** \internal
     * \brief Returns true if the provided local vertex ID is a master vertex.
     *        Returns false otherwise.
     */
    bool l_is_master(lvid_type lvid) const {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid].owner == rpc.procid();
    }

    /** \internal
     * \brief Returns the master procid for vertex lvid.
     */
    procid_t l_master(lvid_type lvid) const {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid].owner;
    }


    /** \internal
     *  \brief Returns a reference to the internal graph representation
     */
    local_graph_type& get_local_graph() {
      return local_graph;
    }

    /** \internal
     *  \brief Returns a const reference to the internal graph representation
     */
    const local_graph_type& get_local_graph() const {
      return local_graph;
    }




    /** \internal
     * This function synchronizes the master vertex data with all the mirrors.
     * This function must be called simultaneously by all machines
     */
    void synchronize() {
      typedef std::pair<vertex_id_type, vertex_data_type> pair_type;
      buffered_exchange<pair_type> vertex_exchange(rpc.dc());
      typename buffered_exchange<pair_type>::buffer_type recv_buffer;
      procid_t sending_proc;
      // Loop over all the local vertex records
      for(lvid_type lvid = 0; lvid < lvid2record.size(); ++lvid) {
        const vertex_record& record = lvid2record[lvid];
        // if this machine is the owner of a record then send the
        // vertex data to all mirrors
        if(record.owner == rpc.procid()) {
          foreach(uint32_t proc, record.mirrors()) {
            const pair_type pair(record.gvid, local_graph.vertex_data(lvid));
            vertex_exchange.send(proc, pair);
          }
        }
        // Receive any vertex data and update local mirrors
        while(vertex_exchange.recv(sending_proc, recv_buffer)) {
          foreach(const pair_type& pair, recv_buffer)  {
            vertex(pair.first).data() = pair.second;
          }
          recv_buffer.clear();
        }
      }
      vertex_exchange.flush();
      while(vertex_exchange.recv(sending_proc, recv_buffer)) {
        foreach(const pair_type& pair, recv_buffer) {
          vertex(pair.first).data() = pair.second;
        }
        recv_buffer.clear();
      }
      ASSERT_TRUE(vertex_exchange.empty());
    } // end of synchronize






    /** \internal
     *  vertex type while provides access to local graph vertices.
     */
    struct local_vertex_type {
      distributed_graph& graph_ref;
      lvid_type lvid;

      local_vertex_type(distributed_graph& graph_ref, lvid_type lvid):
            graph_ref(graph_ref), lvid(lvid) { }

      /// \brief Can be casted from local_vertex_type using an explicit cast
      explicit local_vertex_type(vertex_type v) :graph_ref(v.graph_ref),lvid(v.lvid) { }
      /// \brief Can be casted to vertex_type using an explicit cast
      operator vertex_type() const {
        return vertex_type(graph_ref, lvid);
      }

      bool operator==(local_vertex_type& v) const {
        return lvid == v.lvid;
      }
      
      /// \brief Returns a reference to the data on the local vertex
      const vertex_data_type& data() const {
        return graph_ref.get_local_graph().vertex_data(lvid);
      }

      /// \brief Returns a reference to the data on the local vertex
      vertex_data_type& data() {
        return graph_ref.get_local_graph().vertex_data(lvid);
      }

      /** \brief Returns the number of in edges on the 
       *         local graph of this local vertex
       */
      size_t num_in_edges() const {
        return graph_ref.get_local_graph().num_in_edges(lvid);
      }

      /** \brief Returns the number of in edges on the 
       *         local graph of this local vertex
       */
      size_t num_out_edges() const {
        return graph_ref.get_local_graph().num_out_edges(lvid);
      }

      /// \brief Returns the local ID of this local vertex
      lvid_type id() const {
        return lvid;
      }

      /// \brief Returns the global ID of this local vertex
      vertex_id_type global_id() const {
        return graph_ref.global_vid(lvid);
      }

      /** \brief Returns a list of all in edges on the 
       *         local graph of this local vertex
       */
      local_edge_list_type in_edges() {
        return graph_ref.l_in_edges(lvid);
      }

      /** \brief Returns a list of all out edges on the 
       *         local graph of this local vertex
       */
      local_edge_list_type out_edges() {
        return graph_ref.l_out_edges(lvid);
      }

      /** \brief Returns the owner of this local vertex
       */
      procid_t owner() const {
        return graph_ref.l_get_vertex_record(lvid).owner;
      }

      /** \brief Returns the number of in_edges of this vertex
       *         on the global graph
       */
      size_t global_num_in_edges() const {
        return graph_ref.l_get_vertex_record(lvid).num_in_edges;
      }


      /** \brief Returns the number of out_edges of this vertex
       *         on the global graph
       */
      size_t global_num_out_edges() const {
        return graph_ref.l_get_vertex_record(lvid).num_out_edges;
      }


      /** \brief Returns the set of mirrors of this vertex
       */
      const mirror_type& mirrors() const {
        return graph_ref.l_get_vertex_record(lvid)._mirrors;
      }

      size_t num_mirrors() const {
        return graph_ref.l_get_vertex_record(lvid).num_mirrors();
      }

      /** \brief Returns the vertex record of this
       *         this local vertex
       */
      vertex_record& get_vertex_record() {
        return graph_ref.l_get_vertex_record(lvid);
      }
    };

    
    /** \internal
     *  edge type which provides access to local graph edges */
    class local_edge_type {
    private:
      distributed_graph& graph_ref;
      typename local_graph_type::edge_type e;
    public:
      local_edge_type(distributed_graph& graph_ref,
                      typename local_graph_type::edge_type e):
                                                graph_ref(graph_ref), e(e) { }
                      
      /// \brief Can be converted from edge_type via an explicit cast
      explicit local_edge_type(edge_type ge) :graph_ref(ge.graph_ref),e(ge.e) { }

      /// \brief Can be casted to edge_type using an explicit cast
      operator edge_type() const {
        return edge_type(graph_ref, e);
      }

      /// \brief Returns the source local vertex of the edge
      local_vertex_type source() { return local_vertex_type(graph_ref, e.source().id()); }
      
      /// \brief Returns the target local vertex of the edge
      local_vertex_type target() { return local_vertex_type(graph_ref, e.target().id()); }
      
      
      
      /// \brief Returns a constant reference to the data on the vertex
      const edge_data_type& data() const { return e.data(); }
      
      /// \brief Returns a reference to the data on the vertex
      edge_data_type& data() { return e.data(); }
      
      /// \brief Returns the internal ID of this edge
      edge_id_type id() const { return e.id(); }
    }; 

    /** \internal
     * \brief A functor which converts local_graph_type::edge_type to
     *        local_edge_type 
     */
    struct make_local_edge_type_functor {
      typedef typename local_graph_type::edge_type argument_type;
      typedef local_edge_type result_type;
      distributed_graph& graph_ref;
      make_local_edge_type_functor(distributed_graph& graph_ref):
                                                      graph_ref(graph_ref) { }
      result_type operator() (const argument_type et) const {
        return local_edge_type(graph_ref, et);
      }
    };
    

    /** \internal
     * \brief A list of edges. Used by l_in_edges() and l_out_edges() 
     */
    struct local_edge_list_type {
      make_local_edge_type_functor me_functor;
      typename local_graph_type::edge_list_type elist;
      
      typedef boost::transform_iterator<make_local_edge_type_functor,
                                      typename local_graph_type::edge_list_type::iterator> iterator;
      typedef iterator const_iterator;
      
      local_edge_list_type(distributed_graph& graph_ref,
                           typename local_graph_type::edge_list_type elist) :
                          me_functor(graph_ref), elist(elist) { }
      /// \brief Returns the number of edges in the list
      size_t size() const { return elist.size(); }
      
      /// \brief Random access to the list elements
      local_edge_type operator[](size_t i) const { return me_functor(elist[i]); }
      
      /** \brief Returns an iterator to the beginning of the list. 
       * 
       * Returns an iterator to the beginning of the list. \see end()
       * The iterator_type is local_edge_list_type::iterator.
       * 
       * \code
       * local_edge_list_type::iterator iter = elist.begin();
       * while(iter != elist.end()) {
       *   ... [do stuff] ...
       *   ++iter;
       * }
       * \endcode
       * 
      */
      iterator begin() const { return
          boost::make_transform_iterator(elist.begin(), me_functor); }
          
      /** \brief Returns an iterator to the end of the list. 
       * 
       * Returns an iterator to the end of the list. \see begin()
       * The iterator_type is local_edge_list_type::iterator.
       * 
       * \code
       * local_edge_list_type::iterator iter = elist.begin();
       * while(iter != elist.end()) {
       *   ... [do stuff] ...
       *   ++iter;
       * }
       * \endcode
       * 
      */
      iterator end() const { return
          boost::make_transform_iterator(elist.end(), me_functor); }
          
      /// \brief Returns true if the list is empty
      bool empty() const { return elist.empty(); }
    }; 
    

  private:
      
    // PRIVATE DATA MEMBERS ===================================================> 
    /** The rpc interface for this class */
    mutable dc_dist_object<distributed_graph> rpc;

    bool finalized;

    /** The local graph data */
    local_graph_type local_graph;
    
    /** The map from global vertex ids to vertex records */
    std::vector<vertex_record>  lvid2record;
    
    // boost::unordered_map<vertex_id_type, lvid_type> vid2lvid;
    /** The map from global vertex ids back to local vertex ids */
    typedef cuckoo_map_pow2<vertex_id_type, lvid_type, 3, uint32_t> cuckoo_map_type;
    cuckoo_map_type vid2lvid;

        
    /** The global number of vertices and edges */
    size_t nverts, nedges;

    /** The number of vertices owned by this proc */
    size_t local_own_nverts;

    /** The global number of vertex replica */
    size_t nreplicas;

    /** The beginning edge id for this machine */
    size_t begin_eid;

    /** pointer to the distributed ingress object*/
    idistributed_ingress<VertexData, EdgeData>* ingress_ptr;



    void set_ingress_method(const std::string& method,
        size_t bufsize = 50000, bool usehash = false, bool userecent = false) {

      if(ingress_ptr != NULL) { delete ingress_ptr; ingress_ptr = NULL; }
      if (method == "batch") {
        logstream(LOG_EMPH) << "Use batch ingress, bufsize: " << bufsize
          << ", usehash: " << usehash << ", userecent" << userecent << std::endl;
        ingress_ptr = new distributed_batch_ingress<VertexData, EdgeData>(rpc.dc(), *this,
                                                         bufsize, usehash, userecent);
      } else if (method == "oblivious") {
        logstream(LOG_EMPH) << "Use oblivious ingress, usehash: " << usehash
          << ", userecent: " << userecent << std::endl;
        ingress_ptr = new distributed_oblivious_ingress<VertexData, EdgeData>(rpc.dc(), *this,
                                                             usehash, userecent);
      } else if (method == "identity") {
        logstream(LOG_EMPH) << "Use identity ingress" << std::endl;
        ingress_ptr = new distributed_identity_ingress<VertexData, EdgeData>(rpc.dc(), *this);
      } else {
        ingress_ptr = new distributed_random_ingress<VertexData, EdgeData>(rpc.dc(), *this);
      }
    } // end of set ingress method


    /**
       \internal
       This internal function is used to load a single line from an input stream
     */
    template<typename Fstream>
    bool load_from_stream(const std::string& filename, Fstream& fin, 
                          line_parser_type& line_parser) {
      size_t linecount = 0;
      timer ti; ti.start();
      while(fin.good() && !fin.eof()) {
        std::string line;
        std::getline(fin, line);
        if(!fin.good()) break;
        const bool success = line_parser(*this, filename, line);
        if (!success) {
          logstream(LOG_WARNING) 
            << "Error parsing line " << linecount << " in "
            << filename << ": " << std::endl
            << "\t\"" << line << "\"" << std::endl;  
          return false;
        }
        ++linecount;      
        if (ti.current_time() > 5.0) {
          logstream(LOG_INFO) << linecount << " Lines read" << std::endl;
          ti.start();
        }
      }
      return true;
    } // end of load from stream


    template<typename Fstream, typename Writer>
    void save_vertex_to_stream(vertex_type& vertex, Fstream& fout, Writer writer) {
      fout << writer.save_vertex(vertex);
    } // end of save_vertex_to_stream


    template<typename Fstream, typename Writer>
    void save_edge_to_stream(edge_type& edge, Fstream& fout, Writer writer) {
      std::string ret = writer.save_edge(edge);
      fout << ret;
    } // end of save_edge_to_stream
  



  }; // End of graph
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

