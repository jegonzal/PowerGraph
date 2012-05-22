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

#include <omp.h>
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

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/util/random.hpp>

#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/graph/ingress/idistributed_ingress.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/ingress/distributed_batch_ingress2.hpp>
#include <graphlab/graph/ingress/distributed_oblivious_ingress.hpp>
#include <graphlab/graph/ingress/distributed_random_ingress.hpp>
#include <graphlab/graph/ingress/distributed_identity_ingress.hpp>
#include <graphlab/util/hdfs.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <graphlab/util/cuckoo_map_pow2.hpp>



#include <graphlab/macros_def.hpp>
namespace graphlab { 


  // CLASS GRAPH ==============================================================>  
  template<typename VertexData, typename EdgeData>
  class distributed_graph {
  public:

    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;

    enum SizeType {MAX_MACHINES = 128};
    typedef fixed_dense_bitset<MAX_MACHINES> mirror_type;

    /// The type of the local graph used to store the graph data 
    // typedef graphlab::graph<VertexData, EdgeData> local_graph_type;
    typedef graphlab::graph<VertexData, EdgeData> local_graph_type;

    typedef idistributed_ingress<VertexData, EdgeData> 
    idistributed_ingress_type;

    friend class distributed_ingress_base<VertexData, EdgeData>;

    typedef distributed_random_ingress<VertexData, EdgeData>
        distributed_random_ingress_type;
    friend class distributed_random_ingress<VertexData, EdgeData>;

    typedef distributed_identity_ingress<VertexData, EdgeData>
        distributed_identity_ingress_type;
    friend class distributed_identity_ingress<VertexData, EdgeData>;


    typedef distributed_batch_ingress<VertexData, EdgeData>
        distributed_batch_ingress_type;
    friend class distributed_batch_ingress<VertexData, EdgeData>;

    typedef distributed_oblivious_ingress<VertexData, EdgeData>
        distributed_oblivious_ingress_type;
    friend class distributed_oblivious_ingress<VertexData, EdgeData>;


    /** 
     * The type of the local vertex id and local edge id.
     * While this is the same as the
     * vertex_id_type giving it a separate name will make method calls
     * clearer.
     */
    typedef typename local_graph_type::vertex_id_type  lvid_type;
    typedef typename local_graph_type::edge_id_type    leid_type;

    struct vertex_type;
    struct local_edge_type;
    struct local_edge_list_type;
    typedef bool edge_list_type;  // map it to an impossible type
    typedef bool edge_type;  // map it to an impossible type
    
    struct vertex_type {
      distributed_graph& g;
      vertex_id_type lvid, gvid;
      vertex_type(distributed_graph& g, vertex_id_type lvid):
            g(g), lvid(lvid), gvid(g.global_vid(lvid)) { }

      const vertex_data_type data() const {
        return g.get_local_graph().vertex_data(lvid);
      }

      vertex_data_type data() {
        return g.get_local_graph().vertex_data(lvid);
      }

      size_t num_in_edges() const {
        return g.l_get_vertex_record(lvid).num_in_edges;
      }

      size_t num_out_edges() const {
        return g.l_get_vertex_record(lvid).num_out_edges;
      }
      
      vertex_id_type id() const {
        return gvid;
      }

      vertex_id_type l_id() const {
        return lvid;
      }

      edge_list_type in_edges() __attribute__ ((noreturn)) {
        ASSERT_TRUE(false);
      }

      edge_list_type out_edges() __attribute__ ((noreturn)) {
        ASSERT_TRUE(false);
      }

      local_edge_list_type l_in_edges() {
        return g.l_in_edges(lvid);
      }

      local_edge_list_type l_out_edges() {
        return g.l_out_edges(lvid);
      }
    };

    typedef vertex_type local_vertex_type;
    
    /** This class represents an edge with source() and target()*/
    class local_edge_type {
    private:
      distributed_graph& g;
      typename local_graph_type::edge_type e;
    public:
      local_edge_type(distributed_graph& g,
                      typename local_graph_type::edge_type e): g(g), e(e) { }

      local_vertex_type source() { return local_vertex_type(g, e.source().id()); }
      local_vertex_type target() { return local_vertex_type(g, e.target().id()); }
      
      edge_data_type data() { return e.data(); }
      const edge_data_type data() const { return e.data(); }
      procid_t owner() const { return g.rmi.procid(); }
      edge_id_type id() const { return e.id(); }
    }; 


    struct make_local_edge_type_functor {
      typedef typename local_graph_type::edge_type argument_type;
      typedef local_edge_type result_type;
      distributed_graph& g;
      make_local_edge_type_functor(distributed_graph& g):g(g) { }
      result_type operator() (const argument_type et) const {
        return local_edge_type(g, et);
      }
    };
    

    /** This class represents an edge list stored on local machine*/
    struct local_edge_list_type {
      make_local_edge_type_functor me_functor;
      typename local_graph_type::edge_list_type elist;
      
      typedef boost::transform_iterator<make_local_edge_type_functor,
                                      typename local_graph_type::edge_list_type::iterator> iterator;
      typedef iterator const_iterator;
      
      local_edge_list_type(distributed_graph& g,
                           typename local_graph_type::edge_list_type elist) :
                          me_functor(g), elist(elist) { }
      
      size_t size() const { return elist.size(); }
      local_edge_type operator[](size_t i) const { return me_functor(elist[i]); }
      iterator begin() const { return
          boost::make_transform_iterator(elist.begin(), me_functor); }
      iterator end() const { return
          boost::make_transform_iterator(elist.end(), me_functor); }
      bool empty() const { return elist.empty(); }
    }; 



    /**
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




    /// The master vertex record map
    // typedef boost::unordered_map<vertex_id_type, vertex_record>  vid2record_type;
    typedef std::vector<vertex_record> lvid2record_type;

  private:
      
    // PRIVATE DATA MEMBERS ===================================================> 
    /** The rpc interface for this class */
    mutable dc_dist_object<distributed_graph> rpc;

    bool finalized;

    /** The local graph data */
    local_graph_type local_graph;
    
    /** The map from global vertex ids to vertex records */
    lvid2record_type lvid2record;
    
    /** The map from global vertex ids back to local vertex ids */
    // boost::unordered_map<vertex_id_type, lvid_type> vid2lvid;
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
    idistributed_ingress_type* ingress_ptr; 

  public:

    // CONSTRUCTORS ==========================================================>
    distributed_graph(distributed_control& dc, 
                      const graphlab_options& opts = graphlab_options() ) : 
      rpc(dc, this), finalized(false), vid2lvid(-1),
      nverts(0), nedges(0), local_own_nverts(0), nreplicas(0),
      ingress_ptr(NULL) {
      rpc.barrier();
      std::string ingress_method = "random";
      opts.get_graph_options().get_option("ingress", ingress_method);

      size_t bufsize = 50000;
      bool usehash = false;
      bool userecent = false;
      opts.get_graph_options().get_option("bufsize", bufsize);
      opts.get_graph_options().get_option("usehash", usehash);
      opts.get_graph_options().get_option("userecent", userecent);
      set_ingress_method(ingress_method, bufsize, usehash, userecent);
    }



    // METHODS ===============================================================>


    void set_ingress_method(const std::string& method, 
        size_t bufsize = 50000, bool usehash = false, bool userecent = false) {

      if(ingress_ptr != NULL) { delete ingress_ptr; ingress_ptr = NULL; }
      if (method == "batch") {
        logstream(LOG_INFO) << "Use batch ingress, bufsize: " << bufsize  
          << ", usehash: " << usehash << ", userecent" << userecent << std::endl;
        ingress_ptr = new distributed_batch_ingress_type(rpc.dc(), *this, 
                                                         bufsize, usehash, userecent);
      } else if (method == "oblivious") {
        logstream(LOG_INFO) << "Use oblivious ingress, usehash: " << usehash 
          << ", userecent: " << userecent << std::endl;
        ingress_ptr = new distributed_oblivious_ingress_type(rpc.dc(), *this, 
                                                             usehash, userecent);
      } else if (method == "identity") {
        logstream(LOG_INFO) << "Use identity ingress" << std::endl;
        ingress_ptr = new distributed_identity_ingress_type(rpc.dc(), *this);
      } else {
        ingress_ptr = new distributed_random_ingress_type(rpc.dc(), *this);
      }
    } // end of set ingress method
    

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
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const { return nverts; }

    /** \brief Get the number of edges */
    size_t num_edges() const { return nedges; }

    /** \brief Get the size of replica */
    size_t num_replicas() const { return nreplicas; }


    /** Determine the id of the vertex with highest degree */
    vertex_id_type max_degree_vertex() const { 
      // First compute the locally maximum vertex
      vertex_id_type max_vid = -1;
      size_t max_degree = 0;
      foreach(const vertex_record& vrec, lvid2record) {
        if(vrec.owner == rpc.procid()) {
          const size_t degree = vrec.num_in_edges + vrec.num_out_edges;
          if(degree > max_degree || max_vid == vertex_id_type(-1)) {
            max_vid = vrec.gvid; max_degree = degree; 
          }
        }
      }
      ASSERT_NE(max_vid, vertex_id_type(-1));
      // Exchange with other machines
      typedef std::pair<size_t, vertex_id_type> pair_type;
      std::vector<pair_type> local_bests(rpc.numprocs());
      local_bests[rpc.procid()] = pair_type(max_degree, max_vid);
      rpc.all_gather(local_bests);
      return std::max_element(local_bests.begin(), local_bests.end())->second;      
    } // end of max_degree_vertex

    /** \brief Get the number of vertices local to this proc */
    size_t num_local_vertices() const { return local_graph.num_vertices(); }

    /** \brief Get the number of edges local to this proc */
    size_t num_local_edges() const { return local_graph.num_edges(); }

    /** \brief Get the number of vertices owned by this proc */
    size_t num_local_own_vertices() const { return local_own_nverts; }

    /** \brief get the local vertex id */
    lvid_type local_vid (const vertex_id_type vid) const {
      // typename boost::unordered_map<vertex_id_type, lvid_type>::
      //   const_iterator iter = vid2lvid.find(vid);
      typename cuckoo_map_type::const_iterator iter = vid2lvid.find(vid);
      ASSERT_TRUE(iter != vid2lvid.end());
      return iter->second;
    } // end of local_vertex_id

    vertex_id_type global_vid(const lvid_type lvid) const { 
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid].gvid;
    } // end of global_vertex_id

    const vertex_record& get_vertex_record(const vertex_id_type vid) const {
      // typename boost::unordered_map<vertex_id_type, lvid_type>::
      //   const_iterator iter = vid2lvid.find(vid);
      typename cuckoo_map_type::const_iterator iter = vid2lvid.find(vid);
      ASSERT_TRUE(iter != vid2lvid.end());
      return lvid2record[iter->second];
    }


    vertex_record& l_get_vertex_record(const lvid_type lvid) {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid];
    }

    const vertex_record& l_get_vertex_record(const lvid_type lvid) const {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid];
    }

    bool is_master(const vertex_id_type vid) const {
      typename cuckoo_map_type::const_iterator iter = vid2lvid.find(vid);
      return (iter != vid2lvid.end()) && l_is_master(iter->second);
    }

    bool l_is_master(const lvid_type lvid) const {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid].owner == rpc.procid();
    }


    local_graph_type& get_local_graph() {
      return local_graph;
    }

    const local_graph_type& get_local_graph() const {
      return local_graph;
    }


    vertex_type vertex(vertex_id_type vid) {
      return vertex_type(*this, local_vid(vid));
    }

    const vertex_type vertex(vertex_id_type vid) const {
      return vertex_type(*this, local_vid(vid));
    }

    local_vertex_type l_vertex(vertex_id_type vid) {
      return local_vertex_type(*this, vid);
    }

    const local_vertex_type l_vertex(vertex_id_type vid) const {
      return local_vertex_type(*this, vid);
    }
  
    //! Get all the edge which edge.target() == v
    edge_list_type in_edges(const vertex_id_type vid) const __attribute__((noreturn)) {
      // Not implemented.
      logstream(LOG_WARNING) << "in_edges not implemented. " << std::endl;
      ASSERT_TRUE(false);
    }


    //! Get the number of edges which edge.target() == v
    size_t num_in_edges(const vertex_id_type vid) const {
      return get_vertex_record(vid).num_in_edges;
    }

    //! Get all the edges which edge.source() == v
    edge_list_type out_edges(const vertex_id_type vid) const __attribute__((noreturn)) {
            // Not implemented.
      logstream(LOG_WARNING) << "in_edges not implemented. " << std::endl;
      ASSERT_TRUE(false);
    }

    //! Get the number of edges which edge.source() == v
    size_t num_out_edges(const vertex_id_type vid) const {
      return get_vertex_record(vid).num_out_edges;
    }

    // Get all the local edge which edge.target() == v
    local_edge_list_type l_in_edges(const lvid_type lvid) {
      return local_edge_list_type(*this, local_graph.in_edges(lvid));
    }

    // Get the number of local edges which edge.target() == v
    size_t l_num_in_edges(const lvid_type lvid) const { 
      return local_graph.num_in_edges(lvid);
    }

    // Get all the local edges which edge.source() == v
    local_edge_list_type l_out_edges(const lvid_type lvid) {
      return local_edge_list_type(*this, local_graph.out_edges(lvid));
    }

    // Get the number of local edges which edge.source() == v
    size_t l_num_out_edges(const lvid_type lvid) const {
      return local_graph.num_out_edges(lvid);
    }

    /** \brief Returns a reference to the data stored on the vertex
        v. */
    VertexData& vertex_data(vertex_id_type vid) {
      return local_graph.vertex_data(local_vid(vid));
    }
    
    /** \brief Returns a constant reference to the data stored on the
        vertex v */
    const VertexData& vertex_data(vertex_id_type vid) const {
      return local_graph.vertex_data(local_vid(vid));
    }

    /** \brief Returns a reference to the data stored on the edge
        source->target. */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target){
      return local_graph.edge_data(local_vid(source), local_vid(target));
    }
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, 
                              vertex_id_type target) const {
      return local_graph.edge_data(local_vid(source), local_vid(target));
    }

   
    /** 
     * \brief Creates a vertex containing the vertex data
     */
    void add_vertex(const vertex_id_type& vid, 
                    const VertexData& vdata = VertexData() ) {
      ASSERT_NE(ingress_ptr, NULL);
      ingress_ptr->add_vertex(vid, vdata);
    }

    /**
     * \brief Creates an edge connecting vertex source to vertex target.  
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const EdgeData& edata = EdgeData()) {
      ASSERT_NE(ingress_ptr, NULL);
      ingress_ptr->add_edge(source, target, edata);
    }



    /**
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
        // Be sure to flush on the last vertex
        if(lvid+1 == lvid2record.size()) vertex_exchange.flush();
        // Receive any vertex data and update local mirrors
        while(vertex_exchange.recv(sending_proc, recv_buffer)) {
          foreach(const pair_type& pair, recv_buffer) 
            vertex_data(pair.first) = pair.second;
          recv_buffer.clear();
        }
      }
      ASSERT_TRUE(vertex_exchange.empty());
    } // end of synchronize



    void resize (size_t n) { }

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
    void load(std::string& path, std::string& prefix) {
      std::ostringstream ss;
      ss << prefix << rpc.procid() << ".bin";
      std::string fname = ss.str();

      if (path.substr(path.length()-1, 1) != "/")
        path.append("/");
      fname = path.append(fname);

      logstream(LOG_INFO) << "Load graph from " << fname << std::endl;
      if(boost::starts_with(fname, "hdfs://")) {
        graphlab::hdfs hdfs;
        graphlab::hdfs::fstream in_file(hdfs, fname);
        boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
        fin.push(in_file);
        if(!fin.good()) {
          logstream(LOG_FATAL) << "Error opening file: " << fname << std::endl;
          exit(-1);
        }
        iarchive iarc(fin);
        iarc >> *this;
        fin.pop();
        in_file.close();
      } else {
        std::ifstream fin(fname.c_str());
        iarchive iarc(fin);
        iarc >> *this;
        fin.close();
      }
      logstream(LOG_INFO) << "Finish loading graph from " << fname << std::endl;
    } // end of load

    /** \brief Load part of the distributed graph from a path*/
    void save(std::string& path, std::string& prefix) {
      timer savetime;  savetime.start();
      std::ostringstream ss;
      ss << prefix << rpc.procid() << ".bin";
      std::string fname = ss.str();
      if (path.substr(path.length()-1, 1) != "/")
        path.append("/");
      fname = path.append(fname);
      logstream(LOG_INFO) << "Save graph to " << fname << std::endl;
      if(boost::starts_with(fname, "hdfs://")) {
        graphlab::hdfs hdfs;
        graphlab::hdfs::fstream out_file(hdfs, fname, true);
        boost::iostreams::filtering_stream<boost::iostreams::output> fout;  
        fout.push(out_file);
        if (!fout.good()) {
          logstream(LOG_FATAL) << "Error opening file: " << fname << std::endl;
          exit(-1);
        }
        oarchive oarc(fout);
        oarc << *this;
        fout.pop();
        out_file.close();
      } else {
        std::ofstream fout(fname.c_str());
        oarchive oarc(fout);
        oarc << *this;
        fout.close();
      }
      rpc.full_barrier();
      logstream(LOG_INFO) << "Finish saving graph to " << fname << std::endl;
      std::cout << "Finished saving binary graph: " 
                << savetime.current_time() << std::endl;
    } // end of save






    // Synthetic Generators ===================================================>
    void build_powerlaw(const size_t nverts, const bool in_degree = false, 
                        const double alpha = 2.1, 
                        const size_t truncate = size_t(-1)) {
      std::vector<double> prob(std::min(nverts, truncate), 0);
      std::cout << "constructing pdf" << std::endl;
      for(size_t i = 0; i < prob.size(); ++i) 
        prob[i] = std::pow(double(i+1), -alpha);
      std::cout << "constructing cdf" << std::endl;
      pdf2cdf(prob);
      std::cout << "Building graph" << std::endl;
      size_t target_index = rpc.procid();
      size_t addedvtx = 0;
      for(size_t source = rpc.procid(); source < nverts; 
          source += rpc.numprocs()) {
        const size_t out_degree = sample(prob) + 1;
        for(size_t i = 0; i < out_degree; ++i, ++target_index) {
          size_t target = permutation(nverts, target_index);
          if(source == target) target = permutation(nverts, ++target_index);
          if(in_degree) add_edge(target, source); 
          else add_edge(source, target);
        }
        ++addedvtx;
        if (addedvtx % 10000000 == 0) {
          std::cout << addedvtx << " inserted\n";
        }
      }
    } // end of build powerlaw


    void build_lognormal(const size_t nverts, const bool in_degree = false,  
                         const double mu = 4, const double sigma = 1.3) {
      random::seed(rpc.procid());
      atomic<size_t> target_index = rpc.procid();
      size_t edges_added = 0;
      //#pragma omp parallel for
      for(size_t source = rpc.procid(); source < nverts; 
          source += rpc.numprocs()) {
        const size_t out_degree = 
          std::min(size_t(std::exp(random::gaussian(mu, sigma))), nverts-1);
        for(size_t i = 0; i < out_degree; ++i, ++target_index) {
          size_t target = permutation(nverts, target_index);
          if(source == target) target = permutation(nverts, ++target_index);
          if(in_degree) add_edge(target, source); 
          else add_edge(source, target);
          edges_added++;
        }
        if(edges_added % 1000000 == 0) {
          std::cout << "Edges: " << edges_added 
                    << "\t Vertices:  " << (source/rpc.numprocs())
                    << "\t Ratio:  "  << double(edges_added) / (source/rpc.numprocs())
                    << std::endl;
        }
      }
    } // end of build lognormal

  private:
      
    inline size_t permutation(size_t nverts, size_t x) const  {
      return ((x + rpc.procid()) * 2654435761) % nverts;
    }

    void pdf2cdf(std::vector<double>& pdf) const {
      double Z = 0;
      for(size_t i = 0; i < pdf.size(); ++i) Z += pdf[i];
      for(size_t i = 0; i < pdf.size(); ++i) 
        pdf[i] = pdf[i]/Z + ((i>0)? pdf[i-1] : 0);
    } // end of pdf2cdf
    
    size_t sample(const std::vector<double>& cdf) const {
      return std::upper_bound(cdf.begin(), cdf.end(), 
                              graphlab::random::rand01()) - cdf.begin();  
    } // end of sample




  }; // End of graph

 
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

