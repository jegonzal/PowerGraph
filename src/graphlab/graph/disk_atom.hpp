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

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


#ifndef GRAPHLAB_DISK_ATOM_HPP
#define GRAPHLAB_DISK_ATOM_HPP


#include <sstream>
#include <map>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <boost/iostreams/stream.hpp>
#include <graphlab/logger/logger.hpp>
#include <kchashdb.h>
#include <kccachedb.h>

namespace graphlab {
  
  /**
   * Interface for reading and writing to an atom on disk
   * 
   * The atom serves 2 purposes.
   * First it provides a partition of the graph, including its ghost vertices.
   * Next, it provides an auxiliary hashtable distributed across all atom files, 
   *        identifying the owner of a vertex.
   * 
   * The atom file is a Kyoto Cabinet data store and it contains the following keys:
   * 
   * "vidlist" ==> uint64_t : the vertex after vertex 'vid' in a linked list of vertices \n
   * "numv" ==> uint64_t : the number of vertices in the atom \n
   * "nume" ==> uint64_t : the number of edges in the atom \n
   * "numlocalv" ==> uint64_t : the number of local vertices in the atom. \n
   * "numlocale" ==> uint64_t : the number of local edges in the atom. \n
   * v[vid] ==> archive of owner, vdata : The vertex data of vertex 'vid' \n
   * e[hash] ==> archive of (uint64_t, uint64_t, edata)*: All edges on the edge srcv-->Destv hashing to this entry \n
   * i[vid] ==> uint64_t* : An array of in-vertices of vertex 'vid' \n
   * o[vid] ==> uint64_t* : An array of out-vertices of vertex 'vid' \n
   * c[vid] ==> uint32_t : Color of vertex 'vid' \n
   * 
   * Next the DHT entries are: \n
   * h[vid] ==> uint16_t : The atom file owning vertex vid.  \n
   * These are entirely independent of the previous keys.
   */
  class disk_atom {
  public:
    typedef kyotocabinet::TreeDB storage_type;
  private:

    //! Todo: Fix ugly hack
    typedef graph<bool,bool>::vertex_id_type    vertex_id_type;
    typedef graph<bool,bool>::vertex_color_type vertex_color_type;

    storage_type db;
    kyotocabinet::CacheDB cache;  // a REALLY simple cache of the db.
    // with only one global invalidate flag
    //kyotocabinet::HashDB db;
  
    bool cache_invalid;
  
    atomic<uint64_t> numv;
    atomic<uint64_t> nume;
    atomic<uint64_t> numlocalv;
    atomic<uint64_t> numlocale;
    uint16_t atomid;
    mutex mut[511];
  
    std::string filename;
  
    inline std::string id_to_str(uint64_t i) {
      char c[10]; 
      unsigned char len = compress_int(i, c);
      return std::string(c + 10 - len, (std::streamsize)len); 
    }

    uint32_t rotate_right(uint32_t a, uint32_t shift) {
      return a >> shift | a << (sizeof(a) * 8 - shift);
    }

    // compress into 128K buckets
    inline uint32_t edge_hash_32(uint32_t srcv, uint32_t destv) {
       // note that this is a 32 bit hash!
       srcv ^= rotate_right(srcv, 20) ^ rotate_right(srcv, 12) ^
               rotate_right(srcv, 7) ^ rotate_right(srcv, 4);

       destv ^= rotate_right(destv, 20) ^ rotate_right(destv, 12) ^
               rotate_right(destv, 7) ^ rotate_right(destv, 4);

       return (srcv ^ destv) % 131071;
    }
    
    inline uint32_t edge_hash(uint64_t srcv, uint64_t destv) {
      return edge_hash_32((uint32_t)srcv, (uint32_t)destv);
    }

  public:
   
    /// constructor. Accesses an atom stored at the filename provided
    disk_atom(std::string filename, uint16_t atomid);

  
    ~disk_atom();
  
    /// Increments the number of local edges stored in this atom
    inline void inc_numlocale() {
      numlocale.inc();
    }

    /// Gets the atom ID of this atom
    inline uint16_t atom_id() const {
      return atomid;
    }
  
    inline std::string get_filename() const {
      return filename;
    }
  
    /**
     * \brief Inserts vertex 'vid' into the file without data.
     * If the vertex already exists, it will be overwritten.
     */
    void add_vertex(vertex_id_type vid, uint16_t owner);
  
    /**
     * Reads the entire hash table into cache
     */
    void precache();
  
    /**
     * \brief Inserts vertex 'vid' into the file without data.
     * If the vertex already exists, nothing will be done.
     * Returns true if vertex was added.
     */
    bool add_vertex_skip(vertex_id_type vid, uint16_t owner);
  
  
    /**
     * \brief Inserts vertex 'vid' into the file.
     */
    template <typename T>
    bool add_vertex(vertex_id_type vid, uint16_t owner, 
                    const T &vdata, bool hasdata = true, bool overwrite_if_exists = true) {
      std::stringstream strm;
      oarchive oarc(strm);    
      oarc << owner;
      if (hasdata) oarc << vdata;
      strm.flush();
      if (db.add("v"+id_to_str(vid), strm.str())) {
        uint64_t v64 = (uint64_t)vid;
        db.append("vidlist", 7, (char*)&v64, sizeof(v64));
        cache_invalid = true;
        numv.inc();
        if (owner == atomid) numlocalv.inc();
        return true;
      }
      else if (overwrite_if_exists){
        db.set("v"+id_to_str(vid), strm.str());
        cache_invalid = true;
        return true;
      }
      return false;
    }
  
  
    /**
     * \brief Inserts edge src->target into the file without data. 
     * If the edge already exists it will be overwritten.
     */
    void add_edge(vertex_id_type src, vertex_id_type target);
  
    /**
     * \brief Inserts edge src->target into the file without data. 
     * If the edge already exists, nothing will be done.
     * Returns true if edge was added.
     */
    bool add_edge_skip(vertex_id_type src, vertex_id_type target);
  
    /**
     * \brief Inserts edge src->target into the file. If the edge already exists,
     * it will be overwritten.
     */
    template <typename T>
    bool add_edge(vertex_id_type src, vertex_id_type target, const T &edata,
                  bool hasdata = true, bool overwrite_if_exists = true) {
      bool newedge = true;
      std::string existing_bucket_data;
      std::string key = "e"+id_to_str(edge_hash(src, target));
      size_t mutexid = edge_hash(src, target) % 511;
      mut[mutexid].lock();
      
      std::string edatastring;      
      if (hasdata) {
        std::stringstream edatastream;
        oarchive oarc(edatastream);
        oarc << edata;
        edatastream.flush();
        edatastring = edatastream.str();
      } 
      
      // see if the bucket exists      
      size_t cur_edatalen_offset = 0; // the offset of the edatalen if data already exists
      if (db.get(key, &existing_bucket_data)) {
        // bucket exists. See if the key is in it
        boost::iostreams::stream<boost::iostreams::array_source> istrm(existing_bucket_data.c_str(), 
                                                                 existing_bucket_data.length());   
        iarchive iarc(istrm);
        while(istrm.good()) {
          vertex_id_type src_, target_;
          iarc >> src_ >> target_;
          size_t elen;
          if (src_ == src && target_ == target) {
            newedge = false;
            cur_edatalen_offset = istrm.tellg();
            break;
          }
          iarc >> elen;
          istrm.ignore(elen);
        }
      }
      
      if (newedge) {
        // increment the number of edges
        nume.inc();
        // append to the adjacency entries
        std::string oadj_key = "o"+id_to_str(src);
        uint64_t target64 = (uint64_t)target;
        db.append(oadj_key.c_str(), oadj_key.length(), (char*)&target64, sizeof(target64));
      
        std::string iadj_key = "i"+id_to_str(target);
        uint64_t src64 = (uint64_t)src;
        db.append(iadj_key.c_str(), iadj_key.length(), (char*)&src64, sizeof(src64));
        
        std::stringstream strm;
        oarchive oarc(strm);
        oarc << src << target << edatastring.length();
        strm.flush();
        
        db.append(key, strm.str());
        db.append(key, edatastring);
        cache_invalid = true;
        mut[mutexid].unlock();
        return true;
      }      
      else if (overwrite_if_exists) {
        boost::iostreams::stream<boost::iostreams::array_source> istrm(existing_bucket_data.c_str() + cur_edatalen_offset, 
                                                                        existing_bucket_data.length());   
        iarchive iarc(istrm);
        size_t edatalen;
        iarc >> edatalen;
        size_t length_to_exclude = (size_t)(istrm.tellg()) + edatalen; 
        // get the stored version of the data length
        std::stringstream strm;
        oarchive oarc(strm);
        oarc << edatastring.length();
        strm.flush();
        
        existing_bucket_data.replace(cur_edatalen_offset, length_to_exclude, strm.str() + edatastring);
        db.set(key, existing_bucket_data);
        cache_invalid = true;
        mut[mutexid].unlock();
        return true;
      }
      mut[mutexid].unlock();
      return false;
    }
  
  
    /**
     * \brief Modifies an existing vertex in the file where no data is assigned to the 
     * vertex. User must ensure that the file
     * already contains this vertex. If user is unsure, add_vertex should be used.
     */
    template <typename T>
    void set_vertex(vertex_id_type vid, uint16_t owner) {
      std::stringstream strm;
      oarchive oarc(strm);    
      oarc << owner;
      strm.flush();
      db.set("v"+id_to_str(vid), strm.str());
      cache_invalid = true;
    }
  


    /**
     * \brief Modifies an existing vertex in the file. User must ensure that the file
     * already contains this vertex. If user is unsure, add_vertex should be used.
     */
    template <typename T>
    void set_vertex(vertex_id_type vid, uint16_t owner, const T &vdata) {
      std::stringstream strm;
      oarchive oarc(strm);    
      oarc << owner << vdata;
      strm.flush();
      db.set("v"+id_to_str(vid), strm.str());
      cache_invalid = true;
    }
  
  
    /**
     * \brief Modifies an existing edge in the file. User must ensure that the file
     * already contains this edge. If user is unsure, add_edge should be used.
     */
    template <typename T>
    void set_edge(vertex_id_type src, vertex_id_type target, const T &edata) {
      add_edge(src, target, edata, true, true);
    }
  
  
    /**
     * \brief Reads a vertex from the file returning only the 'owner' of the vertex
     * and not the data.
     * Returns true if vertex exists and false otherwise.
     */
    bool get_vertex(vertex_id_type vid, uint16_t &owner);
  

    /**
     * \brief Reads a vertex from the file returning results in 'owner' and 'vdata'.
     * Returns true if vertex exists and false otherwise.
     * If there is no vertex data stored, vdata will not be modified.
     */
    template <typename T>
    bool get_vertex(vertex_id_type vid, uint16_t &owner, T &vdata) {
      std::string val;
      std::string key = "v"+id_to_str(vid);
      if (cache_invalid || cache.get(key, &val) == false) {
        if (db.get("v"+id_to_str(vid), &val) == false) return false;
      }
    
      boost::iostreams::stream<boost::iostreams::array_source> 
        istrm(val.c_str(), val.length());   
      iarchive iarc(istrm);
      iarc >> owner;
      // try to deserialize vdata if exists
      istrm.peek();
      if (!istrm.eof()) iarc >> vdata;
      return true;
    }
  
    /**
     * \brief Reads a edge from the file returning results in 'owner' and 'vdata'.
     * Returns true if edge exists and false otherwise.
     * If there is no edge data stored, edata will not be modified.
     */
    template <typename T>
    bool get_edge(vertex_id_type src, vertex_id_type target, T &edata) {
      std::string existing_bucket_data;
      std::string key = "e"+id_to_str(edge_hash(src, target));

      // see if the bucket exists      
      if (db.get(key, &existing_bucket_data)) {
        // bucket exists. See if the key is in it
        boost::iostreams::stream<boost::iostreams::array_source> istrm(existing_bucket_data.c_str(), 
                                                                 existing_bucket_data.length());   
        iarchive iarc(istrm);
        while(istrm.good()) {
          vertex_id_type src_, target_;
          size_t elen;
          iarc >> src_ >> target_ >> elen;
          if (src_ == src && target_ == target) {
            if (elen == 0) edata = T();
            iarc >> edata;
            return true;
          }
          else {
            istrm.ignore(elen);
          }
        }
      }
      return false;      
    }


    /**
     * \brief Returns a list of all the vertices in the file
     */
    std::vector<vertex_id_type> enumerate_vertices();
  
    /**
     * \brief Returns a list of all the adjacent atoms in the file
     * and the number of ghost vertices in this atom belonging to the
     * adjacent atom
     */
    std::map<uint16_t, uint32_t> enumerate_adjacent_atoms();
  
    /**
     * \brief Returns the set of incoming vertices of vertex 'vid'
     */
    std::vector<vertex_id_type> get_in_vertices(vertex_id_type vid);
   
   
    /**
     * \brief Returns the set of outgoing vertices of vertex 'vid'
     */
    std::vector<vertex_id_type> get_out_vertices(vertex_id_type vid);


    /**
     * \brief Get the color of the vertex 'vid'.
     * Returns vertex_color_type(-1) if the entry does not exist
     */
    vertex_color_type get_color(vertex_id_type vid);

    /**
     * \brief Sets the color of vertex 'vid'
     */
    void set_color(vertex_id_type vid, vertex_color_type color);
  
    /// Returns the largest color number
    vertex_color_type max_color();
  
    /**
     * \brief Reads from the auxiliary hash table mapping vid ==> owner.
     * Returns (uint16_t)(-1) if the entry does not exist
     */
    uint16_t get_owner(vertex_id_type vid);

    /**
     * \brief Writes to the auxiliary hash table mapping vid ==> owner.
     */
    void set_owner(vertex_id_type vid, uint16_t owner);

    /// \brief empties the atom file
    void clear();
  
    /// \brief Ensures the disk storage is up to date
    void synchronize();
    
    /** \brief Return the total number of vertices stored in this atom, 
     * whether or not the this atom actually owns the vertex.
     */
    inline uint64_t num_vertices() const {
      return numv.value;
    }
  
    /** \brief  Return the total number of edges stored in this atom, whether or 
     * not the this atom actually owns the edge.
     */
    inline uint64_t num_edges() const {
      return nume.value;
    }
  
    /// Number of vertices owned by this atom
    inline uint64_t num_local_vertices() const {
      return numlocalv.value;
    }
  
    /// Number of edges owned by this atom
    inline uint64_t num_local_edges() const {
      return numlocale.value;
    }
  
    /// Returns a reference to the underlying Kyoto Cabinet
    inline storage_type& get_db() {
      return db;
    }
  };

}

#endif

