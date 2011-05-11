/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DISK_ATOM_HPP
#define DISK_ATOM_HPP
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
 * "head_vid" ==> uint64_t : the first vertex in a linked list of vertices \n
 * "tail_vid" ==> uint64_t : the last vertex in a linked list of vertices \n
 * "ll[vid]" ==> uint64_t : the vertex after vertex 'vid' in a linked list of vertices \n
 * "numv" ==> uint64_t : the number of vertices in the atom \n
 * "nume" ==> uint64_t : the number of edges in the atom \n
 * "numlocalv" ==> uint64_t : the number of local vertices in the atom. \n
 * "numlocale" ==> uint64_t : the number of local edges in the atom. \n
 * v[vid] ==> archive of owner, vdata : The vertex data of vertex 'vid' \n
 * e[srcv][destv] ==> archive of edata : The edge on the edge srcv --> destv \n
 * i[vid] ==> uint64_t* : An array of in-vertices of vertex 'vid' \n
 * o[vid] ==> uint64_t* : An array of out-vertices of vertex 'vid' \n
 * c[vid] ==> uint32_t : Color of vertex 'vid' \n
 * 
 * Next the DHT entries are: \n
 * h[vid] ==> uint16_t : The atom file owning vertex vid.  \n
 * These are entirely independent of the previous keys.
 */
class disk_atom{
 private:
  kyotocabinet::HashDB db;
  kyotocabinet::CacheDB cache;  // a REALLY simple cache of the db.
                                // with only one global invalidate flag
  //kyotocabinet::HashDB db;
  
  /// a linked list of vertex IDs stored here
  uint64_t head_vid;
  uint64_t tail_vid;
  bool cache_invalid;
  
  atomic<uint64_t> numv;
  atomic<uint64_t> nume;
  atomic<uint64_t> numlocalv;
  atomic<uint64_t> numlocale;
  uint16_t atomid;
  mutex mut;
  
  std::string filename;
  
  inline std::string id_to_str(uint64_t i) {
    char c[10]; 
    unsigned char len = compress_int(i, c);
    return std::string(c + 10 - len, (std::streamsize)len); 
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
  void add_vertex(vertex_id_t vid, uint16_t owner);
  
  /**
   * Reads the entire hash table into cache
   */
  void precache();
  
  /**
   * \brief Inserts vertex 'vid' into the file without data.
   * If the vertex already exists, nothing will be done.
   * Returns true if vertex was added.
   */
  bool add_vertex_skip(vertex_id_t vid, uint16_t owner);
  
  
  /**
   * \brief Inserts vertex 'vid' into the file. If the vertex already exists,
   * it will be overwritten.
   */
  template <typename T>
  void add_vertex(vertex_id_t vid, uint16_t owner, const T &vdata) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner << vdata;
    strm.flush();
    if (db.add("v"+id_to_str(vid), strm.str())) {
      // first entry in linked list
      mut.lock();
      if (head_vid == (uint64_t)(-1)) {
        head_vid = vid;
        tail_vid = vid;
        mut.unlock();
      }
      else {
        // update linked list
        std::string tail_next_key = "ll" + id_to_str(tail_vid);
        tail_vid = vid;
        mut.unlock();
        db.set(tail_next_key.c_str(), tail_next_key.length(),
               (char*)&vid, sizeof(vid));
        cache_invalid = true;
      }
      numv.inc();
      if (owner == atomid) numlocalv.inc();
    }
    else {
      db.set("v"+id_to_str(vid), strm.str());
      cache_invalid = true;
    }
  }
  
  
  /**
   * \brief Inserts edge src->target into the file without data. 
   * If the edge already exists it will be overwritten.
   */
  void add_edge(vertex_id_t src, vertex_id_t target);
  
  /**
   * \brief Inserts edge src->target into the file without data. 
   * If the edge already exists, nothing will be done.
   * Returns true if edge was added.
   */
  bool add_edge_skip(vertex_id_t src, vertex_id_t target);
  
  /**
   * \brief Inserts edge src->target into the file. If the edge already exists,
   * it will be overwritten.
   */
  template <typename T>
  void add_edge(vertex_id_t src, vertex_id_t target, const T &edata) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << edata;
    strm.flush();
    if (db.add("e"+id_to_str(src)+"_"+id_to_str(target), strm.str())) {
      // increment the number of edges
      nume.inc();
      // append to the adjacency entries
      std::string oadj_key = "o"+id_to_str(src);
      uint64_t target64 = (uint64_t)target;
      db.append(oadj_key.c_str(), oadj_key.length(), (char*)&target64, sizeof(target64));
      
      std::string iadj_key = "i"+id_to_str(target);
      uint64_t src64 = (uint64_t)src;
      db.append(iadj_key.c_str(), iadj_key.length(), (char*)&src64, sizeof(src64));
      cache_invalid = true;
    }
    else {
      db.set("e"+id_to_str(src)+"_"+id_to_str(target), strm.str());
      cache_invalid = true;
    }
  }
  
  
  /**
   * \brief Modifies an existing vertex in the file where no data is assigned to the 
   * vertex. User must ensure that the file
   * already contains this vertex. If user is unsure, add_vertex should be used.
   */
  template <typename T>
  void set_vertex(vertex_id_t vid, uint16_t owner) {
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
  void set_vertex(vertex_id_t vid, uint16_t owner, const T &vdata) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner << vdata;
    strm.flush();
    db.set("v"+id_to_str(vid), strm.str());
    cache_invalid = true;
  }
  

  /**
   * \brief Modifies an existing edge in the file where no data is assigned to the edge. 
   * User must ensure that the file already contains this edge. 
   * If user is unsure, add_edge should be used.
   */
  template <typename T>
  void set_edge(vertex_id_t src, vertex_id_t target) {
    db.set("e"+id_to_str(src)+"_"+id_to_str(target), std::string(""));
    cache_invalid = true;
  }
  
  /**
   * \brief Modifies an existing edge in the file. User must ensure that the file
   * already contains this edge. If user is unsure, add_edge should be used.
   */
  template <typename T>
  void set_edge(vertex_id_t src, vertex_id_t target, const T &edata) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << edata;
    strm.flush();
    db.set("e"+id_to_str(src)+"_"+id_to_str(target), strm.str());
    cache_invalid = true;
  }
  
  
 /**
   * \brief Reads a vertex from the file returning only the 'owner' of the vertex
   * and not the data.
   * Returns true if vertex exists and false otherwise.
   */
  bool get_vertex(vertex_id_t vid, uint16_t &owner);
  

  /**
   * \brief Reads a vertex from the file returning results in 'owner' and 'vdata'.
   * Returns true if vertex exists and false otherwise.
   * If there is no vertex data stored, vdata will not be modified.
   */
  template <typename T>
  bool get_vertex(vertex_id_t vid, uint16_t &owner, T &vdata) {
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
  bool get_edge(vertex_id_t src, vertex_id_t target, T &edata) {
    std::string val;
    std::string key = "e"+id_to_str(src)+"_"+id_to_str(target);
    if (cache_invalid || cache.get(key, &val) == false) {
      if (db.get(key, &val) == false) return false;
    }
    if (val.length() > 0) {
      // there is edata
      boost::iostreams::stream<boost::iostreams::array_source> 
                                  istrm(val.c_str(), val.length());   
      iarchive iarc(istrm);
      iarc >> edata;
    }
    return true;
  }


  /**
   * \brief Returns a list of all the vertices in the file
   */
  std::vector<vertex_id_t> enumerate_vertices();
  
  /**
   * \brief Returns a list of all the adjacent atoms in the file
   * and the number of ghost vertices in this atom belonging to the
   * adjacent atom
   */
  std::map<uint16_t, uint32_t> enumerate_adjacent_atoms();
  
  /**
   * \brief Returns the set of incoming vertices of vertex 'vid'
   */
  std::vector<vertex_id_t> get_in_vertices(vertex_id_t vid);
   
   
  /**
   * \brief Returns the set of outgoing vertices of vertex 'vid'
   */
  std::vector<vertex_id_t> get_out_vertices(vertex_id_t vid);


  /**
   * \brief Get the color of the vertex 'vid'.
   * Returns (uint32_t)(-1) if the entry does not exist
   */
  uint32_t get_color(vertex_id_t vid);

  /**
   * \brief Sets the color of vertex 'vid'
   */
  void set_color(vertex_id_t vid, uint32_t color);
  
  /// Returns the largest color number
  uint32_t max_color();
  
  /**
   * \brief Reads from the auxiliary hash table mapping vid ==> owner.
   * Returns (uint16_t)(-1) if the entry does not exist
   */
  uint16_t get_owner(vertex_id_t vid);

  /**
   * \brief Writes to the auxiliary hash table mapping vid ==> owner.
   */
  void set_owner(vertex_id_t vid, uint16_t owner);

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
  inline kyotocabinet::HashDB& get_db() {
    return db;
  }
};

}

#endif
