#ifndef DISK_ATOM_HPP
#define DISK_ATOM_HPP
#include <sstream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/logger.hpp>
#include <kchashdb.h>

namespace graphlab {
  
/**
 * Interface for reading and writing to an atom on disk
 * 
 * The atom serves 2 purposes.
 * First it provides a partition of the graph, including its ghost vertices.
 * Next, it provides a hashtable distributed across all atom files, identifying 
 *       the owner of a vertex.
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
  //kyotocabinet::HashDB db;
  
  /// a linked list of vertex IDs stored here
  uint64_t head_vid;
  uint64_t tail_vid;
  atomic<uint64_t> numv;
  atomic<uint64_t> nume;
  atomic<uint64_t> numlocalv;
  atomic<uint64_t> numlocale;
  uint16_t atomid;
  mutex mut;
  
  inline std::string id_to_str(uint64_t i) {
    char c[10]; 
    unsigned char len = compress_int(i, c);
    return std::string(c + 10 - len, (std::streamsize)len); 
  }

 public:
   
  /// constructor. Accesses an atom stored at the filename provided
  inline disk_atom(std::string filename, 
                   uint16_t atomid):atomid(atomid) {
    db.tune_buckets(1000 * 1000);
#if __LP64__
    db.tune_map(256 * 1024 * 1024); // 256MB
#endif
    ASSERT_TRUE(db.open(filename));
    // get the pointers to the linked list of vertices
    if (db.get("head_vid", 7, (char*)&head_vid, sizeof(head_vid)) == -1) head_vid = (uint64_t)(-1);
    if (db.get("tail_vid", 7, (char*)&tail_vid, sizeof(tail_vid)) == -1) tail_vid = (uint64_t)(-1);
    if (db.get("numv", 4, (char*)&numv.value, sizeof(numv.value)) == -1) numv = 0;
    if (db.get("nume", 4, (char*)&nume.value, sizeof(nume.value)) == -1) nume = 0;
    if (db.get("numlocalv", 9, (char*)&numlocalv.value, sizeof(numlocalv.value)) == -1) numlocalv = 0;
    if (db.get("numlocale", 9, (char*)&numlocale.value, sizeof(numlocale.value)) == -1) numlocale = 0;
  }

  
  inline uint16_t atom_id() {
    return atomid;
  }
  
  inline ~disk_atom() {
    // before we close the file, we update the linked list
    db.set("head_vid", 7, (char*)&head_vid, sizeof(head_vid));
    db.set("tail_vid", 7, (char*)&tail_vid, sizeof(tail_vid));
    db.set("numv", 4, (char*)&numv.value, sizeof(numv.value));
    db.set("nume", 4, (char*)&nume.value, sizeof(nume.value));
    db.set("numlocalv", 9, (char*)&numlocalv.value, sizeof(numlocalv.value));
    db.set("numlocale", 9, (char*)&numlocale.value, sizeof(numlocale.value));
    db.synchronize();
    db.close();
  }
  
  
  void inc_numlocale() {
    numlocale.inc();
  }
  
  /**
   * Inserts vertex 'vid' into the file without data.
   * If the vertex already exists, it will be overwritten.
   */
  void add_vertex(vertex_id_t vid, uint16_t owner) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner;
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
      }
      numv.inc();
      if (owner == atomid) numlocalv.inc();
    }
    else {
      db.set("v"+id_to_str(vid), strm.str());
    }
  }
  
  /**
   * Inserts vertex 'vid' into the file without data.
   * If the vertex already exists, nothing will be done.
   * Returns true if vertex was added.
   */
  bool add_vertex_skip(vertex_id_t vid, uint16_t owner) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner;
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
      }
      numv.inc();
      if (owner == atomid) numlocalv.inc();
      return true;
    }
    return false;
  }
  
  
  /**
   * Inserts vertex 'vid' into the file. If the vertex already exists,
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
      }
      numv.inc();
      if (owner == atomid) numlocalv.inc();
    }
    else {
      db.set("v"+id_to_str(vid), strm.str());
    }
  }
  
  
  /**
   * Inserts edge src->target into the file without data. 
   * If the edge already exists it will be overwritten.
   */
  void add_edge(vertex_id_t src, vertex_id_t target) {
    if (db.add("e"+id_to_str(src)+"_"+id_to_str(target), std::string(""))) {
      // increment the number of edges
      nume.inc();
      // append to the adjacency entries
      std::string oadj_key = "o"+id_to_str(src);
      uint64_t target64 = (uint64_t)target;
      db.append(oadj_key.c_str(), oadj_key.length(), (char*)&target64, sizeof(target64));
      
      std::string iadj_key = "i"+id_to_str(target);
      uint64_t src64 = (uint64_t)src;
      db.append(iadj_key.c_str(), iadj_key.length(), (char*)&src64, sizeof(src64));
    }
    else {
      db.set("e"+id_to_str(src)+"_"+id_to_str(target), std::string(""));
    }
  }
  
  /**
   * Inserts edge src->target into the file. If the edge already exists,
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
    }
    else {
      db.set("e"+id_to_str(src)+"_"+id_to_str(target), strm.str());
    }
  }
  
  
  /**
   * Modifies an existing vertex in the file where no data is assigned to the 
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
  }
  

  /**
   * Modifies an existing vertex in the file. User must ensure that the file
   * already contains this vertex. If user is unsure, add_vertex should be used.
   */
  template <typename T>
  void set_vertex(vertex_id_t vid, uint16_t owner, const T &vdata) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << owner << vdata;
    strm.flush();
    db.set("v"+id_to_str(vid), strm.str());
  }
  

  /**
   * Modifies an existing edge in the file where no data is assigned to the edge. 
   * User must ensure that the file already contains this edge. 
   * If user is unsure, add_edge should be used.
   */
  template <typename T>
  void set_edge(vertex_id_t src, vertex_id_t target) {
    db.set("e"+id_to_str(src)+"_"+id_to_str(target), std::string(""));
  }
  
  /**
   * Modifies an existing edge in the file. User must ensure that the file
   * already contains this edge. If user is unsure, add_edge should be used.
   */
  template <typename T>
  void set_edge(vertex_id_t src, vertex_id_t target, const T &edata) {
    std::stringstream strm;
    oarchive oarc(strm);    
    oarc << edata;
    strm.flush();
    db.set("e"+id_to_str(src)+"_"+id_to_str(target), strm.str());
  }
  
  
 /**
   * Reads a vertex from the file returning only the 'owner' of the vertex
   * and not the data.
   * Returns true if vertex exists and false otherwise.
   */
  template <typename T>
  bool get_vertex(vertex_id_t vid, uint16_t &owner) {
    std::string val;
    if (db.get("v"+id_to_str(vid), &val) == false) return false;
    
    boost::iostreams::stream<boost::iostreams::array_source> 
                                istrm(val.c_str(), val.length());   
    iarchive iarc(istrm);
    iarc >> owner;
    return true;
  }
  

  /**
   * Reads a vertex from the file returning results in 'owner' and 'vdata'.
   * Returns true if vertex exists and false otherwise.
   * If there is no vertex data stored, vdata will not be modified.
   */
  template <typename T>
  bool get_vertex(vertex_id_t vid, uint16_t &owner, T &vdata) {
    std::string val;
    if (db.get("v"+id_to_str(vid), &val) == false) return false;
    
    boost::iostreams::stream<boost::iostreams::array_source> 
                                istrm(val.c_str(), val.length());   
    iarchive iarc(istrm);
    iarc >> owner;
    // try to deserialize vdata if exists
    if (!istrm.eof()) iarc >> vdata;
    return true;
  }
  
  /**
   * Reads a edge from the file returning results in 'owner' and 'vdata'.
   * Returns true if edge exists and false otherwise.
   * If there is no edge data stored, edata will not be modified.
   */
  template <typename T>
  bool get_edge(vertex_id_t src, vertex_id_t target, T &edata) {
    std::string val;
    if (db.get("e"+id_to_str(src)+"_"+id_to_str(target), &val) == false) return false;
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
   * Returns a list of all the vertices in the file
   */
  inline std::vector<vertex_id_t> enumerate_vertices() {
    std::vector<vertex_id_t> ret;
    if (head_vid == (uint64_t)(-1)) return ret;
    else {
      uint64_t curvid = head_vid;
      while(1) {
        ret.push_back((vertex_id_t)(curvid));
        std::string next_key = "ll" + id_to_str(curvid);
        if (db.get(next_key.c_str(), next_key.length(), (char*)&curvid, sizeof(curvid)) == -1) {
          break;
        }
      }
    }
    return ret;
  }
  
  /**
   * Returns the set of incoming vertices of vertex 'vid'
   */
  inline std::vector<vertex_id_t> get_in_vertices(vertex_id_t vid) {
    std::vector<vertex_id_t> ret;
    std::string val;
    if (db.get("i"+id_to_str(vid), &val)) {
      const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
      ASSERT_TRUE(val.length() % 8 == 0);
      size_t numel = val.length() / 8;
      ret.resize(numel);
      for (size_t i = 0;i < numel; ++i) ret[i] = v[i];
    }
    return ret;    
  }
   
   
  /**
   * Returns the set of outgoing vertices of vertex 'vid'
   */
  inline std::vector<vertex_id_t> get_out_vertices(vertex_id_t vid) {
    std::vector<vertex_id_t> ret;
    std::string val;
    if (db.get("o"+id_to_str(vid), &val)) {
      const uint64_t* v = reinterpret_cast<const uint64_t*>(val.c_str());
      size_t numel = val.length() / 8;
      ASSERT_TRUE(val.length() % 8 == 0);
      ret.resize(numel);
      for (size_t i = 0;i < numel; ++i) ret[i] = v[i];
    }
    return ret;    
  }


  /**
   * Get the color of the vertex 'vid'.
   * Returns (uint32_t)(-1) if the entry does not exist
   */
  uint32_t get_color(vertex_id_t vid) {
    std::string key = "c" + id_to_str(vid);
    uint32_t ret;
    if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1) ret = (uint32_t)(-1);
    return ret;
  }

  /**
   * Sets the color of vertex 'vid'
   */
  void set_color(vertex_id_t vid, uint32_t color) {
    std::string key = "c" + id_to_str(vid);
    db.set(key.c_str(), key.length(), (char*)&color, sizeof(color));
  }
  
  /**
   * Reads from the hash table mapping vid ==> owner.
   * Returns (uint16_t)(-1) if the entry does not exist
   */
  uint16_t get_owner(vertex_id_t vid) {
    std::string key = "h" + id_to_str(vid);
    uint16_t ret;
    if (db.get(key.c_str(), key.length(), (char*)&ret, sizeof(ret)) == -1) ret = (uint16_t)(-1);
    return ret;
  }

  /**
   * Writes to the hash table mapping vid ==> owner.
   */
  void set_owner(vertex_id_t vid, uint16_t owner) {
    std::string key = "h" + id_to_str(vid);
    db.set(key.c_str(), key.length(), (char*)&owner, sizeof(owner));
  }

  /// empties the atom file
  inline void clear() {
    head_vid = (uint64_t)(-1);   
    tail_vid = (uint64_t)(-1);
    numv.value = 0;
    nume.value = 0;
    numlocalv.value = 0;
    numlocale.value = 0;
    db.clear();
  }
  
  inline uint64_t num_vertices() const {
    return numv.value;
  }
  
  inline uint64_t num_edges() const {
    return nume.value;
  }
  
  inline uint64_t num_local_vertices() const {
    return numlocalv.value;
  }
  
  inline uint64_t num_local_edges() const {
    return numlocale.value;
  }
};

}

#endif